#include "migec/assemble.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>

#include "migec/fastq.hpp"
#include "migec/mig_record.hpp"
#include "migec/resource.hpp"
#include "migec/types.hpp"
#include "migec/umi_stats.hpp"

namespace migec {

namespace {

// The value of a SAM-style tag in a FASTQ comment, or an empty view. Tags are TAB separated and
// the comment itself is whatever followed the first space, so a plain split on TAB is enough.
std::string_view tag_value(std::string_view comment, std::string_view key) {
    size_t pos = 0;
    while (pos <= comment.size()) {
        const size_t end = std::min(comment.find('\t', pos), comment.size());
        const std::string_view field = comment.substr(pos, end - pos);
        if (field.size() > key.size() && field.compare(0, key.size(), key) == 0) {
            return field.substr(key.size());
        }
        pos = end + 1;
    }
    return {};
}

int choose_bucket_bits(const std::string& path) {
    std::error_code ec;
    const uintmax_t on_disk = std::filesystem::file_size(path, ec);
    if (ec) return 0;
    // A gzipped FASTQ expands ~3.5x, and a resident record costs about as much again in the
    // string copies and the record vector. 8x on-disk is the working estimate; it decides only
    // how finely the input is cut, and being one bucket out costs nothing but a little IO.
    uint64_t resident = static_cast<uint64_t>(on_disk) * 8;
    int bits = 0;
    while (bits < kMaxBucketBits && (resident >> bits) > kBucketBudgetBytes) ++bits;
    return bits;
}

struct Bucket {
    std::unique_ptr<MigWriter> writer;
    std::string path;
    uint64_t records = 0;
};

// A record held in RAM while its bucket is being sorted. The reader's views die on the next
// next(), so the payload is copied once, here, and never again.
struct Resident {
    uint64_t cell;
    uint64_t umi;
    uint64_t src_index;
    std::string seq, qual;

    // Cell first, so a cell's molecules are contiguous. The identity of a molecule is the WHOLE
    // key: the same UMI in two cells, or in two samples, is two molecules and always was -- a
    // UMI is only ever unique within the compartment it was added to.
    friend bool operator<(const Resident& a, const Resident& b) {
        if (a.cell != b.cell) return a.cell < b.cell;
        if (a.umi != b.umi) return a.umi < b.umi;
        return a.src_index < b.src_index;
    }
};

void bin(std::vector<uint64_t>& histogram, uint64_t n) {
    size_t b = 0;
    while ((n >> b) > 1) ++b;
    if (histogram.size() <= b) histogram.resize(b + 1, 0);
    ++histogram[b];
}

}  // namespace

AssembleStats assemble(const AssembleRequest& request) {
    Stopwatch clock;
    AssembleStats stats;
    stats.sample_id = request.sample_id;

    std::filesystem::path out_dir(request.output_dir);
    std::filesystem::create_directories(out_dir);
    const std::filesystem::path temp_dir = out_dir / ".assemble_buckets";
    std::filesystem::create_directories(temp_dir);

    const int bits = request.bucket_bits > 0 ? request.bucket_bits
                                             : choose_bucket_bits(request.input);
    const size_t n_buckets = static_cast<size_t>(1) << bits;
    stats.buckets = static_cast<int>(n_buckets);

    // ------------------------------------------------------------------ pass 1: partition
    {
        Stopwatch partition_clock;
        std::vector<Bucket> buckets(n_buckets);
        const size_t block_bytes = std::clamp(
            static_cast<size_t>(kWriterBudgetBytes / n_buckets), kMinBlockBytes, kMaxBlockBytes);
        FastqReader reader(request.input);
        FastqRecord rec;
        uint64_t index = 0;
        while (reader.next(rec)) {
            ++stats.reads;
            const std::string_view umi = tag_value(rec.comment, "RX:Z:");
            if (umi.empty()) { ++stats.reads_without_umi; continue; }
            if (stats.umi_length == 0) {
                stats.umi_length = static_cast<int>(umi.size());
                if (stats.umi_length > kMaxBarcodeLen) {
                    throw MigecError("assemble: a " + std::to_string(stats.umi_length) +
                                     " nt UMI does not fit the packed key (max " +
                                     std::to_string(kMaxBarcodeLen) + ")");
                }
            } else if (static_cast<int>(umi.size()) != stats.umi_length) {
                throw MigecError("assemble: read '" + std::string(rec.name) + "' carries a " +
                                 std::to_string(umi.size()) + " nt UMI where the file started " +
                                 "with " + std::to_string(stats.umi_length) +
                                 " -- these are two runs concatenated, not one sample");
            }
            if (stats.sample_id.empty()) {
                const std::string_view bc = tag_value(rec.comment, "BC:Z:");
                stats.sample_id = bc.empty() ? std::string("sample") : std::string(bc);
            }
            const std::string_view cell = tag_value(rec.comment, "CB:Z:");
            if (stats.cell_length == 0 && !cell.empty()) {
                stats.cell_length = static_cast<int>(cell.size());
            }
            bool has_n = false, cell_has_n = false;
            const uint64_t key = pack_barcode(umi, &has_n);
            const uint64_t cell_key = cell.empty() ? 0 : pack_barcode(cell, &cell_has_n);
            // Partition on the cell when there is one: every read of a cell then lands in one
            // bucket, which is what makes a per-cell scope local. Without cells the UMI is the
            // whole key and the partition is on it.
            const uint64_t partition_key = cell.empty() ? key : cell_key;
            Bucket& b = buckets[bucket_of(partition_key, bits)];
            if (!b.writer) {
                MigHeader header;
                header.umi_len = static_cast<uint8_t>(stats.umi_length);
                header.cell_len = static_cast<uint8_t>(stats.cell_length);
                header.bucket_index = static_cast<uint8_t>(bucket_of(partition_key, bits));
                header.bucket_bits = static_cast<uint8_t>(bits);
                header.paired = false;
                header.sample_id = stats.sample_id;
                b.path = (temp_dir /
                          ("bucket." + std::to_string(bucket_of(partition_key, bits)) + ".mig"))
                             .string();
                b.writer = std::make_unique<MigWriter>(b.path, header, block_bytes);
            }
            MigRecord out;
            out.cell = cell_key;
            out.umi = key;
            out.src_index = index++;
            out.flags = static_cast<uint16_t>(kSingleEnd | (has_n ? kUmiHasN : 0) |
                                              (cell_has_n ? kCellHasN : 0));
            out.seq1 = rec.seq;
            out.qual1 = rec.qual;
            b.writer->write(out);
            ++b.records;
        }
        for (Bucket& b : buckets) {
            if (b.writer) b.writer->close();
        }
        stats.partition_seconds = partition_clock.seconds();
        // Bucket paths are regenerated below from the index, so the writers can go now.
    }

    // ------------------------------------------------------------------ pass 2: consense
    const std::filesystem::path fastq_path =
        out_dir / (stats.sample_id + ".consensus.fq" + (request.gzip_level > 0 ? ".gz" : ""));
    FastqWriter writer(fastq_path.string(), request.gzip_level);
    std::FILE* table = std::fopen((out_dir / (stats.sample_id + ".mig.tsv")).string().c_str(), "w");
    if (!table) throw MigecError("assemble: cannot write the per-molecule table");
    std::fprintf(table,
                 "cell\tumi\tcontig\tcontigs\tmolecule\treads\tsupport\tlength\tmean_quality\t"
                 "consensus_error\tlinkage\n");

    double qual_sum = 0.0, err_sum = 0.0;
    uint64_t qual_bases = 0;
    // Counting mode: reads that carried the emitted sequence, over reads in the groups that
    // emitted one. A ratio of two counts rather than a mean of ratios, so a group of 100 reads
    // does not weigh the same as a group of 2.
    uint64_t support_reads = 0, support_of_reads = 0;
    // Per-position base usage over DISTINCT barcodes, which is what the birthday arithmetic needs
    // -- weighting by reads would let one over-amplified molecule set the composition.
    std::vector<std::array<uint64_t, 4>> usage;

    auto emit = [&](uint64_t cell_key, uint64_t key, const std::vector<Resident>& group) {
        ++stats.groups;
        {
            const std::string bases = unpack_barcode(key, stats.umi_length);
            if (usage.empty()) usage.assign(bases.size(), {});
            for (size_t j = 0; j < bases.size(); ++j) {
                const uint8_t c = base_code(bases[j]);
                if (c != kInvalidBase) ++usage[j][c];
            }
        }
        bin(stats.size_histogram, group.size());
        if (group.size() < request.min_reads) {
            stats.reads_dropped += group.size();
            return;
        }
        // The consensus sees at most kMaxReadsPerGroup of them; the molecule keeps all of them as
        // its depth. Records arrive in src_index order, which is input order and carries no
        // relation to the sequence, so the first N are an unbiased sample of the group and cost
        // nothing to take.
        const size_t depth = group.size();
        const size_t used = std::min(depth, kMaxReadsPerGroup);
        if (used < depth) {
            ++stats.groups_capped;
            stats.reads_over_cap += depth - used;
        }
        std::vector<ConsensusRead> reads;
        reads.reserve(used);
        for (size_t i = 0; i < used; ++i) reads.push_back({group[i].seq, group[i].qual});
        const std::vector<Consensus> molecules = assemble_group(reads, request.consensus);
        const uint32_t components = molecules.empty() ? 1 : molecules[0].components;
        if (components > 1) { ++stats.groups_fragmented; stats.contigs += components; }
        if (molecules.size() > components) ++stats.groups_split;
        const std::string umi = unpack_barcode(key, stats.umi_length);
        const std::string cell =
            stats.cell_length ? unpack_barcode(cell_key, stats.cell_length) : std::string();
        // The count is the true depth of the molecule whenever the group produced exactly one --
        // capping the reads that were consensed must never cap the reads that were counted. When
        // the group split or fragmented, each part reports what it actually held, because there is
        // no honest way to divide the reads above the cap between them.
        const bool whole = molecules.size() == 1 && components == 1;
        for (size_t m = 0; m < molecules.size(); ++m) {
            const Consensus& c = molecules[m];
            const uint64_t reported_reads = whole ? depth : c.reads;
            ++stats.molecules;
            double q = 0.0;
            for (char ch : c.qual) q += phred_from_char(ch);
            qual_sum += q;
            qual_bases += c.qual.size();
            err_sum += c.mean_error;
            // The name carries everything a downstream tool needs even when it drops the comment
            // -- dnaio does, so arda only ever sees the name.
            std::string name = stats.sample_id + "." +
                               (cell.empty() ? std::string() : cell + ".") + umi;
            if (c.components > 1) name += ".c" + std::to_string(c.component + 1);
            if (molecules.size() > c.components) name += "." + std::to_string(m + 1);
            std::string tags = "RX:Z:" + umi + "\tBC:Z:" + stats.sample_id;
            if (!cell.empty()) tags += "\tCB:Z:" + cell;
            tags += "\tMI:Z:" + name + "\tcD:i:" + std::to_string(reported_reads);
            writer.write(name, tags, c.seq, c.qual);
            std::fprintf(table, "%s\t%s\t%u\t%u\t%zu\t%" PRIu64 "\t%u\t%zu\t%.2f\t%.3e\t%.2f\n",
                         cell.empty() ? "." : cell.c_str(), umi.c_str(), c.component + 1,
                         c.components, m + 1, reported_reads, c.support, c.seq.size(),
                         c.qual.empty() ? 0.0 : q / c.qual.size(), c.mean_error, c.linkage);
            if (c.support) {
                support_reads += c.support;
                support_of_reads += c.reads;  // the reads the vote was taken over, not the depth
            }
        }
    };

    for (size_t i = 0; i < n_buckets; ++i) {
        const std::filesystem::path path =
            temp_dir / ("bucket." + std::to_string(i) + ".mig");
        if (!std::filesystem::exists(path)) continue;
        std::vector<Resident> records;
        {
            MigReader reader(path.string());
            MigRecord rec;
            while (reader.next(rec)) {
                records.push_back({rec.cell, rec.umi, rec.src_index, std::string(rec.seq1),
                                   std::string(rec.qual1)});
            }
        }
        std::sort(records.begin(), records.end());
        std::vector<Resident> group;
        for (size_t j = 0; j < records.size(); ++j) {
            if (!group.empty() && (records[j].umi != group.front().umi ||
                                   records[j].cell != group.front().cell)) {
                emit(group.front().cell, group.front().umi, group);
                group.clear();
            }
            group.push_back(std::move(records[j]));
        }
        if (!group.empty()) emit(group.front().cell, group.front().umi, group);
        std::filesystem::remove(path);
    }

    std::fclose(table);
    writer.close();
    std::error_code ec;
    std::filesystem::remove(temp_dir, ec);

    // The same barcode_space() checkout reports, on the same barcodes, so the two runs cannot
    // disagree about what the library was.
    UmiComposition comp;
    comp.length = static_cast<int>(usage.size());
    for (const std::array<uint64_t, 4>& u : usage) {
        const double n = static_cast<double>(u[0] + u[1] + u[2] + u[3]);
        comp.freq.push_back({n ? u[0] / n : 0.0, n ? u[1] / n : 0.0, n ? u[2] / n : 0.0,
                             n ? u[3] / n : 0.0});
    }
    stats.space = barcode_space(comp, stats.groups);
    // E[k | k >= 1] for a Poisson-occupied space: how many molecules a group holds on average,
    // given that it holds any. 1.0 for a barcode long enough to be unique.
    //
    // barcode_space() declines lambda past its saturation threshold, because inverting a nearly
    // full space would report "no collisions" for the most collided library there can be. But
    // that is exactly the case this number exists to warn about, so when it is declined the
    // occupancy is clamped just short of full and the result is reported as the lower bound it
    // is -- `saturated` travels with it and says so.
    double occupancy = stats.space.occupancy;
    if (occupancy > 0.99) occupancy = 0.99;
    const double lam = stats.space.saturated ? -std::log1p(-occupancy) : stats.space.lambda;
    stats.expected_molecules_per_group = lam > 0.0 ? lam / (1.0 - std::exp(-lam)) : 1.0;

    stats.mean_quality = qual_bases ? qual_sum / static_cast<double>(qual_bases) : 0.0;
    stats.mean_consensus_error =
        stats.molecules ? err_sum / static_cast<double>(stats.molecules) : 0.0;
    stats.mean_support = support_reads && stats.molecules
                             ? static_cast<double>(support_reads) /
                                   static_cast<double>(support_of_reads)
                             : 0.0;
    stats.wall_seconds = clock.seconds();
    return stats;
}

}  // namespace migec
