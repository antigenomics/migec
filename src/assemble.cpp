#include "migec/assemble.hpp"

#include <algorithm>
#include <cinttypes>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <memory>

#include "migec/fastq.hpp"
#include "migec/mig_record.hpp"
#include "migec/parallel.hpp"
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
    // Never: the floor is a CONSTANT, never a function of the thread count. The buckets are the
    // unit of parallelism in pass 2, so there has to be more than one of them -- but if -t chose
    // how many, then -t would choose the gzip member boundaries too, and two runs at different
    // thread counts would produce byte-different files holding identical records. A fixed floor
    // makes the partition a property of the input alone, and -t only decides who chews on it.
    int floor_bits = kMinBucketBits;

    std::error_code ec;
    const uintmax_t on_disk = std::filesystem::file_size(path, ec);
    if (ec) return floor_bits;
    // MEASURED, not estimated: 8x was a guess and it was 2.4x too small. 4 M reads of 90 nt in a
    // 67 MB gzip hold 1.28 GB across the buckets of pass 2 -- 19x on disk -- because a `Resident`
    // is two heap std::strings (each with its allocator header and its rounded-up bucket) plus
    // three 8-byte keys, not just the 180 bytes of payload. Guessing low is the expensive
    // direction: it picks too few buckets, and pass 2 holds `kBucketConcurrency` of them at once.
    // 20x rounds the measurement up.
    uint64_t resident = static_cast<uint64_t>(on_disk) * 20;
    // Per bucket, so that kBucketConcurrency of them are resident inside the total budget.
    const uint64_t per_bucket = kBucketBudgetBytes / static_cast<uint64_t>(kBucketConcurrency);
    int bits = floor_bits;
    while (bits < kMaxBucketBits && (resident >> bits) > per_bucket) ++bits;
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

    const int threads = worker_count(request.threads, 64);
    stats.threads = threads;
    const int bits = request.bucket_bits > 0 ? request.bucket_bits
                                             : choose_bucket_bits(request.input);
    const size_t n_buckets = static_cast<size_t>(1) << bits;
    stats.buckets = static_cast<int>(n_buckets);

    // ------------------------------------------------------------------ pass 1: partition
    //
    // Measured before it was threaded: 2.07 s of a 2.69 s run on 4 M reads, against a 0.23 s
    // `gzip -dc` floor for the same file. So five sixths of the partition is not the inflate --
    // it is the tag scan, the barcode packing, the record serialisation and the level-1 deflate of
    // each bucket block. All four move to the workers here, and the reader stays where it has to.
    //
    // Ownership, not locking: worker w owns every bucket with `bucket % threads == w` for the whole
    // run, so a bucket file has exactly one writer and no bucket state is ever shared. Records
    // reach a bucket in input order because the chunks are consumed in order and each worker walks
    // its chunk forwards, which is what keeps the bytes identical to the serial version -- and
    // identical at any `-t`, since ownership decides *who* writes a record and never *which* file
    // it lands in or *where* in that file.
    {
        Stopwatch partition_clock;
        std::vector<Bucket> buckets(n_buckets);
        const size_t block_bytes = std::clamp(
            static_cast<size_t>(kWriterBudgetBytes / n_buckets), kMinBlockBytes, kMaxBlockBytes);
        // A CONSTANT, for the same reason the bucket count is: the chunking must not be a function
        // of `-t`. ONE chunk is held rather than one per worker, so this is ~2 MB whatever the
        // thread count.
        //
        // ponytail: 8192 is where a real trade sits, and it is worth naming. `parallel_for` starts
        // and joins its threads per call, so a bigger chunk amortises that: 64 k reads runs 26%
        // faster (2.51 M reads/s against 1.99 M) but costs 16 MB of resident chunk, and that is
        // enough to make PASS 1 the memory peak on a finely partitioned shallow library -- which
        // breaks the property that a finer partition costs less, not more, and
        // `test_shallow_memory_is_still_bounded_by_the_bucket` catches exactly that. The upgrade
        // path is a persistent worker pool rather than a bigger chunk: with the threads started
        // once for the whole pass, chunk size stops buying anything and can stay small.
        constexpr size_t kChunkReads = 8192;
        // What the parallel scan extracts from a record. No strings: the serial pass that follows
        // only validates and counts, so it must not have to touch the comment again.
        struct Parsed {
            uint64_t umi_key = 0, cell_key = 0, src_index = 0;
            uint32_t umi_len = 0, cell_len = 0;
            uint16_t flags = 0;
            uint32_t bucket = 0;
            bool has_umi = false;
        };
        // Sized once and ASSIGNED into, never cleared: `clear()` destroys the four strings of every
        // record, so a fresh chunk costs 4 x kChunkReads allocations and the reader spends its time
        // in malloc rather than in inflate. `assign` reuses the capacity a previous chunk left.
        std::vector<FastqOwned> chunk(kChunkReads);
        std::vector<Parsed> parsed(kChunkReads);
        size_t held = 0;

        FastqReader reader(request.input);
        FastqRecord rec;
        uint64_t index = 0;
        IntakeLimit limit = request.limit;
        bool eof = false, stopped = false;

        while (!eof && !stopped) {
            held = 0;
            while (held < kChunkReads && reader.next(rec)) {
                FastqOwned& slot = chunk[held++];
                slot.name.assign(rec.name);
                slot.comment.assign(rec.comment);
                slot.seq.assign(rec.seq);
                slot.qual.assign(rec.qual);
            }
            if (held < kChunkReads) eof = true;
            if (held == 0) break;
            std::fill(parsed.begin(), parsed.begin() + static_cast<ptrdiff_t>(held), Parsed{});

            parallel_for(held, threads, [&](size_t i, int) {
                const FastqOwned& r = chunk[i];
                Parsed& p = parsed[i];
                const std::string_view umi = tag_value(r.comment, "RX:Z:");
                if (umi.empty()) return;
                p.has_umi = true;
                p.umi_len = static_cast<uint32_t>(umi.size());
                const std::string_view cell = tag_value(r.comment, "CB:Z:");
                p.cell_len = static_cast<uint32_t>(cell.size());
                bool has_n = false, cell_has_n = false;
                p.umi_key = pack_barcode(umi, &has_n);
                p.cell_key = cell.empty() ? 0 : pack_barcode(cell, &cell_has_n);
                p.flags = static_cast<uint16_t>(kSingleEnd | (has_n ? kUmiHasN : 0) |
                                                (cell_has_n ? kCellHasN : 0));
                // Partition on the cell when there is one: every read of a cell then lands in one
                // bucket, which is what makes a per-cell scope local. Without cells the UMI is the
                // whole key and the partition is on it.
                p.bucket = bucket_of(cell.empty() ? p.umi_key : p.cell_key, bits);
            });

            // Serial: everything order-dependent. Counters, the length agreement that says two
            // runs were concatenated, the intake limit, and `src_index` -- which is the read's
            // position in the input and so cannot be handed out by a worker.
            size_t emit = 0;
            for (; emit < held; ++emit) {
                Parsed& p = parsed[emit];
                ++stats.reads;
                if (!p.has_umi) { ++stats.reads_without_umi; continue; }
                if (stats.umi_length == 0) {
                    stats.umi_length = static_cast<int>(p.umi_len);
                    if (stats.umi_length > kMaxBarcodeLen) {
                        throw MigecError("assemble: a " + std::to_string(stats.umi_length) +
                                         " nt UMI does not fit the packed key (max " +
                                         std::to_string(kMaxBarcodeLen) + ")");
                    }
                } else if (static_cast<int>(p.umi_len) != stats.umi_length) {
                    throw MigecError("assemble: read '" + chunk[emit].name + "' carries a " +
                                     std::to_string(p.umi_len) + " nt UMI where the file started " +
                                     "with " + std::to_string(stats.umi_length) +
                                     " -- these are two runs concatenated, not one sample");
                }
                if (stats.sample_id.empty()) {
                    const std::string_view bc = tag_value(chunk[emit].comment, "BC:Z:");
                    stats.sample_id = bc.empty() ? std::string("sample") : std::string(bc);
                }
                if (stats.cell_length == 0 && p.cell_len) {
                    stats.cell_length = static_cast<int>(p.cell_len);
                }
                // Checked after the key is known, so `--limit-umi` counts barcodes, not reads.
                const uint64_t partition_key = p.cell_len ? p.cell_key : p.umi_key;
                if (limit.active() && !limit.admit(stats.reads, partition_key)) {
                    --stats.reads;
                    stats.limited = true;
                    stopped = true;
                    break;
                }
                p.src_index = index++;
                ++buckets[p.bucket].records;
            }

            // The writers are opened here, on the driver, so that a bucket's header carries the
            // sample id and the barcode lengths the serial pass above has just settled.
            for (size_t b = 0; b < n_buckets; ++b) {
                if (buckets[b].records == 0 || buckets[b].writer) continue;
                MigHeader header;
                header.umi_len = static_cast<uint8_t>(stats.umi_length);
                header.cell_len = static_cast<uint8_t>(stats.cell_length);
                header.bucket_index = static_cast<uint8_t>(b);
                header.bucket_bits = static_cast<uint8_t>(bits);
                header.paired = false;
                header.sample_id = stats.sample_id;
                buckets[b].path = (temp_dir / ("bucket." + std::to_string(b) + ".mig")).string();
                buckets[b].writer =
                    std::make_unique<MigWriter>(buckets[b].path, header, block_bytes);
            }

            const size_t written = emit;
            parallel_for(static_cast<size_t>(threads), threads, [&](size_t w, int) {
                for (size_t i = 0; i < written; ++i) {
                    const Parsed& p = parsed[i];
                    if (!p.has_umi || p.bucket % static_cast<uint32_t>(threads) != w) continue;
                    MigRecord out;
                    out.cell = p.cell_key;
                    out.umi = p.umi_key;
                    out.src_index = p.src_index;
                    out.flags = p.flags;
                    out.seq1 = chunk[i].seq;
                    out.qual1 = chunk[i].qual;
                    buckets[p.bucket].writer->write(out);
                }
            });
        }
        for (Bucket& b : buckets) {
            if (b.writer) b.writer->close();
        }
        stats.partition_seconds = partition_clock.seconds();
        // Bucket paths are regenerated below from the index, so the writers can go now.
    }

    // ------------------------------------------------------------------ pass 2: consense
    //
    // One bucket at a time was the memory bound; one bucket per THREAD is the same bound times the
    // thread count, and the buckets are independent by construction -- every read of a barcode is
    // in exactly one of them, because the partition is on that barcode. So this loop parallelises
    // with no shared state at all: each worker owns its output files and its own counters.
    //
    // Never: the output must not depend on -t. Workers write per-bucket temporary files which are
    // concatenated in BUCKET order afterwards, and bucket order is key order, so the consensus
    // FASTQ comes out sorted by barcode whatever the thread count and whatever order the buckets
    // happened to finish in.
    struct BucketOut {
        std::string fastq_path, table_path;
        uint64_t groups = 0, molecules = 0, groups_split = 0, groups_fragmented = 0, contigs = 0;
        uint64_t reads_dropped = 0, groups_capped = 0, reads_over_cap = 0;
        std::vector<uint64_t> size_histogram;
        // Molecules per (depth bin, rounded Phred). Both axes are DISCRETE and small -- a depth
        // bin is a power of two and a Phred is an integer capped at the RT floor -- so the whole
        // joint distribution fits in a fixed grid and there is no reason to thin a scatter plot
        // and hope. `kQualityLevels` covers Q0..Q60; anything above is clamped, and the cap means
        // nothing lands there anyway.
        std::vector<std::array<uint64_t, kQualityLevels>> quality_grid;
        double qual_sum = 0.0, err_sum = 0.0;
        uint64_t qual_bases = 0, support_reads = 0, support_of_reads = 0;
        // Per-position base usage over DISTINCT barcodes, which is what the birthday arithmetic
        // needs -- weighting by reads would let one over-amplified molecule set the composition.
        std::vector<std::array<uint64_t, 4>> usage;
    };
    std::vector<BucketOut> outs(n_buckets);

    // Merged after the loop, in bucket order.
    double qual_sum = 0.0, err_sum = 0.0;
    uint64_t qual_bases = 0, support_reads = 0, support_of_reads = 0;
    std::vector<std::array<uint64_t, kQualityLevels>> quality_grid;
    std::vector<std::array<uint64_t, 4>> usage;

    auto consense_bucket = [&](size_t bucket, BucketOut& out) {
        std::unique_ptr<FastqWriter> writer;
        std::FILE* table = nullptr;
        auto emit = [&](uint64_t cell_key, uint64_t key, const std::vector<Resident>& group) {
        ++out.groups;
        {
            const std::string bases = unpack_barcode(key, stats.umi_length);
            if (out.usage.empty()) out.usage.assign(bases.size(), {});
            for (size_t j = 0; j < bases.size(); ++j) {
                const uint8_t c = base_code(bases[j]);
                if (c != kInvalidBase) ++out.usage[j][c];
            }
        }
        bin(out.size_histogram, group.size());
        if (group.size() < request.min_reads) {
            out.reads_dropped += group.size();
            return;
        }
        // The consensus sees at most kMaxReadsPerGroup of them; the molecule keeps all of them as
        // its depth. Records arrive in src_index order, which is input order and carries no
        // relation to the sequence, so the first N are an unbiased sample of the group and cost
        // nothing to take.
        const size_t depth = group.size();
        const size_t used = std::min(depth, kMaxReadsPerGroup);
        if (used < depth) {
            ++out.groups_capped;
            out.reads_over_cap += depth - used;
        }
        std::vector<ConsensusRead> reads;
        reads.reserve(used);
        for (size_t i = 0; i < used; ++i) reads.push_back({group[i].seq, group[i].qual});
        const std::vector<Consensus> molecules = assemble_group(reads, request.consensus);
        const uint32_t components = molecules.empty() ? 1 : molecules[0].components;
        if (components > 1) { ++out.groups_fragmented; out.contigs += components; }
        if (molecules.size() > components) ++out.groups_split;
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
            ++out.molecules;
            double q = 0.0;
            for (char ch : c.qual) q += phred_from_char(ch);
            out.qual_sum += q;
            out.qual_bases += c.qual.size();
            out.err_sum += c.mean_error;
            if (!c.qual.empty()) {
                size_t depth_bin = 0;
                while ((reported_reads >> depth_bin) > 1) ++depth_bin;
                if (out.quality_grid.size() <= depth_bin) {
                    out.quality_grid.resize(depth_bin + 1, {});
                }
                const double mean_q = q / static_cast<double>(c.qual.size());
                size_t level = static_cast<size_t>(mean_q + 0.5);
                if (level >= kQualityLevels) level = kQualityLevels - 1;
                ++out.quality_grid[depth_bin][level];
            }
            // The name carries everything a downstream tool needs even when it drops the comment
            // -- dnaio does, so arda only ever sees the name.
            std::string name = stats.sample_id + "." +
                               (cell.empty() ? std::string() : cell + ".") + umi;
            if (c.components > 1) name += ".c" + std::to_string(c.component + 1);
            if (molecules.size() > c.components) name += "." + std::to_string(m + 1);
            std::string tags = "RX:Z:" + umi + "\tBC:Z:" + stats.sample_id;
            if (!cell.empty()) tags += "\tCB:Z:" + cell;
            tags += "\tMI:Z:" + name + "\tcD:i:" + std::to_string(reported_reads);
            writer->write(name, tags, c.seq, c.qual);
            std::fprintf(table, "%s\t%s\t%u\t%u\t%zu\t%" PRIu64 "\t%u\t%zu\t%.2f\t%.3e\t%.2f\n",
                         cell.empty() ? "." : cell.c_str(), umi.c_str(), c.component + 1,
                         c.components, m + 1, reported_reads, c.support, c.seq.size(),
                         c.qual.empty() ? 0.0 : q / c.qual.size(), c.mean_error, c.linkage);
            if (c.support) {
                out.support_reads += c.support;
                out.support_of_reads += c.reads;  // the vote's reads, not the molecule's depth
            }
        }
    };

        const std::filesystem::path path =
            temp_dir / ("bucket." + std::to_string(bucket) + ".mig");
        if (!std::filesystem::exists(path)) return;
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

        out.fastq_path = (temp_dir / ("out." + std::to_string(bucket) + ".fq.gz")).string();
        out.table_path = (temp_dir / ("out." + std::to_string(bucket) + ".tsv")).string();
        writer = std::make_unique<FastqWriter>(out.fastq_path, request.gzip_level);
        table = std::fopen(out.table_path.c_str(), "w");
        if (!table) throw MigecError("assemble: cannot write the per-molecule table");

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
        std::fclose(table);
        table = nullptr;
        writer->close();
        writer.reset();
        std::filesystem::remove(path);
    };

    parallel_for(n_buckets, threads, [&](size_t bucket, int) {
        consense_bucket(bucket, outs[bucket]);
    });

    // ------------------------------------------------------------------ merge, in bucket order
    //
    // Bucket order is key order, so concatenating here is what makes the output sorted by barcode
    // AND independent of the thread count. Concatenated gzip members are a valid gzip stream
    // (RFC 1952 s2.2), which is the same property checkout relies on to compress on its workers.
    const std::filesystem::path fastq_path =
        out_dir / (stats.sample_id + ".consensus.fq" + (request.gzip_level > 0 ? ".gz" : ""));
    {
        std::FILE* fq = std::fopen(fastq_path.string().c_str(), "wb");
        if (!fq) throw MigecError("assemble: cannot write " + fastq_path.string());
        std::FILE* tsv =
            std::fopen((out_dir / (stats.sample_id + ".mig.tsv")).string().c_str(), "w");
        if (!tsv) throw MigecError("assemble: cannot write the per-molecule table");
        std::fprintf(tsv,
                     "cell\tumi\tcontig\tcontigs\tmolecule\treads\tsupport\tlength\t"
                     "mean_quality\tconsensus_error\tlinkage\n");
        std::vector<char> buf(1u << 20);
        for (size_t i = 0; i < n_buckets; ++i) {
            const BucketOut& o = outs[i];
            for (const auto& [src, dst] : {std::pair<const std::string&, std::FILE*>(o.fastq_path, fq),
                                           std::pair<const std::string&, std::FILE*>(o.table_path, tsv)}) {
                if (src.empty()) continue;
                std::FILE* in = std::fopen(src.c_str(), "rb");
                if (!in) continue;
                size_t n;
                while ((n = std::fread(buf.data(), 1, buf.size(), in)) > 0) {
                    if (std::fwrite(buf.data(), 1, n, dst) != n) {
                        std::fclose(in);
                        std::fclose(fq);
                        std::fclose(tsv);
                        throw MigecError("assemble: short write merging bucket outputs");
                    }
                }
                std::fclose(in);
                std::filesystem::remove(src);
            }
            stats.groups += o.groups;
            stats.molecules += o.molecules;
            stats.groups_split += o.groups_split;
            stats.groups_fragmented += o.groups_fragmented;
            stats.contigs += o.contigs;
            stats.reads_dropped += o.reads_dropped;
            stats.groups_capped += o.groups_capped;
            stats.reads_over_cap += o.reads_over_cap;
            qual_sum += o.qual_sum;
            err_sum += o.err_sum;
            qual_bases += o.qual_bases;
            support_reads += o.support_reads;
            support_of_reads += o.support_of_reads;
            for (size_t b = 0; b < o.size_histogram.size(); ++b) {
                if (stats.size_histogram.size() <= b) stats.size_histogram.resize(b + 1, 0);
                stats.size_histogram[b] += o.size_histogram[b];
            }
            if (quality_grid.size() < o.quality_grid.size()) {
                quality_grid.resize(o.quality_grid.size(), {});
            }
            for (size_t b = 0; b < o.quality_grid.size(); ++b) {
                for (size_t l = 0; l < kQualityLevels; ++l) {
                    quality_grid[b][l] += o.quality_grid[b][l];
                }
            }
            if (usage.size() < o.usage.size()) usage.resize(o.usage.size(), {});
            for (size_t j = 0; j < o.usage.size(); ++j) {
                for (int c = 0; c < 4; ++c) usage[j][static_cast<size_t>(c)] += o.usage[j][static_cast<size_t>(c)];
            }
        }
        std::fclose(fq);
        std::fclose(tsv);
    }
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

    stats.quality_grid = std::move(quality_grid);
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
