#include "migec/refine.hpp"

#include <algorithm>
#include <cstdio>
#include <filesystem>

#include "migec/fastq.hpp"
#include "migec/resource.hpp"
#include "migec/types.hpp"

namespace migec {

namespace {

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

void bin(std::vector<uint64_t>& histogram, uint64_t n) {
    size_t b = 0;
    while ((n >> b) > 1) ++b;
    if (histogram.size() <= b) histogram.resize(b + 1, 0);
    ++histogram[b];
}

}  // namespace

RefineStats refine(const RefineRequest& request) {
    Stopwatch clock;
    RefineStats stats;
    stats.sample_id = request.sample_id;

    std::filesystem::path out_dir(request.output_dir);
    std::filesystem::create_directories(out_dir);

    // ---------------------------------------------------------------- pass 1: the barcode table
    // Counts only. The entry array is not final until every read has been seen, and the evidence
    // is indexed in parallel with it, so it cannot be filled in the same pass.
    UmiCounts counts(0);
    {
        FastqReader reader(request.input);
        FastqRecord rec;
        bool started = false;
        while (reader.next(rec)) {
            ++stats.reads;
            const std::string_view umi = tag_value(rec.comment, "RX:Z:");
            if (umi.empty()) { ++stats.reads_without_umi; continue; }
            if (!started) {
                stats.umi_length = static_cast<int>(umi.size());
                if (stats.umi_length > kMaxBarcodeLen) {
                    throw MigecError("refine: a " + std::to_string(stats.umi_length) +
                                     " nt UMI does not fit the packed key (max " +
                                     std::to_string(kMaxBarcodeLen) + ")");
                }
                counts = UmiCounts(stats.umi_length);
                started = true;
                if (stats.sample_id.empty()) {
                    const std::string_view bc = tag_value(rec.comment, "BC:Z:");
                    stats.sample_id = bc.empty() ? std::string("sample") : std::string(bc);
                }
            } else if (static_cast<int>(umi.size()) != stats.umi_length) {
                throw MigecError("refine: read '" + std::string(rec.name) + "' carries a " +
                                 std::to_string(umi.size()) + " nt UMI where the file started "
                                 "with " + std::to_string(stats.umi_length) +
                                 " -- these are two runs concatenated, not one sample");
            }
            counts.add(pack_barcode(umi));
        }
    }
    if (!stats.umi_length) throw MigecError("refine: no read carried an RX:Z: tag");

    const std::vector<UmiCounts::Entry>& entries = counts.entries();
    stats.barcodes = entries.size();
    const int L = stats.umi_length;

    auto slot = [&entries](uint64_t key) {
        auto it = std::lower_bound(entries.begin(), entries.end(), key,
                                   [](const UmiCounts::Entry& e, uint64_t k) { return e.key < k; });
        return static_cast<size_t>(it - entries.begin());
    };

    // ---------------------------------------------------------------- pass 2: the evidence
    BarcodeEvidence evidence;
    const int pw = request.use_payload ? request.payload_width : 0;
    if (request.use_quality) {
        evidence.position_error.assign(entries.size() * static_cast<size_t>(L), 0.0f);
    }
    if (pw > 0) {
        evidence.payload_width = pw;
        evidence.payload.assign(entries.size() * static_cast<size_t>(pw), 0);
    }
    if (request.use_quality || pw > 0) {
        FastqReader reader(request.input);
        FastqRecord rec;
        while (reader.next(rec)) {
            const std::string_view umi = tag_value(rec.comment, "RX:Z:");
            if (umi.empty()) continue;
            const size_t i = slot(pack_barcode(umi));
            if (request.use_quality) {
                const std::string_view qx = tag_value(rec.comment, "QX:Z:");
                for (int j = 0; j < L; ++j) {
                    // No QX means the instrument's quality was not carried; fall back to the
                    // global rate by leaving the accumulator at zero, which `correct_umis` reads
                    // as "unknown" only if the whole barcode is zero.
                    const double e = j < static_cast<int>(qx.size())
                                         ? phred_error(phred_from_char(qx[static_cast<size_t>(j)]))
                                         : 0.0;
                    evidence.position_error[i * static_cast<size_t>(L) +
                                            static_cast<size_t>(j)] += static_cast<float>(e);
                }
            }
            if (pw > 0 && evidence.payload[i * static_cast<size_t>(pw)] == 0) {
                // The first read of a barcode is the draft. A modal base per column would need a
                // counter per barcode per column, and the draft is only ever compared to another
                // draft -- it is telling two molecules apart, not calling variants.
                const size_t n = std::min<size_t>(static_cast<size_t>(pw), rec.seq.size());
                std::copy_n(rec.seq.begin(), n,
                            evidence.payload.begin() +
                                static_cast<std::ptrdiff_t>(i * static_cast<size_t>(pw)));
            }
        }
        if (request.use_quality) {
            for (size_t i = 0; i < entries.size(); ++i) {
                for (int j = 0; j < L; ++j) {
                    evidence.position_error[i * static_cast<size_t>(L) + static_cast<size_t>(j)] /=
                        static_cast<float>(entries[i].count);
                }
            }
        }
    }
    stats.table_bytes = static_cast<uint64_t>(entries.size()) *
                        (sizeof(UmiCounts::Entry) +
                         (request.use_quality ? sizeof(float) * static_cast<size_t>(L) : 0) +
                         static_cast<size_t>(pw));

    // ---------------------------------------------------------------- correct
    const CorrectionResult correction = correct_umis(counts, request.correction, evidence);
    stats.merged = correction.merged;
    stats.merged_reads = correction.merged_reads;
    stats.merged_by_payload = correction.merged_by_payload;
    stats.molecules = correction.molecules_observed;
    stats.molecules_corrected = correction.molecules_corrected;
    stats.estimated_error = correction.estimated_error;
    stats.payload_clonality = correction.payload_clonality;
    stats.saturated = correction.saturated;
    for (uint32_t c : correction.corrected) {
        if (c > 0) bin(stats.size_histogram, c);
    }

    // ---------------------------------------------------------------- pass 3: rewrite
    const std::filesystem::path fastq_path =
        out_dir / (stats.sample_id + ".fq" + (request.gzip_level > 0 ? ".gz" : ""));
    {
        FastqWriter writer(fastq_path.string(), request.gzip_level);
        FastqReader reader(request.input);
        FastqRecord rec;
        while (reader.next(rec)) {
            const std::string_view umi = tag_value(rec.comment, "RX:Z:");
            if (umi.empty()) continue;
            const size_t i = slot(pack_barcode(umi));
            const uint32_t root = correction.root[i];
            std::string comment(rec.comment);
            if (root != i) {
                // Rewrite RX to the parent, and record what it was. A correction that cannot be
                // audited is a correction nobody can check.
                const std::string corrected = unpack_barcode(entries[root].key, L);
                const size_t at = comment.find("RX:Z:");
                comment.replace(at + 5, umi.size(), corrected);
                comment += "\tOX:Z:";
                comment += umi;
            }
            writer.write(rec.name, comment, rec.seq, rec.qual);
        }
        writer.close();
    }

    // ---------------------------------------------------------------- tables
    {
        std::FILE* fh =
            std::fopen((out_dir / (stats.sample_id + ".barcodes.tsv")).string().c_str(), "w");
        if (!fh) throw MigecError("refine: cannot write the barcode table");
        std::fprintf(fh, "umi\treads\tcorrected_reads\tparent\n");
        for (size_t i = 0; i < entries.size(); ++i) {
            const std::string umi = unpack_barcode(entries[i].key, L);
            const uint32_t root = correction.root[i];
            std::fprintf(fh, "%s\t%u\t%u\t%s\n", umi.c_str(), entries[i].count,
                         correction.corrected[i],
                         root == i ? "." : unpack_barcode(entries[root].key, L).c_str());
        }
        std::fclose(fh);
    }

    stats.wall_seconds = clock.seconds();
    return stats;
}

}  // namespace migec
