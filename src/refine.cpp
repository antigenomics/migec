#include "migec/refine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
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

// The molecule's whole key: cell barcode then UMI, so that the same UMI in two cells is two
// keys. `cell` is empty on a bulk library and the key is just the UMI.
uint64_t pack_key(std::string_view cell, std::string_view umi) {
    if (cell.empty()) return pack_barcode(umi);
    std::string joined;
    joined.reserve(cell.size() + umi.size());
    joined.append(cell);
    joined.append(umi);
    return pack_barcode(joined);
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
            const std::string_view cell = tag_value(rec.comment, "CB:Z:");
            if (!started) {
                stats.umi_length = static_cast<int>(umi.size());
                stats.cell_length = static_cast<int>(cell.size());
                if (stats.umi_length + stats.cell_length > kMaxBarcodeLen) {
                    throw MigecError(
                        "refine: a " + std::to_string(stats.cell_length) + " nt cell barcode plus "
                        "a " + std::to_string(stats.umi_length) +
                        " nt UMI does not fit the packed key (max " +
                        std::to_string(kMaxBarcodeLen) + ")");
                }
                counts = UmiCounts(stats.umi_length + stats.cell_length);
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
            counts.add(pack_key(cell, umi));
        }
    }
    if (!stats.umi_length) throw MigecError("refine: no read carried an RX:Z: tag");

    const std::vector<UmiCounts::Entry>& entries = counts.entries();
    stats.barcodes = entries.size();
    // The whole barcode: cell then UMI. Correction walks its 3L neighbourhood, so a substitution
    // in either part is found and neither is corrected across the other.
    const int L = stats.umi_length + stats.cell_length;

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
            const std::string_view cell = tag_value(rec.comment, "CB:Z:");
            const size_t i = slot(pack_key(cell, umi));
            if (request.use_quality) {
                // The whole key's quality, cell part first, to line up with the packed key.
                std::string qx(tag_value(rec.comment, "CY:Z:"));
                qx.append(tag_value(rec.comment, "QX:Z:"));
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
            const std::string_view cell = tag_value(rec.comment, "CB:Z:");
            const size_t i = slot(pack_key(cell, umi));
            const uint32_t root = correction.root[i];
            std::string comment(rec.comment);
            if (root != i) {
                // Rewrite RX (and CB, if the substitution landed there) to the parent, and record
                // what they were. A correction that cannot be audited is one nobody can check.
                const std::string whole = unpack_barcode(entries[root].key, L);
                const std::string new_cell = whole.substr(0, static_cast<size_t>(stats.cell_length));
                const std::string new_umi = whole.substr(static_cast<size_t>(stats.cell_length));
                if (new_umi != umi) {
                    const size_t at = comment.find("RX:Z:");
                    comment.replace(at + 5, umi.size(), new_umi);
                    comment += "\tOX:Z:";
                    comment += umi;
                }
                if (!cell.empty() && new_cell != cell) {
                    const size_t at = comment.find("CB:Z:");
                    comment.replace(at + 5, cell.size(), new_cell);
                    comment += "\tOC:Z:";
                    comment += cell;
                }
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
        std::fprintf(fh, "cell\tumi\treads\tcorrected_reads\tparent\n");
        for (size_t i = 0; i < entries.size(); ++i) {
            const std::string whole = unpack_barcode(entries[i].key, L);
            const std::string cell = whole.substr(0, static_cast<size_t>(stats.cell_length));
            const std::string umi = whole.substr(static_cast<size_t>(stats.cell_length));
            const uint32_t root = correction.root[i];
            std::fprintf(fh, "%s\t%s\t%u\t%u\t%s\n", cell.empty() ? "." : cell.c_str(),
                         umi.c_str(), entries[i].count, correction.corrected[i],
                         root == i ? "." : unpack_barcode(entries[root].key, L).c_str());
        }
        std::fclose(fh);
    }

    // ---------------------------------------------------------------- diagnostics
    // Three tables, all drawn from what is already held. Plots are never made in here: a figure
    // has to be redrawable from a committed TSV.
    {
        // The barcode-rank curve, Cell Ranger's plot. Log-spaced ranks, because the full curve is
        // one row per barcode and that is hundreds of millions of them for a figure that is read
        // on a log axis anyway.
        std::vector<uint32_t> sorted;
        sorted.reserve(entries.size());
        for (uint32_t c : correction.corrected) {
            if (c > 0) sorted.push_back(c);
        }
        std::sort(sorted.begin(), sorted.end(), std::greater<uint32_t>());
        uint64_t total_reads = 0;
        for (uint32_t c : sorted) total_reads += c;

        std::FILE* fh =
            std::fopen((out_dir / (stats.sample_id + ".rank.tsv")).string().c_str(), "w");
        if (!fh) throw MigecError("refine: cannot write the rank table");
        std::fprintf(fh, "rank\treads\tcumulative_reads\tcumulative_fraction\n");
        uint64_t cum = 0;
        size_t next = 0;
        for (size_t i = 0; i < sorted.size(); ++i) {
            cum += sorted[i];
            if (i == next || i + 1 == sorted.size()) {
                std::fprintf(fh, "%zu\t%u\t%llu\t%.6f\n", i + 1, sorted[i],
                             static_cast<unsigned long long>(cum),
                             total_reads ? static_cast<double>(cum) /
                                               static_cast<double>(total_reads) : 0.0);
                next = std::max(next + 1, static_cast<size_t>(static_cast<double>(next) * 1.05));
            }
        }
        std::fclose(fh);
    }
    {
        // Per MIG-size bin: how much of it was error, and how diverse the sequence is there.
        // Error children pile up at low counts; finding them at high counts means the correction
        // is merging real molecules, and the sequence entropy is what says whether a bin holds
        // one artefact repeated or a real population.
        const size_t nbins = stats.size_histogram.size() + 1;
        std::vector<uint64_t> barcodes(nbins, 0), merged_in(nbins, 0), reads_in(nbins, 0);
        std::vector<std::vector<std::array<uint32_t, 4>>> comp(nbins);
        for (size_t i = 0; i < entries.size(); ++i) {
            size_t b = 0;
            while ((entries[i].count >> b) > 1) ++b;
            if (b >= nbins) b = nbins - 1;
            ++barcodes[b];
            reads_in[b] += entries[i].count;
            if (correction.root[i] != i) ++merged_in[b];
            if (pw > 0) {
                if (comp[b].empty()) comp[b].assign(static_cast<size_t>(pw), {});
                for (int j = 0; j < pw; ++j) {
                    const uint8_t code =
                        base_code(evidence.payload[i * static_cast<size_t>(pw) +
                                                   static_cast<size_t>(j)]);
                    if (code != kInvalidBase) ++comp[b][static_cast<size_t>(j)][code];
                }
            }
        }
        std::FILE* fh =
            std::fopen((out_dir / (stats.sample_id + ".bins.tsv")).string().c_str(), "w");
        if (!fh) throw MigecError("refine: cannot write the bin table");
        std::fprintf(fh, "min_reads\tmax_reads\tbarcodes\treads\tmerged\tfraction_erroneous\t"
                         "payload_entropy_bits\n");
        for (size_t b = 0; b < nbins; ++b) {
            if (!barcodes[b]) continue;
            double entropy = 0.0;
            int scored = 0;
            for (const std::array<uint32_t, 4>& col : comp[b]) {
                const double n_col = col[0] + col[1] + col[2] + col[3];
                if (n_col < 2) continue;
                double h = 0.0;
                for (uint32_t v : col) {
                    if (!v) continue;
                    const double q = static_cast<double>(v) / n_col;
                    h -= q * std::log2(q);
                }
                entropy += h;
                ++scored;
            }
            std::fprintf(fh, "%llu\t%llu\t%llu\t%llu\t%llu\t%.6f\t%.4f\n",
                         1ull << b, (1ull << (b + 1)) - 1,
                         static_cast<unsigned long long>(barcodes[b]),
                         static_cast<unsigned long long>(reads_in[b]),
                         static_cast<unsigned long long>(merged_in[b]),
                         static_cast<double>(merged_in[b]) / static_cast<double>(barcodes[b]),
                         scored ? entropy / scored : 0.0);
        }
        std::fclose(fh);
    }

    stats.wall_seconds = clock.seconds();
    return stats;
}

}  // namespace migec
