#include "migec/refine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <filesystem>
#include <unordered_map>

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

// The same key, with the cell barcode replaced by its whitelist-corrected form.
uint64_t pack_key_snapped(uint64_t snapped_cell, int cell_length, std::string_view umi) {
    if (!cell_length) return pack_barcode(umi);
    std::string joined;
    joined.reserve(static_cast<size_t>(cell_length) + umi.size());
    joined.append(unpack_barcode(snapped_cell, cell_length));
    joined.append(umi);
    return pack_barcode(joined);
}

// Surviving molecules per power-of-two bin of their CORRECTED count. A barcode's bin moves when
// its children are folded in, so this cannot be read off the pre-correction histogram -- and it is
// computed once rather than per bin, because a scan per bin is O(barcodes x bins) and the barcodes
// are the thing that scales.
std::vector<uint64_t> molecules_per_bin(const CorrectionResult& correction, size_t nbins) {
    std::vector<uint64_t> out(nbins, 0);
    for (uint32_t c : correction.corrected) {
        if (!c) continue;
        size_t idx = 0;
        while ((c >> idx) > 1) ++idx;
        if (idx >= nbins) idx = nbins - 1;
        ++out[idx];
    }
    return out;
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

    // ------------------------------------------------------- pass 0: snap cells to the whitelist
    // Before anything else, because every downstream key contains the cell barcode: correcting it
    // afterwards would mean rebuilding the whole table.
    std::unordered_map<uint64_t, uint64_t> cell_remap;  // observed key -> whitelist key
    if (!request.cell_whitelist.empty()) {
        const Whitelist list = Whitelist::load(request.cell_whitelist);
        // Distinct observed cell barcodes and their read counts. Bounded by the cells, not the
        // reads -- and on a droplet run the ambient barcodes dominate that count, which is
        // exactly the population the background prior is about.
        std::unordered_map<uint64_t, uint32_t> seen;
        std::unordered_map<uint64_t, std::string> quals;
        uint64_t total_reads = 0;
        {
            FastqReader reader(request.input);
            FastqRecord rec;
            while (reader.next(rec)) {
                const std::string_view cell = tag_value(rec.comment, "CB:Z:");
                if (cell.empty()) continue;
                if (static_cast<int>(cell.size()) != list.length()) {
                    throw MigecError(
                        "refine: the reads carry a " + std::to_string(cell.size()) +
                        " nt cell barcode and the whitelist holds " +
                        std::to_string(list.length()) + " nt entries");
                }
                const uint64_t key = pack_barcode(cell);
                ++seen[key];
                ++total_reads;
                if (quals.find(key) == quals.end()) {
                    quals.emplace(key, std::string(tag_value(rec.comment, "CY:Z:")));
                }
            }
        }
        stats.whitelist.barcodes = seen.size();

        std::vector<uint32_t> counts(list.size(), 0);
        std::vector<uint64_t> off;
        for (const auto& kv : seen) {
            const size_t at = list.index_of(kv.first);
            if (at != static_cast<size_t>(-1)) {
                counts[at] = kv.second;
                ++stats.whitelist.exact;
            } else {
                off.push_back(kv.first);
            }
        }
        // Barcodes with no distance-1 entry cannot be single substitutions of anything on the
        // list, so their reads measure how much of this library is genuinely off-list.
        uint64_t far_reads = 0;
        std::vector<uint64_t> near;
        for (uint64_t key : off) {
            bool has_neighbour = false;
            for (int j = 0; j < list.length() && !has_neighbour; ++j) {
                const int shift = 62 - 2 * j;
                const uint64_t cur = (key >> shift) & 3u;
                for (uint64_t b = 0; b < 4; ++b) {
                    if (b == cur) continue;
                    if (list.contains((key & ~(uint64_t{3} << shift)) | (b << shift))) {
                        has_neighbour = true;
                        break;
                    }
                }
            }
            if (has_neighbour) {
                near.push_back(key);
            } else {
                ++stats.whitelist.far;
                far_reads += seen[key];
            }
        }
        stats.whitelist.background_prior =
            request.whitelist.background_prior >= 0.0
                ? request.whitelist.background_prior
                : Whitelist::measure_background(far_reads, total_reads, off.size());

        for (uint64_t key : near) {
            const std::string observed = unpack_barcode(key, list.length());
            const std::string fixed =
                list.correct(observed, quals[key], counts, request.whitelist,
                             stats.whitelist.background_prior);
            if (fixed.empty()) {
                ++stats.whitelist.off_list;
            } else {
                cell_remap[key] = pack_barcode(fixed);
                ++stats.whitelist.corrected;
                stats.whitelist.reads_corrected += seen[key];
            }
        }
        stats.whitelist.off_list += stats.whitelist.far;
    }

    auto snap = [&cell_remap](uint64_t key) {
        if (cell_remap.empty()) return key;
        auto it = cell_remap.find(key);
        return it == cell_remap.end() ? key : it->second;
    };

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
            counts.add(cell.empty()
                           ? pack_barcode(umi)
                           : pack_key_snapped(snap(pack_barcode(cell)), stats.cell_length, umi));
        }
    }
    if (!stats.umi_length) throw MigecError("refine: no read carried an RX:Z: tag");

    const std::vector<UmiCounts::Entry>& entries = counts.entries();
    stats.barcodes = entries.size();
    // The whole barcode: cell then UMI. Correction walks its 3L neighbourhood, so a substitution
    // in either part is found and neither is corrected across the other.
    const int L = stats.umi_length + stats.cell_length;

    // Index of `key`, or entries.size() when it is absent. Returning the lower bound unchecked
    // would be a silent misattribution: a key past the end indexes one element off the evidence
    // arrays, and a key that merely sorts next to a real one writes another barcode's reads into
    // its slot. Pass 1 inserts every key the later passes look up, so a miss should be
    // unreachable -- this is what makes "should be" checkable instead of assumed.
    auto slot = [&entries](uint64_t key) -> size_t {
        auto it = std::lower_bound(entries.begin(), entries.end(), key,
                                   [](const UmiCounts::Entry& e, uint64_t k) { return e.key < k; });
        if (it == entries.end() || it->key != key) return entries.size();
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
            const size_t i = slot(cell.empty()
                                     ? pack_barcode(umi)
                                     : pack_key_snapped(snap(pack_barcode(cell)),
                                                        stats.cell_length, umi));
            if (i >= entries.size()) continue;
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
            const uint64_t snapped = cell.empty() ? 0 : snap(pack_barcode(cell));
            const size_t i = slot(cell.empty()
                                     ? pack_barcode(umi)
                                     : pack_key_snapped(snapped, stats.cell_length, umi));
            if (i >= entries.size()) {
                // Unreachable unless pass 1 and pass 3 disagree about the key. Pass the read
                // through untouched rather than dropping it: emitting a read whose barcode was
                // not corrected is recoverable, losing it silently is not.
                writer.write(rec.name, rec.comment, rec.seq, rec.qual);
                continue;
            }
            const uint32_t root = correction.root[i];
            std::string comment(rec.comment);
            // The final barcode is the root's, which already carries the whitelist snap because
            // the snap happened before the table was built. Compare against what the READ said,
            // not against the root only: a barcode can be snapped by the whitelist without being
            // merged by the posterior, and rewriting it only in the second case would leave the
            // CB tag disagreeing with the key everything else was grouped on.
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

    // ---------------------------------------------------------------- cell calling
    if (stats.cell_length > 0) {
        // Molecules per cell, after correction. Molecules, never reads: a cell is a set of
        // captured molecules, and read depth is amplification.
        std::vector<std::pair<uint64_t, uint32_t>> per_cell;  // (cell key, molecules)
        {
            // A barcode packs MSB-first, so base 0 of the CELL is in the top bits and the whole
            // key occupies bits 63 down to 64-2L. The cell is therefore the top 2*cell_length
            // bits -- not "everything above the UMI", which is only the same thing when cell and
            // UMI happen to fill all 32 bases.
            const int shift = 64 - 2 * stats.cell_length;
            uint64_t current = 0;
            uint32_t n = 0;
            bool started = false;
            // entries are sorted by the packed key, and the cell occupies the HIGH bits, so a
            // cell's molecules are already contiguous -- no grouping pass and no map.
            for (size_t i = 0; i < entries.size(); ++i) {
                if (correction.corrected[i] == 0) continue;
                const uint64_t cell_key = entries[i].key >> shift;
                if (!started || cell_key != current) {
                    if (started) per_cell.emplace_back(current, n);
                    current = cell_key;
                    n = 0;
                    started = true;
                }
                ++n;
            }
            if (started) per_cell.emplace_back(current, n);
        }
        stats.cells_observed = per_cell.size();

        std::vector<uint32_t> sizes;
        sizes.reserve(per_cell.size());
        for (const auto& kv : per_cell) sizes.push_back(kv.second);
        std::sort(sizes.begin(), sizes.end(), std::greater<uint32_t>());

        if (!sizes.empty()) {
            // OrdMag: the 99th percentile of the top `expect_cells`, over ten.
            const size_t top = std::min<size_t>(static_cast<size_t>(std::max(1, request.expect_cells)),
                                                sizes.size());
            const size_t at = static_cast<size_t>(0.01 * static_cast<double>(top));
            const uint32_t robust_max = sizes[std::min(at, top - 1)];
            stats.cell_threshold = std::max<uint32_t>(1, robust_max / 10);

            // The knee, reported next to it rather than instead of it: the rank furthest from the
            // chord joining the first and last points of the log-log curve. It is a description of
            // the data; OrdMag is the rule that makes the call, and when the two disagree badly
            // that is worth seeing.
            if (sizes.size() > 2) {
                const double x0 = 0.0, y0 = std::log10(static_cast<double>(sizes.front()));
                const double x1 = std::log10(static_cast<double>(sizes.size()));
                const double y1 = std::log10(static_cast<double>(sizes.back()));
                double best = -1.0;
                for (size_t i = 1; i + 1 < sizes.size(); ++i) {
                    const double x = std::log10(static_cast<double>(i + 1));
                    const double y = std::log10(static_cast<double>(sizes[i]));
                    const double d = std::fabs((y1 - y0) * x - (x1 - x0) * y + x1 * y0 - y1 * x0);
                    if (d > best) {
                        best = d;
                        stats.knee_rank = i + 1;
                        stats.knee_molecules = sizes[i];
                    }
                }
            }
            for (uint32_t v : sizes) {
                if (v >= stats.cell_threshold) {
                    ++stats.cells_called;
                    stats.molecules_in_called += v;
                }
            }
        }

        std::FILE* fh =
            std::fopen((out_dir / (stats.sample_id + ".cells.tsv")).string().c_str(), "w");
        if (!fh) throw MigecError("refine: cannot write the cell table");
        std::fprintf(fh, "cell\tmolecules\tcalled\n");
        for (const auto& kv : per_cell) {
            std::fprintf(fh, "%s\t%u\t%d\n",
                         unpack_barcode(kv.first << (64 - 2 * stats.cell_length),
                                        stats.cell_length).c_str(),
                         kv.second, kv.second >= stats.cell_threshold ? 1 : 0);
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
        // What is left over: a surviving barcode that still looks like a child of a surviving
        // neighbour is one the posterior declined to merge. Counting those per bin gives a
        // residual false-molecule rate measured on this library rather than derived.
        //
        // ⚠ "Much larger neighbour" alone is not the test. At 1-3 reads per UMI nothing is 20x
        // anything, so a count-ratio criterion reports zero residual in precisely the regime where
        // the residual is worst -- the same trap the correction posterior itself fell into. The
        // payload is what still separates them at one read: a neighbour whose reads agree on the
        // molecule is a child whatever the counts say.
        std::vector<uint64_t> suspected(nbins, 0);
        {
            const int width = L;
            for (size_t i = 0; i < entries.size(); ++i) {
                if (correction.corrected[i] == 0) continue;
                const uint32_t mine = correction.corrected[i];
                bool looks_like_a_child = false;
                for (int j = 0; j < width && !looks_like_a_child; ++j) {
                    const int shift = 62 - 2 * j;
                    const uint64_t cur = (entries[i].key >> shift) & 3u;
                    for (uint64_t b = 0; b < 4; ++b) {
                        if (b == cur) continue;
                        const uint64_t want =
                            (entries[i].key & ~(uint64_t{3} << shift)) | (b << shift);
                        const size_t at = slot(want);
                        if (at >= entries.size() || correction.corrected[at] == 0) continue;
                        if (correction.corrected[at] >= 20u * mine) {
                            looks_like_a_child = true;
                            break;
                        }
                        // ...and payload agreement is only evidence when two unrelated barcodes
                        // do NOT agree anyway. On a clonal library they do -- measured here at
                        // 0.80 on an HIV amplicon -- and using it would call almost every
                        // singleton a residual child. `correct_umis` already discounts payload
                        // agreement by exactly this number; the residual estimate has to as well.
                        if (pw > 0 && correction.payload_clonality < 0.5 &&
                            correction.corrected[at] >= mine) {
                            int mism = 0, cmp = 0;
                            for (int k = 0; k < pw; ++k) {
                                const char x = evidence.payload[i * static_cast<size_t>(pw) +
                                                                static_cast<size_t>(k)];
                                const char y = evidence.payload[at * static_cast<size_t>(pw) +
                                                                static_cast<size_t>(k)];
                                if (!x || !y || x == 'N' || y == 'N') continue;
                                ++cmp;
                                mism += x != y;
                            }
                            if (cmp >= 8 && mism * 20 <= cmp) {
                                looks_like_a_child = true;
                                break;
                            }
                        }
                    }
                }
                if (!looks_like_a_child) continue;
                size_t bidx = 0;
                while ((mine >> bidx) > 1) ++bidx;
                if (bidx >= nbins) bidx = nbins - 1;
                ++suspected[bidx];
                ++stats.suspected_residual;
            }
        }
        const std::vector<uint64_t> surviving_in = molecules_per_bin(correction, nbins);
        // The smallest size at which the residual rate is acceptable. Reported, never applied.
        for (size_t b = 0; b < nbins; ++b) {
            if (!barcodes[b]) continue;
            const uint64_t surviving = surviving_in[b];
            const double fdr = surviving ? static_cast<double>(suspected[b]) /
                                               static_cast<double>(surviving) : 0.0;
            if (b == 0) stats.residual_fdr_at_one = fdr;
            if (!stats.mig_size_threshold && fdr <= request.target_fdr) {
                stats.mig_size_threshold = static_cast<uint32_t>(1ull << b);
            }
        }

        std::FILE* fh =
            std::fopen((out_dir / (stats.sample_id + ".bins.tsv")).string().c_str(), "w");
        if (!fh) throw MigecError("refine: cannot write the bin table");
        std::fprintf(fh, "min_reads\tmax_reads\tbarcodes\treads\tmerged\tfraction_erroneous\t"
                         "molecules\tsuspected_residual\tresidual_fdr\tpayload_entropy_bits\n");
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
            const uint64_t surviving = surviving_in[b];
            std::fprintf(fh, "%llu\t%llu\t%llu\t%llu\t%llu\t%.6f\t%llu\t%llu\t%.6f\t%.4f\n",
                         1ull << b, (1ull << (b + 1)) - 1,
                         static_cast<unsigned long long>(barcodes[b]),
                         static_cast<unsigned long long>(reads_in[b]),
                         static_cast<unsigned long long>(merged_in[b]),
                         static_cast<double>(merged_in[b]) / static_cast<double>(barcodes[b]),
                         static_cast<unsigned long long>(surviving),
                         static_cast<unsigned long long>(suspected[b]),
                         surviving ? static_cast<double>(suspected[b]) /
                                         static_cast<double>(surviving) : 0.0,
                         scored ? entropy / scored : 0.0);
        }
        std::fclose(fh);
    }

    stats.wall_seconds = clock.seconds();
    return stats;
}

}  // namespace migec
