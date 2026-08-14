#include "migec/refine.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <filesystem>
#include <map>
#include <unordered_map>

#include "migec/fastq.hpp"
#include "migec/mig_record.hpp"
#include "migec/parallel.hpp"
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

// One read, as refine's passes use it, from either input.
//
// The strings are views: into the FASTQ record on the FASTQ route, and into the walker's own
// scratch on the `.mig` route -- valid until the next read either way. Keeping the passes on
// STRINGS rather than on packed keys is deliberate: it is the same code for both routes, so a
// `.mig` run cannot drift from a FASTQ run by packing the key a second, subtly different way.
struct SourceRead {
    std::string_view umi, cell;   // the barcodes
    std::string_view qx, cy;      // ...and their own quality, empty when the input has none
    std::string_view seq, qual;
    std::string_view name, comment;  // FASTQ only; a `.mig` record carries no read name
    uint64_t src_index = 0;
    uint16_t flags = 0;
};

// Walks the input once, whichever it is. `fn` returning false stops the walk, which is how the
// intake limit ends a pass.
void for_each_read(const RefineRequest& request, const std::function<bool(const SourceRead&)>& fn) {
    if (request.mig_inputs.empty()) {
        FastqReader reader(request.input);
        FastqRecord rec;
        SourceRead sr;
        while (reader.next(rec)) {
            sr.umi = tag_value(rec.comment, "RX:Z:");
            sr.cell = tag_value(rec.comment, "CB:Z:");
            sr.qx = tag_value(rec.comment, "QX:Z:");
            sr.cy = tag_value(rec.comment, "CY:Z:");
            sr.seq = rec.seq;
            sr.qual = rec.qual;
            sr.name = rec.name;
            sr.comment = rec.comment;
            if (!fn(sr)) return;
        }
        return;
    }
    // Bucket order is key order, and the buckets were handed over sorted, so this walk is stable
    // across passes -- which is what lets pass 3 line up with the table pass 1 built.
    std::string umi, cell;  // scratch, assigned into rather than rebuilt
    SourceRead sr;
    for (const std::string& path : request.mig_inputs) {
        MigReader reader(path);
        const MigHeader& h = reader.header();
        MigRecord rec;
        while (reader.next(rec)) {
            umi = unpack_barcode(rec.umi, h.umi_len);
            cell = h.cell_len ? unpack_barcode(rec.cell, h.cell_len) : std::string();
            sr.umi = umi;
            sr.cell = cell;
            // Empty on a v1 file, which is exactly how a FASTQ with no QX tag reads: the
            // posterior falls back to the library's global rate rather than to a worse number.
            sr.qx = rec.qual_umi;
            sr.cy = rec.qual_cell;
            sr.seq = rec.seq1;
            sr.qual = rec.qual1;
            sr.name = {};
            sr.comment = {};
            sr.src_index = rec.src_index;
            sr.flags = rec.flags;
            if (!fn(sr)) return;
        }
    }
}

}  // namespace

RefineStats refine(const RefineRequest& request) {
    Stopwatch clock;
    RefineStats stats;
    stats.sample_id = request.sample_id;

    std::filesystem::path out_dir(request.output_dir);
    std::filesystem::create_directories(out_dir);

    // `.mig` input: read the partition's shape off the files rather than choosing one. The
    // buckets ARE the partition and the output keeps it, so `assemble` can take refine's output
    // exactly as it takes checkout's.
    int mig_bits = 0;
    uint64_t written_reads = 0;
    if (!request.mig_inputs.empty()) {
        bool first = true;
        for (const std::string& path : request.mig_inputs) {
            MigReader probe(path);
            const MigHeader& h = probe.header();
            if (first) {
                mig_bits = h.bucket_bits;
                if (stats.sample_id.empty()) stats.sample_id = h.sample_id;
                first = false;
            } else if (h.bucket_bits != mig_bits) {
                throw MigecError("refine: '" + path + "' is cut into 2^" +
                                 std::to_string(h.bucket_bits) + " buckets and the first file into 2^" +
                                 std::to_string(mig_bits) + " -- these are two partitions, not one");
            } else if (h.sample_id != stats.sample_id) {
                // Never: two samples in one table would correct one sample's barcode into the
                // other's. A UMI repeats across samples by design.
                throw MigecError("refine: '" + path + "' holds sample '" + h.sample_id +
                                 "' but the run is refining '" + stats.sample_id +
                                 "' -- refine is a per-sample stage");
            }
        }
    }

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
        for_each_read(request, [&](const SourceRead& r) {
            if (r.cell.empty()) return true;
            if (static_cast<int>(r.cell.size()) != list.length()) {
                throw MigecError(
                    "refine: the reads carry a " + std::to_string(r.cell.size()) +
                    " nt cell barcode and the whitelist holds " +
                    std::to_string(list.length()) + " nt entries");
            }
            const uint64_t key = pack_barcode(r.cell);
            ++seen[key];
            ++total_reads;
            if (quals.find(key) == quals.end()) quals.emplace(key, std::string(r.cy));
            return true;
        });
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

    Stopwatch phase;
    // ---------------------------------------------------------------- pass 1: the barcode table
    // Counts only. The entry array is not final until every read has been seen, and the evidence
    // is indexed in parallel with it, so it cannot be filled in the same pass.
    UmiCounts counts(0);
    {
        bool started = false;
        IntakeLimit limit = request.limit;
        for_each_read(request, [&](const SourceRead& r) {
            ++stats.reads;
            const std::string_view umi = r.umi;
            if (umi.empty()) { ++stats.reads_without_umi; return true; }
            // After the barcode is known, so --limit-umi counts barcodes. Pass 3 stops at the same
            // read count, so the rewritten file is exactly the prefix the table was built from --
            // rewriting reads whose barcode was never in the table would emit uncorrected reads
            // beside corrected ones with nothing saying which was which.
            if (limit.active() && !limit.admit(stats.reads, pack_barcode(umi))) {
                --stats.reads;
                stats.limited = true;
                return false;
            }
            const std::string_view cell = r.cell;
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
                    const std::string_view bc = tag_value(r.comment, "BC:Z:");
                    stats.sample_id = bc.empty() ? std::string("sample") : std::string(bc);
                }
            } else if (static_cast<int>(umi.size()) != stats.umi_length) {
                throw MigecError("refine: read '" + std::string(r.name) + "' carries a " +
                                 std::to_string(umi.size()) + " nt UMI where the file started "
                                 "with " + std::to_string(stats.umi_length) +
                                 " -- these are two runs concatenated, not one sample");
            }
            counts.add(cell.empty()
                           ? pack_barcode(umi)
                           : pack_key_snapped(snap(pack_barcode(cell)), stats.cell_length, umi));
            return true;
        });
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
        for_each_read(request, [&](const SourceRead& r) {
            const std::string_view umi = r.umi;
            if (umi.empty()) return true;
            const std::string_view cell = r.cell;
            const size_t i = slot(cell.empty()
                                     ? pack_barcode(umi)
                                     : pack_key_snapped(snap(pack_barcode(cell)),
                                                        stats.cell_length, umi));
            if (i >= entries.size()) return true;
            if (request.use_quality) {
                // The whole key's quality, cell part first, to line up with the packed key.
                std::string qx(r.cy);
                qx.append(r.qx);
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
                const size_t n = std::min<size_t>(static_cast<size_t>(pw), r.seq.size());
                std::copy_n(r.seq.begin(), n,
                            evidence.payload.begin() +
                                static_cast<std::ptrdiff_t>(i * static_cast<size_t>(pw)));
            }
            return true;
        });
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

    stats.table_seconds = phase.seconds();

    // ---------------------------------------------------------------- correct
    phase = Stopwatch();
    stats.threads = worker_count(request.correction.threads, counts.entries().size());
    const CorrectionResult correction = correct_umis(counts, request.correction, evidence);
    stats.correct_seconds = phase.seconds();
    phase = Stopwatch();
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
    //
    // Chunked and parallel, the same shape checkout uses: a round of chunks is read serially,
    // rewritten AND compressed on the workers, and appended in chunk order. Never: the chunk size is
    // a constant, so the member boundaries -- and therefore the bytes -- do not depend on -t.
    const std::filesystem::path fastq_path =
        out_dir / (stats.sample_id + ".fq" + (request.gzip_level > 0 ? ".gz" : ""));
    if (!request.mig_inputs.empty()) {
        // `.mig` in, `.mig` out. The corrected barcode is a different key, and a key decides its
        // bucket, so a corrected read can belong in a bucket other than the one it arrived in --
        // the output is re-partitioned on the NEW key rather than copied bucket for bucket, or it
        // would not be a partition any more and `assemble` would group across it.
        //
        // Note: the audit trail is `<sample>.barcodes.tsv`, which carries every barcode with its
        // parent. A `.mig` record has no room for the pre-correction barcode the way a FASTQ
        // comment has OX:Z:, and adding two u64 columns to every record to carry it would cost
        // 16 bytes a READ for something that is one row a BARCODE in the table.
        //
        // Note: serial, and deliberately so for now. The rewrite is a decompress, a table lookup
        // and a recompress per record; the FASTQ route threads it because it was measured at 83%
        // of that run's wall clock. This one is measured in tests/benchmark/test_refine_speed.py
        // and threads when a number says to, not before.
        const size_t n_buckets = static_cast<size_t>(1) << mig_bits;
        std::vector<std::unique_ptr<MigWriter>> writers(n_buckets);
        std::vector<std::string> paths(n_buckets);
        const size_t block_bytes = std::clamp<size_t>(32u << 20 >> mig_bits, 256u << 10, 4u << 20);
        std::string whole;
        for (const std::string& path : request.mig_inputs) {
            MigReader reader(path);
            const MigHeader& h = reader.header();
            MigRecord rec;
            while (reader.next(rec)) {
                std::string umi = unpack_barcode(rec.umi, h.umi_len);
                std::string cell = h.cell_len ? unpack_barcode(rec.cell, h.cell_len) : std::string();
                const uint64_t snapped = cell.empty() ? 0 : snap(pack_barcode(cell));
                const size_t i = slot(cell.empty()
                                          ? pack_barcode(umi)
                                          : pack_key_snapped(snapped, stats.cell_length, umi));
                MigRecord out = rec;
                if (i < entries.size()) {
                    whole = unpack_barcode(entries[correction.root[i]].key, L);
                    const std::string new_cell = whole.substr(0, static_cast<size_t>(stats.cell_length));
                    const std::string new_umi = whole.substr(static_cast<size_t>(stats.cell_length));
                    if (new_umi != umi) {
                        out.umi = pack_barcode(new_umi);
                        out.flags = static_cast<uint16_t>(out.flags | kUmiCorrected);
                    }
                    if (!cell.empty() && new_cell != cell) {
                        out.cell = pack_barcode(new_cell);
                        out.flags = static_cast<uint16_t>(out.flags | kCellCorrected);
                    }
                }
                const uint32_t b =
                    bucket_of(h.cell_len ? out.cell : out.umi, static_cast<int>(mig_bits));
                if (!writers[b]) {
                    MigHeader oh;
                    oh.umi_len = static_cast<uint8_t>(stats.umi_length);
                    oh.cell_len = static_cast<uint8_t>(stats.cell_length);
                    oh.bucket_index = static_cast<uint8_t>(b);
                    oh.bucket_bits = static_cast<uint8_t>(mig_bits);
                    oh.paired = h.paired;
                    // The quality stored is the READ's own, of the bases it actually carried --
                    // not the parent's. kUmiCorrected says the barcode was replaced; the quality
                    // describes what this read saw, which is the only thing it can describe.
                    oh.barcode_quality = h.barcode_quality;
                    oh.sample_id = stats.sample_id;
                    paths[b] = (out_dir / (stats.sample_id + "." + bucket_suffix(b) + ".mig"))
                                   .string();
                    writers[b] = std::make_unique<MigWriter>(paths[b], oh, block_bytes);
                }
                writers[b]->write(out);
                ++written_reads;
            }
        }
        for (size_t b = 0; b < n_buckets; ++b) {
            if (!writers[b]) continue;
            writers[b]->close();
            stats.mig_paths.push_back(paths[b]);
        }
    } else {
        // Reads per chunk. A CONSTANT, so the gzip member boundaries -- and therefore the bytes
        // -- do not depend on --threads. It is deliberately small: the rewrite holds one chunk of
        // records plus one output buffer PER WORKER, so this number is multiplied by the thread
        // count. 4096 x ~250 B is ~1 MB a thread each way; 16384 was four times that and showed
        // up as 100 MB of fixed cost at 32 threads before anything had been read.
        constexpr size_t kChunkReads = 4096;
        const int rewrite_threads = worker_count(request.correction.threads, 64);
        // Sized once and ASSIGNED into, never cleared -- the same trade assemble's partition
        // makes. `clear()` destroys four std::string per record, so every round would cost
        // 4 x kChunkReads allocations and the reader would sit in malloc instead of inflate.
        std::vector<std::vector<FastqOwned>> chunks(static_cast<size_t>(rewrite_threads),
                                                    std::vector<FastqOwned>(kChunkReads));
        std::vector<size_t> held(static_cast<size_t>(rewrite_threads), 0);
        std::vector<std::string> plain(static_cast<size_t>(rewrite_threads));
        std::vector<std::string> packed(static_cast<size_t>(rewrite_threads));
        std::vector<std::string> scratch_comment(static_cast<size_t>(rewrite_threads));
        std::vector<std::string> scratch_whole(static_cast<size_t>(rewrite_threads));

        std::FILE* fh = std::fopen(fastq_path.string().c_str(), "wb");
        if (!fh) throw MigecError("refine: cannot write " + fastq_path.string());
        FastqReader reader(request.input);
        FastqRecord rec;
        bool eof = false;
        uint64_t written = 0;

        auto rewrite_one = [&](const FastqOwned& r, std::string& dst, std::string& comment,
                               std::string& whole) {
            const std::string_view umi = tag_value(r.comment, "RX:Z:");
            if (umi.empty()) return;
            const std::string_view cell = tag_value(r.comment, "CB:Z:");
            const uint64_t snapped = cell.empty() ? 0 : snap(pack_barcode(cell));
            const size_t i = slot(cell.empty()
                                     ? pack_barcode(umi)
                                     : pack_key_snapped(snapped, stats.cell_length, umi));
            if (i >= entries.size()) {
                // Unreachable unless pass 1 and pass 3 disagree about the key. Pass the read
                // through untouched rather than dropping it: emitting a read whose barcode was
                // not corrected is recoverable, losing it silently is not.
                append_fastq(dst, r.name, r.comment, r.seq, r.qual);
                return;
            }
            const uint32_t root = correction.root[i];
            comment.assign(r.comment);
            // The final barcode is the root's, which already carries the whitelist snap because
            // the snap happened before the table was built. Compare against what the READ said,
            // not against the root only: a barcode can be snapped by the whitelist without being
            // merged by the posterior, and rewriting it only in the second case would leave the
            // CB tag disagreeing with the key everything else was grouped on.
            //
            // `comment` and `whole` are the worker's own scratch, reused read after read, and the
            // two halves are views into `whole`: a fresh string plus two substr is four
            // allocations per read on a path every read takes.
            whole = unpack_barcode(entries[root].key, L);
            const std::string_view new_cell =
                std::string_view(whole).substr(0, static_cast<size_t>(stats.cell_length));
            const std::string_view new_umi =
                std::string_view(whole).substr(static_cast<size_t>(stats.cell_length));
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
            append_fastq(dst, r.name, comment, r.seq, r.qual);
        };

        while (!eof) {
            size_t filled = 0;
            for (; filled < chunks.size(); ++filled) {
                std::vector<FastqOwned>& c = chunks[filled];
                size_t n = 0;
                while (n < kChunkReads && (!stats.limited || written < stats.reads) &&
                       reader.next(rec)) {
                    ++written;
                    FastqOwned& slot = c[n++];
                    slot.name.assign(rec.name);
                    slot.comment.assign(rec.comment);
                    slot.seq.assign(rec.seq);
                    slot.qual.assign(rec.qual);
                }
                held[filled] = n;
                if (n == 0) { eof = true; break; }
                if (n < kChunkReads) { ++filled; eof = true; break; }
            }
            if (filled == 0) break;

            parallel_for(filled, rewrite_threads, [&](size_t t, int) {
                std::string& dst = plain[t];
                dst.clear();
                for (size_t k = 0; k < held[t]; ++k) {
                    rewrite_one(chunks[t][k], dst, scratch_comment[t], scratch_whole[t]);
                }
                if (request.gzip_level > 0) {
                    gzip_member(dst, packed[t], request.gzip_level);
                } else {
                    packed[t] = dst;
                }
            });

            for (size_t t = 0; t < filled; ++t) {
                if (packed[t].empty()) continue;
                if (std::fwrite(packed[t].data(), 1, packed[t].size(), fh) != packed[t].size()) {
                    std::fclose(fh);
                    throw MigecError("refine: short write");
                }
            }
        }
        std::fclose(fh);
    }
    stats.rewrite_seconds = phase.seconds();

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

        // Cell Ranger's barcode-rank plot, which is the figure everyone already knows how to
        // read: cells sorted by how many DISTINCT UMIs they carry, on log-log, with the knee
        // where real cells stop and ambient starts. Never reads -- an over-amplified molecule
        // would put an empty droplet high up the curve, which is exactly the artefact the plot
        // exists to make visible.
        //
        // Log-spaced ranks, because one row per barcode is hundreds of millions of rows for a
        // figure that is read on a log axis anyway. The first and last rank are always emitted so
        // the curve's ends are exact.
        if (!sizes.empty()) {
            std::FILE* rf =
                std::fopen((out_dir / (stats.sample_id + ".cell_rank.tsv")).string().c_str(), "w");
            if (!rf) throw MigecError("refine: cannot write the cell rank table");
            std::fprintf(rf, "rank\tumis\tcalled\tcumulative_umis\tcumulative_fraction\n");
            uint64_t total = 0;
            for (uint32_t v : sizes) total += v;
            uint64_t cum = 0;
            size_t next = 0;
            for (size_t i = 0; i < sizes.size(); ++i) {
                cum += sizes[i];
                if (i == next || i + 1 == sizes.size()) {
                    std::fprintf(rf, "%zu\t%u\t%d\t%llu\t%.6f\n", i + 1, sizes[i],
                                 sizes[i] >= stats.cell_threshold ? 1 : 0,
                                 static_cast<unsigned long long>(cum),
                                 total ? static_cast<double>(cum) / static_cast<double>(total)
                                       : 0.0);
                    next = std::max(next + 1,
                                    static_cast<size_t>(static_cast<double>(next) * 1.05));
                }
            }
            std::fclose(rf);
        }
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
        // The EXACT MIG size spectrum: one row per distinct size, not a power-of-two bin. Two
        // things need it and neither can be had from the binned table:
        //
        //   * the size histogram with BOTH series -- how many molecules were seen n times, and how
        //     many reads those molecules account for. The two peak in different places whenever
        //     the library is over-sequenced, and a bin four wide hides that;
        //   * the rank/Zipf curve, which is the cumulative count of this table. Power-of-two bins
        //     turn it into four steps.
        //
        // It costs one row per distinct depth, which is bounded by the deepest molecule and is a
        // few thousand rows on any real library -- not one row per molecule.
        std::map<uint32_t, std::pair<uint64_t, uint64_t>> spectrum;  // size -> (molecules, reads)
        for (uint32_t c : correction.corrected) {
            if (c == 0) continue;
            auto& e = spectrum[c];
            ++e.first;
            e.second += c;
        }
        std::FILE* sf =
            std::fopen((out_dir / (stats.sample_id + ".sizes.tsv")).string().c_str(), "w");
        if (!sf) throw MigecError("refine: cannot write the size spectrum");
        std::fprintf(sf, "size\tlog1p_size\tmolecules\treads\n");
        for (const auto& [size, counts] : spectrum) {
            std::fprintf(sf, "%u\t%.6f\t%llu\t%llu\n", size,
                         std::log1p(static_cast<double>(size)),
                         static_cast<unsigned long long>(counts.first),
                         static_cast<unsigned long long>(counts.second));
        }
        std::fclose(sf);
    }
    {
        // ------------------------------------------- barcode errors against the parent's depth
        //
        // The error rate measured at every amplification depth, instead of once for the library.
        // A parent carrying c reads had c*L barcode bases for the instrument and the polymerase to
        // get wrong, and its error children are what they got wrong, so two estimators of the same
        // eps fall out of the same row:
        //
        //   distinct children   u(c) = 3L (1 - exp(-c eps / 3))   -> eps = -(3/c) ln(1 - u/3L)
        //   reads in children   r(c) = c L eps                    -> eps = r / (c L)
        //
        // They are the check on each other, not two opinions. The first counts distinct
        // NEIGHBOURS and a barcode has only 3L of them, so it saturates; the second counts READS,
        // of which there is no ceiling, so it does not. Where the two part company on the figure
        // is where the neighbourhood filled, read off the data rather than predicted. Past
        // saturation the first is left blank rather than reported as a small number, because
        // inverting a full neighbourhood returns "no errors" for the most error-ridden library
        // there can be.
        //
        // Never: BOTH are bounded by the merges correction actually made, so neither is a
        // saturation-free estimator and this table must not be read as one. Measured against an
        // injected rate on simulated libraries (`tests/synthetic/test_umi_errors.py`), as a
        // fraction of the truth:
        //
        //     occupancy    0.2%   2.3%   9.8%    33%    100%
        //     distance-1   0.97   0.96   0.76   0.45   0.001
        //     children     0.99   0.95   0.88   0.62   0.00
        //
        // So it is the better of the two wherever either works, and at 100% they both go to zero
        // for the same reason: on a full barcode space `correct_umis` refuses to merge -- and it
        // is right to, because a distance-1 neighbour there is more likely a real molecule than a
        // child. `saturated` is what says the answer is a floor. Read that flag; do not read this
        // table instead of it.
        //
        // Never: this counts children that were FOUND, so at 1-3 reads/UMI it is a lower bound --
        // 80% of barcode errors there have no sequenced parent to be merged into and are
        // unreachable in principle. That is why `error_at_depth` is read where correction is
        // near-complete instead of averaged over every molecule, and why the depth it was read at
        // travels with it.
        struct Row { uint64_t parents = 0, children = 0, child_reads = 0; };
        std::map<uint32_t, Row> by_depth;  // the PARENT's own reads -> what its children came to
        for (size_t i = 0; i < entries.size(); ++i) {
            const uint32_t root = correction.root[i];
            if (root == i) {
                ++by_depth[entries[i].count].parents;
            } else {
                Row& r = by_depth[entries[root].count];
                ++r.children;
                r.child_reads += entries[i].count;
            }
        }

        const double Lf = static_cast<double>(L);
        const double ceiling = 3.0 * Lf;  // distinct distance-1 neighbours a barcode can have
        // Deep enough that a child almost always has its parent in the table, so the count is not
        // a lower bound worth quoting. Below this the estimate is reported per row but not taken.
        constexpr uint32_t kTrustedDepth = 10;
        uint64_t all_parent_reads = 0, all_child_reads = 0;
        uint64_t deep_parent_reads = 0, deep_child_reads = 0;
        for (const auto& [depth, row] : by_depth) {
            all_parent_reads += static_cast<uint64_t>(depth) * row.parents;
            all_child_reads += row.child_reads;
            if (depth >= kTrustedDepth) {
                deep_parent_reads += static_cast<uint64_t>(depth) * row.parents;
                deep_child_reads += row.child_reads;
            }
        }
        stats.error_from_children =
            all_parent_reads ? static_cast<double>(all_child_reads) /
                                   (static_cast<double>(all_parent_reads) * Lf) : 0.0;
        if (deep_parent_reads) {
            stats.error_at_depth = static_cast<double>(deep_child_reads) /
                                   (static_cast<double>(deep_parent_reads) * Lf);
            stats.error_depth = kTrustedDepth;
            if (stats.error_at_depth > 0.0) {
                stats.error_phred = -10.0 * std::log10(stats.error_at_depth);
            }
        }

        // `neighbours` and `estimate` are constant down the column, and that is deliberate: the
        // figure draws the saturation ceiling and the library's own estimate as reference lines,
        // and a figure that cannot be redrawn from the committed table alone is a figure that will
        // one day disagree with the report. Two repeated columns are cheaper than a gnuplot script
        // that shells out to recover them.
        std::FILE* ef =
            std::fopen((out_dir / (stats.sample_id + ".umi_errors.tsv")).string().c_str(), "w");
        if (!ef) throw MigecError("refine: cannot write the barcode error table");
        std::fprintf(ef, "parent_reads\tparents\tchild_barcodes\tchild_reads\t"
                         "children_per_parent\treads_per_parent\tneighbours\t"
                         "error_from_variants\terror_from_reads\tphred_from_reads\testimate\n");
        for (const auto& [depth, row] : by_depth) {
            if (!row.parents) continue;
            const double c = static_cast<double>(depth);
            const double np = static_cast<double>(row.parents);
            const double u = static_cast<double>(row.children) / np;
            const double r = static_cast<double>(row.child_reads) / np;
            const double eps_reads = r / (c * Lf);
            const bool invertible = u < ceiling;
            const double eps_variants = invertible ? -(3.0 / c) * std::log1p(-u / ceiling) : 0.0;
            std::fprintf(ef, "%u\t%llu\t%llu\t%llu\t%.6f\t%.6f\t%.1f\t", depth,
                         static_cast<unsigned long long>(row.parents),
                         static_cast<unsigned long long>(row.children),
                         static_cast<unsigned long long>(row.child_reads), u, r, ceiling);
            if (invertible) {
                std::fprintf(ef, "%.6e\t", eps_variants);
            } else {
                std::fprintf(ef, ".\t");  // the neighbourhood is full; inverting it says nothing
            }
            if (eps_reads > 0.0) {
                std::fprintf(ef, "%.6e\t%.2f\t", eps_reads, -10.0 * std::log10(eps_reads));
            } else {
                std::fprintf(ef, "0.000000e+00\t.\t");
            }
            std::fprintf(ef, "%.6e\n", stats.error_at_depth);
        }
        std::fclose(ef);
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
        // Note: "Much larger neighbour" alone is not the test. At 1-3 reads per UMI nothing is 20x
        // anything, so a count-ratio criterion reports zero residual in precisely the regime where
        // the residual is worst -- the same trap the correction posterior itself fell into. The
        // payload is what still separates them at one read: a neighbour whose reads agree on the
        // molecule is a child whatever the counts say.
        std::vector<uint64_t> suspected(nbins, 0);
        {
            const int width = L;
            // Threaded, and it is the whole reason this is not the serial tail any more: the scan
            // is 3L binary searches per surviving barcode, which measured 0.53 s of a 2.17 s run on
            // 2 M reads -- a quarter of refine, on one core, after everything else had been
            // parallelised. It reads the barcode table and the payload draft and writes nothing
            // shared, so a per-worker bin counter merged afterwards is the whole change; the
            // counters are integers, so the merge is order-independent and `-t` still changes
            // nothing but the clock.
            std::vector<std::vector<uint64_t>> per_worker(
                static_cast<size_t>(stats.threads), std::vector<uint64_t>(nbins, 0));
            parallel_for(entries.size(), stats.threads, [&](size_t i, int w) {
                if (correction.corrected[i] == 0) return;
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
                if (!looks_like_a_child) return;
                size_t bidx = 0;
                while ((mine >> bidx) > 1) ++bidx;
                if (bidx >= nbins) bidx = nbins - 1;
                ++per_worker[static_cast<size_t>(w)][bidx];
            });
            for (const std::vector<uint64_t>& mine : per_worker) {
                for (size_t b = 0; b < nbins; ++b) {
                    suspected[b] += mine[b];
                    stats.suspected_residual += mine[b];
                }
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
