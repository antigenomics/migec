// Python bindings. Deliberately thin: everything here is either a whole-file operation (so the
// per-call overhead is irrelevant) or a small utility used by tests and notebooks. Nothing in a
// per-read loop crosses this boundary.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include "migec/assemble.hpp"
#include "migec/checkout.hpp"
#include "migec/fastq.hpp"
#include "migec/mig_record.hpp"
#include "migec/pattern.hpp"
#include "migec/refine.hpp"
#include "migec/resource.hpp"
#include "migec/subsample.hpp"
#include "migec/suggest.hpp"
#include "migec/types.hpp"
#include "migec/umi_stats.hpp"
#include "migec/version.hpp"

namespace py = pybind11;
using namespace migec;

namespace {

// One decoded record, owning its strings, for the Python side. Views into the reader arena would
// dangle the moment the iterator advanced.
struct PyMigRecord {
    std::string cell, umi;
    uint64_t src_index;
    uint16_t flags;
    uint8_t umi_minq, cell_minq;
    std::string seq1, qual1, seq2, qual2;
};

class PyMigFile {
public:
    explicit PyMigFile(const std::string& path) : reader_(path) {}

    py::dict header() const {
        const MigHeader& h = reader_.header();
        py::dict d;
        d["format_version"] = h.format_version;
        d["umi_len"] = h.umi_len;
        d["cell_len"] = h.cell_len;
        d["bucket_index"] = h.bucket_index;
        d["bucket_bits"] = h.bucket_bits;
        d["paired"] = h.paired;
        d["sample_id"] = h.sample_id;
        d["provenance"] = h.provenance;
        d["quality_calibration"] = h.quality_calibration;
        return d;
    }

    std::vector<PyMigRecord> read_all() {
        const MigHeader& h = reader_.header();
        std::vector<PyMigRecord> out;
        MigRecord r;
        while (reader_.next(r)) {
            PyMigRecord p;
            p.cell = unpack_barcode(r.cell, h.cell_len);
            p.umi = unpack_barcode(r.umi, h.umi_len);
            p.src_index = r.src_index;
            p.flags = r.flags;
            p.umi_minq = r.umi_minq;
            p.cell_minq = r.cell_minq;
            p.seq1 = std::string(r.seq1);
            p.qual1 = std::string(r.qual1);
            p.seq2 = std::string(r.seq2);
            p.qual2 = std::string(r.qual2);
            out.push_back(std::move(p));
        }
        return out;
    }

private:
    MigReader reader_;
};

// Runs a whole FASTQ (or a pair) through checkout in C++, writing one output per sample.
// Whole-file, so the pybind boundary is crossed once per run rather than once per read.
py::dict py_run_checkout(const std::string& in_path, const std::string& in_path2,
                         const std::vector<std::string>& sample_ids,
                         const std::vector<std::string>& patterns,
                          const std::vector<std::string>& slaves, const std::string& out_prefix,
                         const std::string& trim, int min_umi_quality, bool write_unmatched,
                         int threads, int max_offset, uint64_t limit_reads) {
    if (sample_ids.size() != patterns.size()) {
        throw MigecError("run_checkout: sample_ids and patterns have different lengths");
    }
    PatternSet set;
    for (size_t i = 0; i < patterns.size(); ++i) {
        set.add(sample_ids[i], patterns[i], i < slaves.size() ? slaves[i] : std::string());
    }

    CheckoutParams params;
    params.match.max_offset = max_offset;
    params.min_umi_quality = min_umi_quality;
    if (trim == "none") {
        params.trim = TrimMode::kNone;
    } else if (trim == "pattern") {
        params.trim = TrimMode::kPattern;
    } else {
        throw MigecError("run_checkout: trim must be 'none' or 'pattern', got '" + trim + "'");
    }

    CheckoutRequest req;
    req.r1 = in_path;
    req.r2 = in_path2;
    req.out_prefix = out_prefix;
    req.write_unmatched = write_unmatched;
    req.threads = threads;
    req.limit_reads = limit_reads;

    // The clock covers the per-sample statistics below as well as the driver. They are serial,
    // single-threaded, and on a 12 nt UMI cost about 2 us per read -- four times the matching. A
    // reads/s figure that stopped at the driver would be measuring the matcher, not checkout.
    Stopwatch clock;
    CheckoutStats stats;
    {
        py::gil_scoped_release release;
        stats = run_checkout(set, params, req);
    }

    const CheckoutCounters& c = stats.counters;
    py::dict out;
    out["total"] = c.total;
    out["assigned"] = c.assigned;
    out["unmatched"] = c.unmatched;
    out["ambiguous"] = c.ambiguous;
    out["short_payload"] = c.short_payload;
    out["bad_umi"] = c.bad_umi;
    out["normalised"] = c.normalised;
    out["paired"] = !in_path2.empty();
    out["threads"] = stats.threads;
    {
        const QualityCalibration& c = stats.counters.calibration;
        py::dict cal;
        cal["fitted"] = c.fitted;
        cal["quality_independent"] = c.quality_independent;
        cal["slope"] = c.slope;
        cal["bases"] = c.bases;
        cal["positions_dropped"] = c.positions_dropped;
        py::list pos;
        for (size_t i = 0; i < c.by_position.size(); ++i) {
            uint64_t n = 0, bad = 0;
            for (size_t q = 0; q < 61; ++q) {
                n += c.by_position[i][q][0] + c.by_position[i][q][1];
                bad += c.by_position[i][q][1];
            }
            if (!n) continue;
            py::dict e;
            e["position"] = i;
            e["bases"] = n;
            e["mismatches"] = bad;
            e["rate"] = static_cast<double>(bad) / static_cast<double>(n);
            e["used"] = i < c.position_used.size() && c.position_used[i] != 0;
            pos.append(e);
        }
        cal["per_position"] = pos;
        py::list per_q;
        for (size_t q = 0; q < c.counts.size(); ++q) {
            const uint64_t n = c.counts[q][0] + c.counts[q][1];
            if (!n) continue;
            py::dict e;
            e["phred"] = q;
            e["bases"] = n;
            e["mismatches"] = c.counts[q][1];
            e["observed"] = static_cast<double>(c.counts[q][1]) / static_cast<double>(n);
            e["nominal"] = phred_error(static_cast<uint8_t>(q));
            e["fitted"] = c.error(static_cast<int>(q));
            per_q.append(e);
        }
        cal["per_phred"] = per_q;
        out["quality_calibration"] = cal;
    }
    out["match_seconds"] = stats.wall_seconds;
    out["trimmed_bases"] = c.trimmed_bases;
    out["peak_rss_bytes"] = stats.peak_rss_bytes;
    out["umi_memory_bytes"] = stats.umi_memory_bytes;

    const std::vector<UmiCounts>& umi_counts = stats.umi_counts;
    py::list per_sample;
    // Distinct sample ids, not input rows: two rows may declare the same sample.
    for (size_t i = 0; i < stats.sample_ids.size(); ++i) {
        const CoverageHistogram h = umi_counts[i].histogram();
        const UmiComposition comp = umi_counts[i].composition(false);
        py::dict s;
        s["sample_id"] = stats.sample_ids[i];
        s["reads"] = stats.sample_reads[i];
        s["umis"] = umi_counts[i].distinct();
        s["mean_reads_per_umi"] = h.mean_reads_per_umi();
        s["reads_in_migs_ge5"] = h.reads_in_migs_at_least(5);
        s["over_sequenced"] = h.over_sequenced();
        s["hist_reads"] = h.reads;
        s["hist_units"] = h.units;
        s["umi_length"] = comp.length;
        s["total_entropy"] = comp.total_entropy();
        s["total_information"] = comp.total_information();
        s["effective_length"] = comp.effective_length();
        s["effective_space"] = comp.effective_space();
        py::list per_pos;
        for (int j = 0; j < comp.length; ++j) {
            py::dict d;
            d["position"] = j;
            d["A"] = comp.freq[static_cast<size_t>(j)][0];
            d["C"] = comp.freq[static_cast<size_t>(j)][1];
            d["G"] = comp.freq[static_cast<size_t>(j)][2];
            d["T"] = comp.freq[static_cast<size_t>(j)][3];
            d["entropy"] = comp.entropy(j);
            d["information"] = comp.information(j);
            d["collision"] = comp.collision(j);
            per_pos.append(d);
        }
        s["composition"] = per_pos;

        // The payload length distribution after trimming, and the mean. A pattern placed one base
        // off trims one base off every read, and this is the only place that shows.
        py::list lens;
        double len_sum = 0.0, len_n = 0.0;
        if (i < stats.sample_payload_len.size()) {
            for (size_t L = 0; L < kPayloadHistLen; ++L) {
                const uint64_t n = stats.sample_payload_len[i][L];
                if (!n) continue;
                py::dict e;
                e["length"] = L;
                e["reads"] = n;
                e["at_least"] = L + 1 == kPayloadHistLen;  // the catch-all bin
                lens.append(e);
                len_sum += static_cast<double>(L * n);
                len_n += static_cast<double>(n);
            }
        }
        s["payload_lengths"] = lens;
        s["mean_payload_length"] = len_n ? len_sum / len_n : 0.0;

        const CorrectionResult cr = correct_umis(umi_counts[i]);
        s["umi_error_rate"] = cr.estimated_error;

        // The barcode space, as the birthday problem sees it: how full it is, how many molecules
        // that implies, and how many MIGs are really two molecules pooled.
        const BarcodeSpace bs = barcode_space(comp, umi_counts[i].distinct());
        py::dict space;
        space["length"] = bs.length;
        space["nominal_space"] = bs.nominal_space;
        space["effective_space"] = bs.effective_space;
        space["effective_length"] = bs.effective_length;
        space["bias_loss"] = bs.bias_loss;
        space["observed"] = bs.observed;
        space["occupancy"] = bs.occupancy;
        space["lambda"] = bs.lambda;
        space["molecules"] = bs.molecules;
        space["hidden"] = bs.hidden;
        space["p_multi"] = bs.p_multi;
        space["saturated"] = bs.saturated;
        s["barcode_space"] = space;

        // ...and the error budget: what the reported Phred and the polymerase predict, against
        // what the distance-1 excess actually found.
        const ErrorBudget eb =
            error_budget(comp, stats.sample_phred[i], cr.estimated_error, umi_counts[i].distinct());
        py::dict budget;
        budget["from_phred"] = eb.from_phred;
        budget["mean_phred"] = eb.mean_phred;
        budget["from_polymerase"] = eb.from_polymerase;
        budget["predicted"] = eb.predicted;
        budget["estimated"] = eb.estimated;
        budget["ratio"] = eb.ratio;
        budget["barcodes_with_error"] = eb.barcodes_with_error;
        budget["neighbour_occupancy"] = eb.neighbour_occupancy;
        budget["estimate_unreliable"] = eb.estimate_unreliable;
        s["error_budget"] = budget;

        py::list phred;
        for (int q = 0; q <= 60; ++q) {
            if (stats.sample_phred[i][static_cast<size_t>(q)]) {
                py::dict d;
                d["phred"] = q;
                d["bases"] = stats.sample_phred[i][static_cast<size_t>(q)];
                phred.append(d);
            }
        }
        s["umi_phred"] = phred;
        s["umis_merged"] = cr.merged;
        s["reads_merged"] = cr.merged_reads;
        s["molecules_observed"] = cr.molecules_observed;
        s["molecules_corrected"] = cr.molecules_corrected;
        s["saturated"] = cr.saturated;
        per_sample.append(s);
    }
    out["samples"] = per_sample;
    out["wall_seconds"] = clock.seconds();
    out["reads_per_second"] =
        clock.seconds() > 0.0 ? static_cast<double>(c.total) / clock.seconds() : 0.0;
    return out;
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "migec native core: FASTQ IO, the .mig intermediate format, barcode primitives.";
    m.attr("__version__") = MIGEC_VERSION;
    m.attr("MIG_FORMAT_VERSION") = kMigFormatVersion;
    // The two measured constants, exported so the CLI's help text and defaults read them rather
    // than repeating them. X2 fitted the floor, X3 the split threshold; both live in
    // consensus.hpp next to the comment that says what measured them.
    m.attr("RT_FLOOR") = ConsensusParams{}.rt_floor;
    m.attr("LINKAGE_THRESHOLD") = ConsensusParams{}.linkage_threshold;
    // Also measured, and also in one place: level 6 spent 83% of refine's wall clock compressing.
    m.attr("GZIP_LEVEL") = AssembleRequest{}.gzip_level;

    py::register_exception<MigecError>(m, "MigecError", PyExc_RuntimeError);

    py::class_<PyMigRecord>(m, "MigRecord")
        .def_readonly("cell", &PyMigRecord::cell)
        .def_readonly("umi", &PyMigRecord::umi)
        .def_readonly("src_index", &PyMigRecord::src_index)
        .def_readonly("flags", &PyMigRecord::flags)
        .def_readonly("umi_minq", &PyMigRecord::umi_minq)
        .def_readonly("cell_minq", &PyMigRecord::cell_minq)
        .def_readonly("seq1", &PyMigRecord::seq1)
        .def_readonly("qual1", &PyMigRecord::qual1)
        .def_readonly("seq2", &PyMigRecord::seq2)
        .def_readonly("qual2", &PyMigRecord::qual2)
        .def("__repr__", [](const PyMigRecord& r) {
            return "<MigRecord cell=" + r.cell + " umi=" + r.umi + " len=" +
                   std::to_string(r.seq1.size()) + "/" + std::to_string(r.seq2.size()) + ">";
        });

    py::class_<PyMigFile>(m, "MigFile",
                          "Reads a .mig file. read_all() materialises every record, so it is for "
                          "tests, notebooks and small slices -- not for a production pass.")
        .def(py::init<const std::string&>(), py::arg("path"))
        .def_property_readonly("header", &PyMigFile::header)
        .def("read_all", &PyMigFile::read_all);

    m.def(
        "pack_barcode",
        [](const std::string& seq) {
            bool has_n = false;
            uint64_t p = pack_barcode(seq, &has_n);
            return py::make_tuple(p, has_n);
        },
        py::arg("seq"),
        "Pack up to 32 bases into a uint64, 2 bits per base, base 0 in the high bits so that "
        "packed order equals lexicographic order. Returns (packed, has_n); an N is stored as A "
        "and reported through the flag.");

    m.def("unpack_barcode", &unpack_barcode, py::arg("packed"), py::arg("length"));

    m.def("bucket_of", &bucket_of, py::arg("key"), py::arg("bits"),
          "Range-partition bucket for a packed barcode: the top `bits` bits. Order preserving, "
          "unlike a hash, which is what lets a barcode and its 1-mismatch neighbours be corrected "
          "within one bucket.");

    m.def(
        "reverse_complement",
        [](std::string seq, std::string qual) {
            reverse_complement(seq, qual);
            return py::make_tuple(seq, qual);
        },
        py::arg("seq"), py::arg("qual") = std::string(),
        "Reverse-complement a sequence and reverse its quality string together. They always "
        "travel together because reversing one and not the other is the classic silent bug.");

    m.def(
        "count_fastq",
        [](const std::string& path) {
            FastqReader r(path);
            FastqRecord rec;
            uint64_t n = 0;
            {
                py::gil_scoped_release release;
                while (r.next(rec)) ++n;
            }
            return n;
        },
        py::arg("path"), "Count records in a (optionally gzipped) FASTQ. Validates as it goes.");

    m.def("run_checkout", &py_run_checkout, py::arg("in_path"), py::arg("in_path2") = std::string(),
          py::arg("sample_ids") = std::vector<std::string>(),
          py::arg("patterns") = std::vector<std::string>(),
          py::arg("slaves") = std::vector<std::string>(), py::arg("out_prefix") = std::string(),
          py::arg("trim") = "pattern", py::arg("min_umi_quality") = 0,
          py::arg("write_unmatched") = false, py::arg("threads") = 0,
          py::arg("max_offset") = -1, py::arg("limit_reads") = uint64_t{0},
          "Demultiplex a FASTQ (or an R1/R2 pair) by barcode pattern, extract and trim the UMI, "
          "and write one gzipped FASTQ per sample with RX/QX/BC tags in the header. Returns a "
          "summary dict including the per-sample coverage histogram, UMI composition, correction "
          "statistics, and the wall clock, throughput and peak RSS of the run. The whole file is "
          "processed in C++ with the GIL released, across `threads` workers (0 = one per core); "
          "the output is byte-identical whatever that is set to.");

    m.def(
        "correct_umis",
        [](const std::vector<std::string>& umis, const std::vector<std::string>& payloads,
           const std::vector<std::string>& umi_quals, double sequencing_error,
           double polymerase_error, int pcr_cycles, double min_posterior,
           double max_child_fraction) {
            if (umis.empty()) throw MigecError("correct_umis: no UMIs given");
            const int len = static_cast<int>(umis.front().size());
            if (!payloads.empty() && payloads.size() != umis.size()) {
                throw MigecError("correct_umis: payloads and umis must be one per read");
            }
            if (!umi_quals.empty() && umi_quals.size() != umis.size()) {
                throw MigecError("correct_umis: umi_quals and umis must be one per read");
            }
            UmiCounts counts(len);
            for (const std::string& u : umis) {
                if (static_cast<int>(u.size()) != len) {
                    throw MigecError("correct_umis: mixed UMI lengths (" +
                                     std::to_string(len) + " and " + std::to_string(u.size()) +
                                     ")");
                }
                counts.add(pack_barcode(u));
            }

            // Per-barcode evidence, accumulated over that barcode's reads. Built here rather than
            // asked of the caller: it is derived from the same reads and getting it out of step
            // with `counts` would silently mis-attribute one barcode's reads to another.
            const std::vector<UmiCounts::Entry>& entries = counts.entries();
            BarcodeEvidence evidence;
            std::vector<uint32_t> seen(entries.size(), 0);
            auto slot = [&entries](uint64_t key) {
                auto it = std::lower_bound(
                    entries.begin(), entries.end(), key,
                    [](const UmiCounts::Entry& e, uint64_t k) { return e.key < k; });
                return static_cast<size_t>(it - entries.begin());
            };
            if (!umi_quals.empty()) {
                evidence.position_error.assign(entries.size() * static_cast<size_t>(len), 0.0f);
                for (size_t r = 0; r < umis.size(); ++r) {
                    const size_t i = slot(pack_barcode(umis[r]));
                    for (int j = 0; j < len && j < static_cast<int>(umi_quals[r].size()); ++j) {
                        evidence.position_error[i * static_cast<size_t>(len) +
                                                static_cast<size_t>(j)] +=
                            static_cast<float>(phred_error(phred_from_char(umi_quals[r][j])));
                    }
                }
                for (size_t i = 0; i < entries.size(); ++i) {
                    for (int j = 0; j < len; ++j) {
                        evidence.position_error[i * static_cast<size_t>(len) +
                                                static_cast<size_t>(j)] /=
                            static_cast<float>(entries[i].count);
                    }
                }
            }
            if (!payloads.empty()) {
                size_t width = payloads.front().size();
                for (const std::string& p : payloads) width = std::min(width, p.size());
                evidence.payload_width = static_cast<int>(width);
                // A modal base per column would need a counter per barcode per column; the first
                // read of each barcode is a draft, and a draft is all the LR needs -- it is
                // comparing two sequences, not calling variants.
                evidence.payload.assign(entries.size() * width, 0);
                for (size_t r = 0; r < umis.size(); ++r) {
                    const size_t i = slot(pack_barcode(umis[r]));
                    if (seen[i]++) continue;
                    std::copy_n(payloads[r].begin(), width, evidence.payload.begin() +
                                static_cast<std::ptrdiff_t>(i * width));
                }
            }
            CorrectionParams params;
            params.sequencing_error = sequencing_error;
            params.polymerase_error = polymerase_error;
            params.pcr_cycles = pcr_cycles;
            params.min_posterior = min_posterior;
            params.max_child_fraction = max_child_fraction;
            CorrectionResult r;
            {
                py::gil_scoped_release release;
                r = correct_umis(counts, params, evidence);
            }
            py::dict d;
            d["estimated_error"] = r.estimated_error;
            d["merged"] = r.merged;
            d["merged_reads"] = r.merged_reads;
            d["molecules_observed"] = r.molecules_observed;
            d["molecules_corrected"] = r.molecules_corrected;
            d["saturated"] = r.saturated;
            d["payload_clonality"] = r.payload_clonality;
            d["merged_by_payload"] = r.merged_by_payload;
            // The merge map, as barcodes. This is what a caller needs to apply the correction to
            // reads, and what any evaluation against a known truth has to have.
            py::list merges;
            py::dict corrected;
            for (size_t i = 0; i < entries.size(); ++i) {
                const std::string umi = unpack_barcode(entries[i].key, len);
                if (r.root[i] != i) {
                    py::dict e;
                    e["from"] = umi;
                    e["to"] = unpack_barcode(entries[r.root[i]].key, len);
                    e["reads"] = entries[i].count;
                    merges.append(e);
                } else {
                    corrected[py::str(umi)] = r.corrected[i];
                }
            }
            d["merges"] = merges;
            d["corrected"] = corrected;
            return d;
        },
        py::arg("umis"), py::arg("payloads") = std::vector<std::string>(),
        py::arg("umi_quals") = std::vector<std::string>(), py::arg("sequencing_error") = -1.0,
        py::arg("polymerase_error") = 1e-5, py::arg("pcr_cycles") = 25,
        py::arg("min_posterior") = 0.95, py::arg("max_child_fraction") = 0.5,
        "Fold error-child barcodes into their parents. `umis` is one entry per read. Returns the "
        "merge map, the corrected per-barcode counts, the error rate used, and the "
        "collision-corrected molecule count. A barcode with no plausible parent is always kept.");

    m.def(
        "subsample",
        [](const std::string& input, const std::string& output, uint32_t per_10k, bool by_cell,
           int gzip_level) {
            SubsampleRequest req;
            req.input = input;
            req.output = output;
            req.per_10k = per_10k;
            req.by_cell = by_cell;
            req.gzip_level = gzip_level;
            SubsampleStats st;
            {
                py::gil_scoped_release release;
                st = subsample(req);
            }
            py::dict d;
            d["reads"] = st.reads;
            d["reads_kept"] = st.reads_kept;
            d["barcodes"] = st.barcodes_seen;
            d["reads_without_umi"] = st.reads_without_umi;
            d["reads_per_barcode_median"] = st.reads_per_barcode_median;
            d["reads_per_barcode_max"] = st.reads_per_barcode_max;
            py::list examples;
            for (const auto& [barcode, depth] : st.examples) {
                examples.append(py::make_tuple(barcode, depth));
            }
            d["examples"] = examples;
            d["wall_seconds"] = st.wall_seconds;
            return d;
        },
        py::arg("input"), py::arg("output"), py::arg("per_10k") = 100u,
        py::arg("by_cell") = true, py::arg("gzip_level") = AssembleRequest{}.gzip_level,
        "Keep ALL the reads of a hash-selected fraction of the barcodes. Never a fraction of the "
        "reads, and never the first N barcodes: both destroy the MIG size distribution the "
        "fixture exists to preserve.");

    m.def(
        "refine",
        [](const std::string& input, const std::string& out_dir, const std::string& sample_id,
           bool use_quality, bool use_payload, int payload_width, double min_posterior,
           int expect_cells, const std::string& cell_whitelist, double min_whitelist_posterior,
           double target_fdr, int gzip_level, int threads, uint64_t limit_reads,
           uint64_t limit_umis) {
            RefineRequest req;
            req.input = input;
            req.output_dir = out_dir;
            req.sample_id = sample_id;
            req.use_quality = use_quality;
            req.use_payload = use_payload;
            req.payload_width = payload_width;
            req.correction.min_posterior = min_posterior;
            req.expect_cells = expect_cells;
            req.cell_whitelist = cell_whitelist;
            req.whitelist.min_posterior = min_whitelist_posterior;
            req.target_fdr = target_fdr;
            req.gzip_level = gzip_level;
            req.correction.threads = threads;
            req.limit.reads = limit_reads;
            req.limit.umis = limit_umis;
            RefineStats st;
            {
                py::gil_scoped_release release;
                st = refine(req);
            }
            py::dict d;
            d["sample_id"] = st.sample_id;
            d["reads"] = st.reads;
            d["reads_without_umi"] = st.reads_without_umi;
            d["barcodes"] = st.barcodes;
            d["merged"] = st.merged;
            d["merged_reads"] = st.merged_reads;
            d["merged_by_payload"] = st.merged_by_payload;
            d["molecules"] = st.molecules;
            d["molecules_corrected"] = st.molecules_corrected;
            d["estimated_error"] = st.estimated_error;
            d["payload_clonality"] = st.payload_clonality;
            d["saturated"] = st.saturated;
            d["umi_length"] = st.umi_length;
            d["cell_length"] = st.cell_length;
            d["cells_observed"] = st.cells_observed;
            d["cells_called"] = st.cells_called;
            d["molecules_in_called"] = st.molecules_in_called;
            d["cell_threshold"] = st.cell_threshold;
            d["knee_rank"] = st.knee_rank;
            d["knee_molecules"] = st.knee_molecules;
            d["mig_size_threshold"] = st.mig_size_threshold;
            d["residual_fdr_at_one"] = st.residual_fdr_at_one;
            d["suspected_residual"] = st.suspected_residual;
            d["error_from_children"] = st.error_from_children;
            d["error_at_depth"] = st.error_at_depth;
            d["error_depth"] = st.error_depth;
            d["error_phred"] = st.error_phred;
            py::dict wl;
            wl["barcodes"] = st.whitelist.barcodes;
            wl["exact"] = st.whitelist.exact;
            wl["corrected"] = st.whitelist.corrected;
            wl["off_list"] = st.whitelist.off_list;
            wl["reads_corrected"] = st.whitelist.reads_corrected;
            wl["far"] = st.whitelist.far;
            wl["background_prior"] = st.whitelist.background_prior;
            d["whitelist"] = wl;
            d["table_bytes"] = st.table_bytes;
            d["wall_seconds"] = st.wall_seconds;
            d["limited"] = st.limited;
            d["table_seconds"] = st.table_seconds;
            d["correct_seconds"] = st.correct_seconds;
            d["rewrite_seconds"] = st.rewrite_seconds;
            d["threads"] = st.threads;
            d["peak_rss_bytes"] = peak_rss_bytes();
            py::list hist;
            for (size_t b = 0; b < st.size_histogram.size(); ++b) {
                py::dict e;
                e["min_reads"] = static_cast<uint64_t>(1) << b;
                e["max_reads"] = (static_cast<uint64_t>(1) << (b + 1)) - 1;
                e["molecules"] = st.size_histogram[b];
                hist.append(e);
            }
            d["coverage"] = hist;
            return d;
        },
        py::arg("input"), py::arg("out_dir"), py::arg("sample_id") = std::string(),
        py::arg("use_quality") = true, py::arg("use_payload") = true,
        py::arg("payload_width") = 32, py::arg("min_posterior") = 0.95,
        py::arg("expect_cells") = 3000, py::arg("cell_whitelist") = std::string(),
        py::arg("min_whitelist_posterior") = 0.975, py::arg("target_fdr") = 0.05,
        py::arg("gzip_level") = RefineRequest{}.gzip_level, py::arg("threads") = 0,
        py::arg("limit_reads") = uint64_t{0}, py::arg("limit_umis") = uint64_t{0},
        "Correct barcode errors and rewrite the reads with the corrected barcode. Holds the "
        "barcode table, never the reads. A merged read keeps what it was in an OX:Z: tag, so the "
        "correction can be audited.");

    m.def(
        "assemble",
        [](const std::string& input, const std::string& out_dir, const std::string& sample_id,
           double rt_floor, double linkage_threshold, bool contig, bool fast, uint32_t min_reads,
           int gzip_level, int bucket_bits, int threads, uint64_t limit_reads,
           uint64_t limit_umis) {
            AssembleRequest req;
            req.input = input;
            req.output_dir = out_dir;
            req.sample_id = sample_id;
            req.consensus.rt_floor = rt_floor;
            req.consensus.linkage_threshold = linkage_threshold;
            req.consensus.contig = contig;
            req.consensus.fast = fast;
            req.min_reads = min_reads;
            req.gzip_level = gzip_level;
            req.bucket_bits = bucket_bits;
            req.threads = threads;
            req.limit.reads = limit_reads;
            req.limit.umis = limit_umis;
            AssembleStats st;
            {
                py::gil_scoped_release release;
                st = assemble(req);
            }
            py::dict d;
            d["sample_id"] = st.sample_id;
            d["reads"] = st.reads;
            d["reads_without_umi"] = st.reads_without_umi;
            d["reads_dropped"] = st.reads_dropped;
            d["groups"] = st.groups;
            d["molecules"] = st.molecules;
            d["groups_split"] = st.groups_split;
            d["groups_fragmented"] = st.groups_fragmented;
            d["contigs"] = st.contigs;
            d["contig_mode"] = contig;
            d["fast_mode"] = fast;
            d["mean_support"] = st.mean_support;
            d["groups_capped"] = st.groups_capped;
            d["reads_over_cap"] = st.reads_over_cap;
            d["max_reads_per_group"] = static_cast<uint64_t>(kMaxReadsPerGroup);
            d["cell_length"] = st.cell_length;
            d["expected_molecules_per_group"] = st.expected_molecules_per_group;
            py::dict sp;
            sp["effective_space"] = st.space.effective_space;
            sp["effective_length"] = st.space.effective_length;
            sp["occupancy"] = st.space.occupancy;
            sp["lambda"] = st.space.lambda;
            sp["molecules"] = st.space.molecules;
            sp["p_multi"] = st.space.p_multi;
            sp["saturated"] = st.space.saturated;
            d["barcode_space"] = sp;
            d["umi_length"] = st.umi_length;
            d["buckets"] = st.buckets;
            d["threads"] = st.threads;
            d["limited"] = st.limited;
            d["mean_quality"] = st.mean_quality;
            d["mean_consensus_error"] = st.mean_consensus_error;
            d["quality_cap"] = -10.0 * std::log10(rt_floor);
            d["wall_seconds"] = st.wall_seconds;
            d["partition_seconds"] = st.partition_seconds;
            d["peak_rss_bytes"] = peak_rss_bytes();
            py::list hist;
            for (size_t b = 0; b < st.size_histogram.size(); ++b) {
                py::dict e;
                e["min_reads"] = static_cast<uint64_t>(1) << b;
                e["max_reads"] = (static_cast<uint64_t>(1) << (b + 1)) - 1;
                e["groups"] = st.size_histogram[b];
                hist.append(e);
            }
            d["quality_grid"] = [&] {
                py::list rows;
                for (size_t b = 0; b < st.quality_grid.size(); ++b) {
                    for (size_t q = 0; q < kQualityLevels; ++q) {
                        if (!st.quality_grid[b][q]) continue;
                        py::dict e;
                        e["min_reads"] = static_cast<uint64_t>(1) << b;
                        e["max_reads"] = (static_cast<uint64_t>(1) << (b + 1)) - 1;
                        e["quality"] = static_cast<uint64_t>(q);
                        e["molecules"] = st.quality_grid[b][q];
                        rows.append(e);
                    }
                }
                return rows;
            }();
            d["coverage"] = hist;
            return d;
        },
        py::arg("input"), py::arg("out_dir"), py::arg("sample_id") = std::string(),
        // Defaulted off the header, never retyped: rt_floor is X2's measurement and
        // linkage_threshold is X3's, and a second copy of a measured constant is a copy that
        // drifts from the measurement that justifies it.
        py::arg("rt_floor") = ConsensusParams{}.rt_floor,
        py::arg("linkage_threshold") = ConsensusParams{}.linkage_threshold,
        py::arg("contig") = false, py::arg("fast") = false,
        py::arg("min_reads") = 1u, py::arg("gzip_level") = AssembleRequest{}.gzip_level,
        py::arg("bucket_bits") = 0, py::arg("threads") = 0,
        py::arg("limit_reads") = uint64_t{0}, py::arg("limit_umis") = uint64_t{0},
        "Collapse the reads of each UMI into a consensus. Reads are range partitioned into .mig "
        "buckets and sorted one bucket at a time, so nothing scales with the library. Emitted "
        "quality is capped at -10 log10(rt_floor), the RT/first-cycle-PCR error that no consensus "
        "removes; a group is split into two molecules only when the co-segregation of its minor "
        "alleles exceeds `linkage_threshold`, which is a measured false-positive point.");

    m.def(
        "suggest",
        [](const std::string& path, int cycles, uint64_t max_reads, double umi_deviation) {
            Suggestion sg;
            {
                py::gil_scoped_release release;
                sg = suggest_pattern(profile_cycles(path, cycles, max_reads), umi_deviation);
            }
            py::dict d;
            d["pattern"] = sg.pattern;
            d["umi_length"] = sg.umi_length;
            d["anchor_length"] = sg.anchor_length;
            d["note"] = sg.note;
            d["reads"] = sg.profile.reads;
            d["read_length"] = sg.profile.read_length;
            py::list cyc;
            for (const CycleStats& c : sg.profile.cycles) {
                const std::array<double, 4> f = c.frequencies();
                py::dict e;
                e["cycle"] = c.cycle;
                e["A"] = f[0];
                e["C"] = f[1];
                e["G"] = f[2];
                e["T"] = f[3];
                e["entropy"] = c.entropy();
                e["collision"] = c.collision();
                e["consensus"] = std::string(1, c.consensus());
                e["consensus_fraction"] = c.consensus_fraction();
                e["deviation"] = c.deviation_from_uniform() / 2.0;
                e["mean_phred"] = c.mean_phred();
                cyc.append(e);
            }
            d["cycles"] = cyc;
            py::list km;
            for (const KmerHit& h : sg.profile.kmers) {
                py::dict e;
                e["kmer"] = h.kmer;
                e["count"] = h.count;
                e["expected"] = h.expected;
                e["ratio"] = h.ratio;
                e["mean_position"] = h.mean_position;
                e["read_fraction"] = h.read_fraction;
                km.append(e);
            }
            d["kmers"] = km;
            d["kmers_scanned"] = sg.profile.kmers_scanned;
            py::list segs;
            for (const Segment& sm : sg.segments) {
                py::dict e;
                e["kind"] = sm.kind == CycleKind::kUmi ? "umi"
                            : sm.kind == CycleKind::kConstant ? "constant" : "variable";
                e["begin"] = sm.begin;
                e["end"] = sm.end;
                e["length"] = sm.length();
                e["consensus"] = sm.consensus;
                e["mean_deviation"] = sm.mean_deviation;
                segs.append(e);
            }
            d["segments"] = segs;
            return d;
        },
        py::arg("path"), py::arg("cycles") = 60, py::arg("max_reads") = 200000,
        py::arg("umi_deviation") = 0.18,
        "Read the barcode layout off the reads: per-cycle base composition, segmented into UMI "
        "(all four bases near 1/4), constant (one base near 1) and variable runs, and a "
        "paste-ready pattern.");

    m.def("peak_rss_bytes", &peak_rss_bytes,
          "Peak resident set size of this process in bytes, 0 if the platform will not say.");
    m.def("hardware_threads", &hardware_threads);

    m.def(
        "match_pattern",
        [](const std::string& spec, const std::string& seq, const std::string& qual,
           int max_offset) {
            BarcodePattern p = BarcodePattern::compile(spec);
            MatchParams mp;
            mp.max_offset = max_offset;
            PatternMatch m2 = p.match(seq, qual, mp);
            py::dict d;
            d["found"] = m2.found;
            d["offset"] = m2.offset;
            d["score"] = m2.score;
            d["margin"] = m2.margin;
            d["umi"] = m2.umi;
            d["umi_qual"] = m2.umi_qual;
            // The cell barcode is captured separately by `X` positions, and was missing here
            // until the audit: a pattern could extract one and no Python caller could see it.
            d["cell"] = m2.cell;
            d["cell_qual"] = m2.cell_qual;
            d["payload_begin"] = m2.payload_begin;
            return d;
        },
        py::arg("spec"), py::arg("seq"), py::arg("qual") = std::string(),
        py::arg("max_offset") = -1,
        "Match one barcode pattern against one read. For inspecting a pattern interactively; the "
        "per-read path in a real run stays in C++.");

    m.def(
        "umi_statistics",
        [](const std::vector<std::string>& umis) {
            if (umis.empty()) throw MigecError("umi_statistics: no UMIs given");
            UmiCounts counts(static_cast<int>(umis[0].size()));
            for (const auto& u : umis) {
                if (u.size() != umis[0].size()) {
                    throw MigecError("umi_statistics: UMIs of differing length (" +
                                     std::to_string(umis[0].size()) + " vs " +
                                     std::to_string(u.size()) + ")");
                }
                counts.add(pack_barcode(u));
            }
            const CoverageHistogram h = counts.histogram();
            const UmiComposition comp = counts.composition(false);
            py::dict d;
            d["distinct"] = counts.distinct();
            d["total"] = counts.total();
            d["mean_reads_per_umi"] = h.mean_reads_per_umi();
            d["hist_reads"] = h.reads;
            d["hist_units"] = h.units;
            d["total_entropy"] = comp.total_entropy();
            d["total_information"] = comp.total_information();
            d["effective_length"] = comp.effective_length();
            d["effective_space"] = comp.effective_space();
            py::list per_pos;
            for (int j = 0; j < comp.length; ++j) {
                py::dict e;
                e["position"] = j;
                e["A"] = comp.freq[static_cast<size_t>(j)][0];
                e["C"] = comp.freq[static_cast<size_t>(j)][1];
                e["G"] = comp.freq[static_cast<size_t>(j)][2];
                e["T"] = comp.freq[static_cast<size_t>(j)][3];
                e["entropy"] = comp.entropy(j);
                e["information"] = comp.information(j);
                e["collision"] = comp.collision(j);
                per_pos.append(e);
            }
            d["composition"] = per_pos;
            return d;
        },
        py::arg("umis"),
        "Coverage histogram, per-position composition, entropy and information for a list of "
        "observed UMI strings (one entry per read).");
}
