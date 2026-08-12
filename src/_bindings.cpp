// Python bindings. Deliberately thin: everything here is either a whole-file operation (so the
// per-call overhead is irrelevant) or a small utility used by tests and notebooks. Nothing in a
// per-read loop crosses this boundary.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "migec/checkout.hpp"
#include "migec/fastq.hpp"
#include "migec/mig_record.hpp"
#include "migec/pattern.hpp"
#include "migec/resource.hpp"
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
                         const std::vector<std::string>& patterns, const std::string& out_prefix,
                         const std::string& trim, int min_umi_quality, bool write_unmatched,
                         int threads) {
    if (sample_ids.size() != patterns.size()) {
        throw MigecError("run_checkout: sample_ids and patterns have different lengths");
    }
    PatternSet set;
    for (size_t i = 0; i < patterns.size(); ++i) set.add(sample_ids[i], patterns[i]);

    CheckoutParams params;
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
    out["match_seconds"] = stats.wall_seconds;
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

        const CorrectionResult cr = correct_umis(umi_counts[i]);
        s["umi_error_rate"] = cr.estimated_error;
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
          py::arg("patterns") = std::vector<std::string>(), py::arg("out_prefix") = std::string(),
          py::arg("trim") = "pattern", py::arg("min_umi_quality") = 0,
          py::arg("write_unmatched") = false, py::arg("threads") = 0,
          "Demultiplex a FASTQ (or an R1/R2 pair) by barcode pattern, extract and trim the UMI, "
          "and write one gzipped FASTQ per sample with RX/QX/BC tags in the header. Returns a "
          "summary dict including the per-sample coverage histogram, UMI composition, correction "
          "statistics, and the wall clock, throughput and peak RSS of the run. The whole file is "
          "processed in C++ with the GIL released, across `threads` workers (0 = one per core); "
          "the output is byte-identical whatever that is set to.");

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
