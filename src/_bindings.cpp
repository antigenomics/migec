// Python bindings. Deliberately thin: everything here is either a whole-file operation (so the
// per-call overhead is irrelevant) or a small utility used by tests and notebooks. Nothing in a
// per-read loop crosses this boundary.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <string>
#include <vector>

#include "migec/fastq.hpp"
#include "migec/mig_record.hpp"
#include "migec/types.hpp"
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
}
