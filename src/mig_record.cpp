// .mig reader/writer. See include/migec/mig_record.hpp for why the layout is what it is.
//
// ponytail: blocks are compressed with zlib deflate level 1, not zstd-1. We already link zlib for
// gzipped FASTQ, and libzstd would mean a second system dependency to satisfy in manylinux and
// macOS wheels for a measured difference of a few percent on this payload. The block header
// carries a codec byte, so switching to zstd later is a new codec id, not a format break.

#include "migec/mig_record.hpp"

#include <zlib.h>

#include <cstdio>
#include <cstring>
#include <vector>

#include "migec/types.hpp"

namespace migec {
namespace {

constexpr uint8_t kCodecNone = 0;
constexpr uint8_t kCodecDeflate1 = 1;

// Fixed part of a record, as stored. Sequence and quality live in the column-major payload.
#pragma pack(push, 1)
struct StoredRecord {
    uint64_t cell;
    uint64_t umi;
    uint64_t src_index;
    uint16_t flags;
    uint8_t umi_minq;
    uint8_t cell_minq;
    uint32_t len1;
    uint32_t len2;
};
struct BlockHeader {
    uint32_t n_records;
    uint32_t raw_bytes;
    uint32_t stored_bytes;
    uint32_t crc32;  // of the uncompressed payload
    uint8_t codec;
    uint8_t reserved[3];
};
#pragma pack(pop)

void put_u16(std::string& s, uint16_t v) { s.append(reinterpret_cast<const char*>(&v), 2); }
void put_u32(std::string& s, uint32_t v) { s.append(reinterpret_cast<const char*>(&v), 4); }
void put_u64(std::string& s, uint64_t v) { s.append(reinterpret_cast<const char*>(&v), 8); }
void put_str(std::string& s, const std::string& v) {
    put_u32(s, static_cast<uint32_t>(v.size()));
    s += v;
}

std::string serialize_header(const MigHeader& h) {
    std::string s(kMigMagic, 4);
    put_u16(s, h.format_version);
    s.push_back(static_cast<char>(h.umi_len));
    s.push_back(static_cast<char>(h.cell_len));
    s.push_back(static_cast<char>(h.bucket_index));
    s.push_back(static_cast<char>(h.bucket_bits));
    s.push_back(static_cast<char>(h.paired ? 1 : 0));
    s.push_back(static_cast<char>(h.barcode_quality ? 1 : 0));  // v2
    put_str(s, h.sample_id);
    put_str(s, h.provenance);
    put_u32(s, static_cast<uint32_t>(h.quality_calibration.size()));
    for (float f : h.quality_calibration) s.append(reinterpret_cast<const char*>(&f), 4);
    return s;
}

class Reader {
public:
    explicit Reader(const std::string& path) : path_(path) {
        f_ = std::fopen(path.c_str(), "rb");
        if (!f_) throw MigecError("mig_reader: cannot open " + path);
    }
    ~Reader() { if (f_) std::fclose(f_); }

    void must_read(void* dst, size_t n, const char* what) {
        if (std::fread(dst, 1, n, f_) != n) {
            throw MigecError(std::string("mig_reader: truncated file, wanted ") + what + " at offset " +
                             std::to_string(offset()) + " in " + path_);
        }
    }
    long offset() { return std::ftell(f_); }
    std::FILE* f() { return f_; }

private:
    std::string path_;
    std::FILE* f_ = nullptr;
};

}  // namespace

// ------------------------------------------------------------------------------------------- //

struct MigWriter::Impl {
    std::FILE* f = nullptr;
    std::string path;   // for the error message; a path in a message is what makes it actionable
    MigHeader header;
    size_t block_bytes;
    uint64_t n_written = 0;

    std::vector<StoredRecord> pending;
    std::string seq1, seq2, qual1, qual2;
    // The barcode quality columns, fixed width per record. Empty unless the header says the file
    // carries them.
    std::string qumi, qcell;
    std::string zbuf;

    void flush() {
        if (pending.empty()) return;
        // Column-major: the fixed records, then all seq1, all seq2, all qual1, all qual2.
        std::string raw;
        raw.reserve(pending.size() * sizeof(StoredRecord) + seq1.size() + seq2.size() +
                    qual1.size() + qual2.size() + qumi.size() + qcell.size());
        raw.append(reinterpret_cast<const char*>(pending.data()),
                   pending.size() * sizeof(StoredRecord));
        raw += seq1;
        raw += seq2;
        raw += qual1;
        raw += qual2;
        raw += qumi;   // empty unless header.barcode_quality
        raw += qcell;

        uLongf bound = compressBound(static_cast<uLong>(raw.size()));
        zbuf.resize(bound);
        uLongf out_len = bound;
        int rc = compress2(reinterpret_cast<Bytef*>(zbuf.data()), &out_len,
                           reinterpret_cast<const Bytef*>(raw.data()),
                           static_cast<uLong>(raw.size()), 1);
        if (rc != Z_OK) throw MigecError("mig_writer: zlib compress failed");

        BlockHeader bh{};
        bh.n_records = static_cast<uint32_t>(pending.size());
        bh.raw_bytes = static_cast<uint32_t>(raw.size());
        bh.stored_bytes = static_cast<uint32_t>(out_len);
        bh.crc32 = static_cast<uint32_t>(
            crc32(0L, reinterpret_cast<const Bytef*>(raw.data()), static_cast<uInt>(raw.size())));
        bh.codec = kCodecDeflate1;
        if (std::fwrite(&bh, sizeof(bh), 1, f) != 1 ||
            std::fwrite(zbuf.data(), 1, out_len, f) != out_len) {
            throw MigecError("mig_writer: short write");
        }

        pending.clear();
        seq1.clear(); seq2.clear(); qual1.clear(); qual2.clear();
        qumi.clear(); qcell.clear();
    }
};

MigWriter::MigWriter(const std::string& path, const MigHeader& header, size_t block_bytes)
    : impl_(new Impl) {
    impl_->f = std::fopen(path.c_str(), "wb");
    if (!impl_->f) throw MigecError("mig_writer: cannot open " + path + " for writing");
    impl_->path = path;
    impl_->header = header;
    impl_->header.format_version = kMigFormatVersion;
    impl_->block_bytes = block_bytes;
    const std::string hdr = serialize_header(impl_->header);
    if (std::fwrite(hdr.data(), 1, hdr.size(), impl_->f) != hdr.size()) {
        throw MigecError("mig_writer: cannot write header to " + path);
    }
}

MigWriter::~MigWriter() {
    if (impl_ && impl_->f) {
        try { close(); } catch (...) {}  // a throwing destructor is worse than a lost error
    }
}
MigWriter::MigWriter(MigWriter&&) noexcept = default;
MigWriter& MigWriter::operator=(MigWriter&&) noexcept = default;

void MigWriter::write(const MigRecord& rec) {
    Impl& im = *impl_;
    if (rec.qual1.size() != rec.seq1.size() || rec.qual2.size() != rec.seq2.size()) {
        throw MigecError("mig_writer: sequence and quality lengths differ");
    }
    StoredRecord sr{};
    sr.cell = rec.cell;
    sr.umi = rec.umi;
    sr.src_index = rec.src_index;
    sr.flags = rec.flags;
    sr.umi_minq = rec.umi_minq;
    sr.cell_minq = rec.cell_minq;
    sr.len1 = static_cast<uint32_t>(rec.seq1.size());
    sr.len2 = static_cast<uint32_t>(rec.seq2.size());
    im.pending.push_back(sr);
    im.seq1.append(rec.seq1);
    im.seq2.append(rec.seq2);
    im.qual1.append(rec.qual1);
    im.qual2.append(rec.qual2);
    if (im.header.barcode_quality) {
        // Fixed width, so a record that disagrees would shift every column after it and be read
        // back as another record's quality. Refused here, where the record is attributable.
        if (rec.qual_umi.size() != im.header.umi_len ||
            rec.qual_cell.size() != im.header.cell_len) {
            throw MigecError("mig_writer: this file carries barcode quality, so every record needs "
                             "a " + std::to_string(im.header.umi_len) + " nt UMI quality and a " +
                             std::to_string(im.header.cell_len) + " nt cell quality (got " +
                             std::to_string(rec.qual_umi.size()) + " and " +
                             std::to_string(rec.qual_cell.size()) + ")");
        }
        im.qumi.append(rec.qual_umi);
        im.qcell.append(rec.qual_cell);
    }
    ++im.n_written;

    if (im.pending.size() * sizeof(StoredRecord) + im.seq1.size() + im.seq2.size() +
            im.qual1.size() + im.qual2.size() + im.qumi.size() + im.qcell.size() >=
        im.block_bytes) {
        im.flush();
    }
}

void MigWriter::close() {
    Impl& im = *impl_;
    if (!im.f) return;
    im.flush();
    // Footer: a zero-record block header acts as the terminator, then the record count and the
    // magic again, so a reader can validate the file was closed cleanly.
    BlockHeader term{};
    term.codec = kCodecNone;
    // Never: the footer is what the reader checks to know the file was closed cleanly, and fclose
    // is where a full disk shows up -- the block writes above went into the C library's buffer.
    // Ignoring either reports a successful run over a file the reader will refuse.
    const bool wrote_footer = std::fwrite(&term, sizeof(term), 1, im.f) == 1 &&
                              std::fwrite(&im.n_written, sizeof(im.n_written), 1, im.f) == 1 &&
                              std::fwrite(kMigMagic, 1, 4, im.f) == 4;
    const bool closed = std::fclose(im.f) == 0;
    im.f = nullptr;
    if (!wrote_footer || !closed) {
        throw MigecError("mig_writer: could not finish " + im.path + " -- is the disk full?");
    }
}

uint64_t MigWriter::records_written() const { return impl_->n_written; }

// ------------------------------------------------------------------------------------------- //

struct MigReader::Impl {
    std::unique_ptr<Reader> r;
    MigHeader header;
    uint64_t declared = 0;

    std::string raw;                 // uncompressed payload arena for the current block
    std::vector<StoredRecord> recs;
    size_t pos = 0;                  // index into recs
    size_t off_seq1 = 0, off_seq2 = 0, off_qual1 = 0, off_qual2 = 0;
    size_t off_qumi = 0, off_qcell = 0;
    bool eof = false;

    bool load_block() {
        BlockHeader bh{};
        size_t got = std::fread(&bh, 1, sizeof(bh), r->f());
        if (got == 0) {
            // No footer at all: the writer never closed. That is a truncation, not an EOF.
            throw MigecError("mig_reader: file ends without a footer (writer did not close)");
        }
        if (got != sizeof(bh)) throw MigecError("mig_reader: truncated block header");
        if (bh.n_records == 0) {  // terminator
            uint64_t n = 0;
            char magic[4] = {0, 0, 0, 0};
            if (std::fread(&n, 1, 8, r->f()) == 8 && std::fread(magic, 1, 4, r->f()) == 4 &&
                std::memcmp(magic, kMigMagic, 4) == 0) {
                declared = n;
            } else {
                throw MigecError("mig_reader: corrupt footer");
            }
            eof = true;
            return false;
        }

        std::string stored(bh.stored_bytes, '\0');
        if (std::fread(stored.data(), 1, bh.stored_bytes, r->f()) != bh.stored_bytes) {
            throw MigecError("mig_reader: truncated block payload (declared " +
                             std::to_string(bh.stored_bytes) + " bytes)");
        }
        raw.assign(bh.raw_bytes, '\0');
        if (bh.codec == kCodecDeflate1) {
            uLongf out_len = bh.raw_bytes;
            int rc = uncompress(reinterpret_cast<Bytef*>(raw.data()), &out_len,
                                reinterpret_cast<const Bytef*>(stored.data()),
                                static_cast<uLong>(bh.stored_bytes));
            if (rc != Z_OK || out_len != bh.raw_bytes) {
                throw MigecError("mig_reader: block decompression failed");
            }
        } else if (bh.codec == kCodecNone) {
            raw = stored;
        } else {
            throw MigecError("mig_reader: unknown block codec " + std::to_string(bh.codec));
        }
        uint32_t crc = static_cast<uint32_t>(
            crc32(0L, reinterpret_cast<const Bytef*>(raw.data()), static_cast<uInt>(raw.size())));
        if (crc != bh.crc32) throw MigecError("mig_reader: block CRC mismatch (corrupt file)");

        // Never: BEFORE the memcpy. `n_records` is a u32 out of the file and the CRC above only
        // proves the payload is what the writer wrote -- a crafted header declaring more records
        // than the payload holds would read past the decompressed buffer.
        if (static_cast<size_t>(bh.n_records) * sizeof(StoredRecord) > raw.size()) {
            throw MigecError("mig_reader: block declares " + std::to_string(bh.n_records) +
                             " records, which do not fit its " + std::to_string(raw.size()) +
                             " byte payload");
        }
        recs.resize(bh.n_records);
        std::memcpy(recs.data(), raw.data(), bh.n_records * sizeof(StoredRecord));
        size_t o = bh.n_records * sizeof(StoredRecord);
        size_t l1 = 0, l2 = 0;
        for (const auto& s : recs) { l1 += s.len1; l2 += s.len2; }
        off_seq1 = o;
        off_seq2 = off_seq1 + l1;
        off_qual1 = off_seq2 + l2;
        off_qual2 = off_qual1 + l1;
        off_qumi = off_qual2 + l2;
        const size_t bq = header.barcode_quality ? bh.n_records : 0;
        off_qcell = off_qumi + bq * header.umi_len;
        if (off_qcell + bq * header.cell_len != raw.size()) {
            throw MigecError("mig_reader: block size inconsistent");
        }
        pos = 0;
        return true;
    }
};

MigReader::MigReader(const std::string& path) : impl_(new Impl) {
    impl_->r.reset(new Reader(path));
    char magic[4];
    impl_->r->must_read(magic, 4, "magic");
    if (std::memcmp(magic, kMigMagic, 4) != 0) throw MigecError("mig_reader: not a .mig file: " + path);
    auto& h = impl_->header;
    impl_->r->must_read(&h.format_version, 2, "format version");
    if (h.format_version < kMigMinReadableVersion || h.format_version > kMigFormatVersion) {
        throw MigecError("mig_reader: unsupported format version " +
                         std::to_string(h.format_version) + " (this build reads " +
                         std::to_string(kMigMinReadableVersion) + ".." +
                         std::to_string(kMigFormatVersion) + ")");
    }
    uint8_t b[5];
    impl_->r->must_read(b, 5, "header fields");
    h.umi_len = b[0]; h.cell_len = b[1]; h.bucket_index = b[2]; h.bucket_bits = b[3];
    h.paired = b[4] != 0;
    // v2 added one byte here. An older file simply does not carry the barcode quality columns,
    // and every reader of them checks the flag rather than the version.
    h.barcode_quality = false;
    if (h.format_version >= 2) {
        uint8_t bq = 0;
        impl_->r->must_read(&bq, 1, "barcode quality flag");
        h.barcode_quality = bq != 0;
    }
    auto read_str = [&](std::string& out) {
        uint32_t n = 0;
        impl_->r->must_read(&n, 4, "string length");
        out.assign(n, '\0');
        if (n) impl_->r->must_read(out.data(), n, "string body");
    };
    read_str(h.sample_id);
    read_str(h.provenance);
    // Never: `bucket_bits` and `bucket_index` are one byte each and go straight into a shift and a
    // vector index in every stage that reads a partition. `1 << 200` is undefined behaviour and
    // `bucket_paths[200]` on a two-bucket vector is an out-of-bounds write, both reachable from a
    // file that merely got corrupted in the right byte. The format's own limit is 8 bits (256
    // buckets); anything past it is not a migec partition.
    if (h.bucket_bits > 8) {
        throw MigecError("mig_reader: " + path + " declares 2^" + std::to_string(h.bucket_bits) +
                         " buckets, past the format's limit of 2^8");
    }
    if (h.bucket_index >= (1u << h.bucket_bits) && !(h.bucket_bits == 0 && h.bucket_index == 0)) {
        throw MigecError("mig_reader: " + path + " calls itself bucket " +
                         std::to_string(h.bucket_index) + " of 2^" +
                         std::to_string(h.bucket_bits) + ", which does not exist");
    }
    uint32_t nq = 0;
    impl_->r->must_read(&nq, 4, "calibration length");
    h.quality_calibration.resize(nq);
    if (nq) impl_->r->must_read(h.quality_calibration.data(), nq * 4, "calibration table");
}

MigReader::~MigReader() = default;
MigReader::MigReader(MigReader&&) noexcept = default;
MigReader& MigReader::operator=(MigReader&&) noexcept = default;

const MigHeader& MigReader::header() const { return impl_->header; }
uint64_t MigReader::records_declared() const { return impl_->declared; }

bool MigReader::next(MigRecord& out) {
    Impl& im = *impl_;
    while (im.pos >= im.recs.size()) {
        if (im.eof) return false;
        if (!im.load_block()) return false;
    }
    const StoredRecord& sr = im.recs[im.pos];
    out.cell = sr.cell;
    out.umi = sr.umi;
    out.src_index = sr.src_index;
    out.flags = sr.flags;
    out.umi_minq = sr.umi_minq;
    out.cell_minq = sr.cell_minq;
    out.seq1 = std::string_view(im.raw.data() + im.off_seq1, sr.len1);
    out.seq2 = std::string_view(im.raw.data() + im.off_seq2, sr.len2);
    out.qual1 = std::string_view(im.raw.data() + im.off_qual1, sr.len1);
    out.qual2 = std::string_view(im.raw.data() + im.off_qual2, sr.len2);
    if (im.header.barcode_quality) {
        out.qual_umi = std::string_view(im.raw.data() + im.off_qumi, im.header.umi_len);
        out.qual_cell = std::string_view(im.raw.data() + im.off_qcell, im.header.cell_len);
        im.off_qumi += im.header.umi_len;
        im.off_qcell += im.header.cell_len;
    } else {
        out.qual_umi = {};
        out.qual_cell = {};
    }
    im.off_seq1 += sr.len1;
    im.off_seq2 += sr.len2;
    im.off_qual1 += sr.len1;
    im.off_qual2 += sr.len2;
    ++im.pos;
    return true;
}

}  // namespace migec
