#include "migec/fastq.hpp"

#include <zlib.h>

#include <cstring>
#include <string>
#include <vector>

#include "migec/types.hpp"

namespace migec {
namespace {

constexpr size_t kChunk = 1u << 20;

// gzFile handles plain files too -- zlib falls back to a straight read when the magic is absent,
// so there is no second code path for uncompressed input.
struct GzIn {
    gzFile f = nullptr;
    explicit GzIn(const std::string& path) {
        f = gzopen(path.c_str(), "rb");
        if (!f) throw MigecError("fastq_reader: cannot open " + path);
        gzbuffer(f, 1u << 18);
    }
    ~GzIn() { if (f) gzclose(f); }
};

}  // namespace

struct FastqReader::Impl {
    std::string path;
    std::unique_ptr<GzIn> in;
    std::string buf;
    size_t pos = 0;     // parse cursor into buf
    size_t filled = 0;  // valid bytes in buf
    bool eof = false;
    uint64_t n = 0;

    // Guarantees at least one more byte in buf if the stream has one, compacting what is already
    // consumed. Returns false at true end of stream.
    bool fill() {
        if (eof) return false;
        if (pos > 0) {
            buf.erase(0, pos);
            filled -= pos;
            pos = 0;
        }
        size_t old = filled;
        buf.resize(old + kChunk);
        int got = gzread(in->f, buf.data() + old, static_cast<unsigned>(kChunk));
        if (got < 0) throw MigecError("fastq_reader: read error in " + path);
        filled = old + static_cast<size_t>(got);
        buf.resize(filled);
        if (got == 0) eof = true;
        return got > 0;
    }

    // Scans forward from `from` for one line, WITHOUT consuming it. Returns npos if the buffer
    // does not yet contain a full line.
    size_t scan_line(size_t from, std::string_view& out) const {
        const char* base = buf.data();
        const void* nl = std::memchr(base + from, '\n', filled - from);
        if (!nl) return std::string::npos;
        size_t idx = static_cast<size_t>(static_cast<const char*>(nl) - base);
        size_t len = idx - from;
        if (len && base[idx - 1] == '\r') --len;  // tolerate CRLF
        out = std::string_view(base + from, len);
        return idx + 1;
    }

    // Brings a whole 4-line record into the buffer and returns views into it. Every view is
    // produced only after the last possible fill(), because fill() compacts the buffer and would
    // otherwise leave the header view dangling -- which is the classic way this parser breaks on
    // exactly one record in ten million, at a buffer boundary.
    bool record(std::string_view lines[4]) {
        while (true) {
            size_t p = pos;
            bool complete = true;
            for (int i = 0; i < 4; ++i) {
                size_t nxt = scan_line(p, lines[i]);
                if (nxt == std::string::npos) { complete = false; break; }
                p = nxt;
            }
            if (complete) {
                pos = p;
                return true;
            }
            if (!fill()) {
                // Stream ended. Either cleanly between records, or mid-record.
                if (pos >= filled) return false;
                // Tolerate a missing final newline, but only for a record that is otherwise whole.
                size_t q = pos;
                for (int i = 0; i < 4; ++i) {
                    size_t nxt = scan_line(q, lines[i]);
                    if (nxt == std::string::npos) {
                        if (i == 3 && q < filled) {
                            lines[3] = std::string_view(buf.data() + q, filled - q);
                            pos = filled;
                            return true;
                        }
                        throw MigecError("fastq_reader: truncated record " + std::to_string(n + 1) +
                                         " in " + path);
                    }
                    q = nxt;
                }
                pos = q;
                return true;
            }
        }
    }
};

FastqReader::FastqReader(const std::string& path) : impl_(new Impl) {
    impl_->path = path;
    impl_->in.reset(new GzIn(path));
}
FastqReader::~FastqReader() = default;
FastqReader::FastqReader(FastqReader&&) noexcept = default;
FastqReader& FastqReader::operator=(FastqReader&&) noexcept = default;

uint64_t FastqReader::records_read() const { return impl_->n; }

bool FastqReader::next(FastqRecord& out) {
    Impl& im = *impl_;
    std::string_view lines[4];
    if (!im.record(lines)) return false;

    std::string_view head = lines[0];
    if (head.empty() || head[0] != '@') {
        throw MigecError("fastq_reader: expected '@' at record " + std::to_string(im.n + 1) +
                         " of " + im.path);
    }
    head.remove_prefix(1);
    size_t sp = head.find_first_of(" \t");
    if (sp == std::string_view::npos) {
        out.name = head;
        out.comment = std::string_view();
    } else {
        out.name = head.substr(0, sp);
        out.comment = head.substr(sp + 1);
    }

    const std::string_view seq = lines[1], plus = lines[2], qual = lines[3];
    if (plus.empty() || plus[0] != '+') {
        throw MigecError("fastq_reader: expected '+' in record " + std::to_string(im.n + 1) +
                         " of " + im.path);
    }
    if (seq.size() != qual.size()) {
        throw MigecError("fastq_reader: sequence and quality lengths differ in record " +
                         std::to_string(im.n + 1) + " of " + im.path);
    }
    out.seq = seq;
    out.qual = qual;
    ++im.n;
    return true;
}

// ------------------------------------------------------------------------------------------- //

struct FastqWriter::Impl {
    gzFile gz = nullptr;
    std::FILE* raw = nullptr;
    std::string buf;

    void put(std::string_view s) { buf.append(s); }
    void flush(bool force) {
        if (!force && buf.size() < (1u << 20)) return;
        if (buf.empty()) return;
        if (gz) {
            if (gzwrite(gz, buf.data(), static_cast<unsigned>(buf.size())) == 0) {
                throw MigecError("fastq_writer: gzip write failed");
            }
        } else if (std::fwrite(buf.data(), 1, buf.size(), raw) != buf.size()) {
            throw MigecError("fastq_writer: short write");
        }
        buf.clear();
    }
};

FastqWriter::FastqWriter(const std::string& path, int gzip_level) : impl_(new Impl) {
    const bool gzipped = path.size() > 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
    if (gzipped) {
        std::string mode = "wb" + std::to_string(gzip_level);
        impl_->gz = gzopen(path.c_str(), mode.c_str());
        if (!impl_->gz) throw MigecError("fastq_writer: cannot open " + path);
    } else {
        impl_->raw = std::fopen(path.c_str(), "wb");
        if (!impl_->raw) throw MigecError("fastq_writer: cannot open " + path);
    }
}
FastqWriter::~FastqWriter() {
    try { close(); } catch (...) {}
}
FastqWriter::FastqWriter(FastqWriter&&) noexcept = default;
FastqWriter& FastqWriter::operator=(FastqWriter&&) noexcept = default;

void FastqWriter::write(std::string_view name, std::string_view comment, std::string_view seq,
                        std::string_view qual) {
    if (seq.size() != qual.size()) {
        throw MigecError("fastq_writer: sequence and quality lengths differ for read " +
                         std::string(name));
    }
    Impl& im = *impl_;
    im.put("@");
    im.put(name);
    if (!comment.empty()) {
        im.put(" ");
        im.put(comment);
    }
    im.put("\n");
    im.put(seq);
    im.put("\n+\n");
    im.put(qual);
    im.put("\n");
    im.flush(false);
}

void FastqWriter::close() {
    Impl& im = *impl_;
    if (!im.gz && !im.raw) return;
    im.flush(true);
    if (im.gz) { gzclose(im.gz); im.gz = nullptr; }
    if (im.raw) { std::fclose(im.raw); im.raw = nullptr; }
}

}  // namespace migec
