// FASTQ IO. Plain or gzipped, decided by the first two bytes rather than the file name, because
// half the world writes gzipped data to a name ending in .fastq.

#ifndef MIGEC_FASTQ_HPP
#define MIGEC_FASTQ_HPP

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_set>
#include <string_view>

namespace migec {

// A record's fields are views into the reader's buffer and are invalidated by the next next().
struct FastqRecord {
    std::string_view name;     // without '@', up to the first whitespace
    std::string_view comment;  // whatever followed that whitespace, may be empty
    std::string_view seq;
    std::string_view qual;
};

// The value of a SAM-style tag in a FASTQ comment, or an empty view. Tags are TAB separated and
// the comment itself is whatever followed the first space, so a plain split on TAB is enough.
// This is the inter-stage contract (docs/formats.rst): refine, assemble and subsample all read
// the same tags out of the same comments, so they read them with the same function.
inline std::string_view tag_value(std::string_view comment, std::string_view key) {
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

class FastqReader {
public:
    explicit FastqReader(const std::string& path);
    ~FastqReader();
    FastqReader(FastqReader&&) noexcept;
    FastqReader& operator=(FastqReader&&) noexcept;

    // False at clean end of file. Throws MigecError on a malformed record or a truncated final
    // record -- silently dropping a partial record is how a run quietly loses its last reads.
    bool next(FastqRecord& out);

    uint64_t records_read() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class FastqWriter {
public:
    // gzip is chosen by a ".gz" suffix on `path`.
    explicit FastqWriter(const std::string& path, int gzip_level = 6);
    ~FastqWriter();
    FastqWriter(FastqWriter&&) noexcept;
    FastqWriter& operator=(FastqWriter&&) noexcept;

    // `comment` is written after a single space. Callers putting SAM tags there must use TABs
    // between the tags themselves -- that is what bwa -C and minimap2 -y require.
    void write(std::string_view name, std::string_view comment, std::string_view seq,
               std::string_view qual);
    void close();

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Bounded intake: stop after so many reads, or after so many distinct barcodes.
//
// This is for getting an answer out of a 400 GB run in a minute -- a smoke test, a pattern check,
// a look at the first thing that comes out. Never: it is NOT a sample. The first N reads of a
// FASTQ are one tile of one flowcell, and the first N barcodes are the barcodes that happen to
// sort early in the file, so nothing measured under a limit describes the library. `subsample`
// exists for that, and takes whole barcodes by hash so the MIG size distribution survives.
//
// The barcode set is bounded by `umis` itself -- the caller's own number -- so the memory this
// costs is the memory the caller asked for.
struct IntakeLimit {
    uint64_t reads = 0;  // 0 means no limit
    uint64_t umis = 0;

    bool active() const { return reads || umis; }
    // True while the stage should keep reading. `key` is the barcode this read carries.
    bool admit(uint64_t seen_reads, uint64_t key) {
        if (reads && seen_reads > reads) return false;
        if (!umis) return true;
        if (keys_.size() >= umis && keys_.find(key) == keys_.end()) return false;
        keys_.insert(key);
        return true;
    }

private:
    std::unordered_set<uint64_t> keys_;
};

// One FASTQ record, owning its bytes. The reader's views die on the next next(), so anything that
// hands a record to another thread has to own it.
struct FastqOwned {
    std::string name, comment, seq, qual;
};

// Appends one record in FASTQ form. Same layout FastqWriter emits, so a buffer built here and a
// file written there are byte-identical.
void append_fastq(std::string& dst, std::string_view name, std::string_view comment,
                  std::string_view seq, std::string_view qual);

// A complete gzip member for `in`, at `level`. Concatenated members are themselves a valid gzip
// stream (RFC 1952 s2.2), which is what lets every stage compress on its workers and leave the
// serial path with nothing to do but append bytes. zlib at its default level 6 manages ~7 MB/s on
// read payload against ~137 at level 1, so this is never called at 6 on a hot path.
void gzip_member(std::string_view in, std::string& out, int level);

}  // namespace migec

#endif  // MIGEC_FASTQ_HPP
