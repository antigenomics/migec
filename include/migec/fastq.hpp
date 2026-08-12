// FASTQ IO. Plain or gzipped, decided by the first two bytes rather than the file name, because
// half the world writes gzipped data to a name ending in .fastq.

#ifndef MIGEC_FASTQ_HPP
#define MIGEC_FASTQ_HPP

#include <memory>
#include <string>
#include <string_view>

namespace migec {

// A record's fields are views into the reader's buffer and are invalidated by the next next().
struct FastqRecord {
    std::string_view name;     // without '@', up to the first whitespace
    std::string_view comment;  // whatever followed that whitespace, may be empty
    std::string_view seq;
    std::string_view qual;
};

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

}  // namespace migec

#endif  // MIGEC_FASTQ_HPP
