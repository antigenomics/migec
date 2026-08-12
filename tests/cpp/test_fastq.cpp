#include "doctest.h"

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <string>

#include "migec/fastq.hpp"
#include "migec/types.hpp"

using namespace migec;

namespace {

struct TempFile {
    std::string path;
    explicit TempFile(const char* suffix = ".fq") {
        char tmpl[] = "/tmp/migec_test_XXXXXX";
        int fd = mkstemp(tmpl);
        if (fd >= 0) ::close(fd);
        std::remove(tmpl);  // the suffix makes a different name; gzopen needs to create it
        path = std::string(tmpl) + suffix;
    }
    ~TempFile() { std::remove(path.c_str()); }
};

void write_text(const std::string& path, const std::string& body) {
    std::FILE* f = std::fopen(path.c_str(), "wb");
    REQUIRE(f != nullptr);
    std::fwrite(body.data(), 1, body.size(), f);
    std::fclose(f);
}

}  // namespace

TEST_CASE("fastq parses names, comments and CRLF") {
    TempFile tf;
    write_text(tf.path,
               "@read1 some comment here\nACGT\n+\nIIII\n"
               "@read2\nGGGG\n+\n!!!!\n"
               "@read3\tafter a tab\r\nTTTT\r\n+\r\n####\r\n");
    FastqReader r(tf.path);
    FastqRecord rec;

    REQUIRE(r.next(rec));
    CHECK(rec.name == "read1");
    CHECK(rec.comment == "some comment here");
    CHECK(rec.seq == "ACGT");
    CHECK(rec.qual == "IIII");

    REQUIRE(r.next(rec));
    CHECK(rec.name == "read2");
    CHECK(rec.comment.empty());

    REQUIRE(r.next(rec));
    CHECK(rec.name == "read3");
    CHECK(rec.comment == "after a tab");
    CHECK(rec.seq == "TTTT");

    CHECK_FALSE(r.next(rec));
    CHECK(r.records_read() == 3);
}

TEST_CASE("a record straddling the read buffer survives") {
    // The reader refills in 1 MB chunks and compacts as it goes; a record that spans a refill is
    // where a view-based parser dangles. Write well over a chunk so many records straddle.
    TempFile tf;
    std::string body;
    const int n = 20000;
    for (int i = 0; i < n; ++i) {
        std::string seq(100 + (i % 37), "ACGT"[i % 4]);
        body += "@r" + std::to_string(i) + " c" + std::to_string(i) + "\n";
        body += seq + "\n+\n" + std::string(seq.size(), 'I') + "\n";
    }
    write_text(tf.path, body);

    FastqReader r(tf.path);
    FastqRecord rec;
    int i = 0;
    while (r.next(rec)) {
        REQUIRE(rec.name == "r" + std::to_string(i));
        REQUIRE(rec.comment == "c" + std::to_string(i));
        REQUIRE(rec.seq.size() == static_cast<size_t>(100 + (i % 37)));
        REQUIRE(rec.seq.size() == rec.qual.size());
        REQUIRE(rec.seq[0] == "ACGT"[i % 4]);
        ++i;
    }
    CHECK(i == n);
}

TEST_CASE("malformed records raise rather than silently stopping") {
    {
        TempFile tf;
        write_text(tf.path, "not a fastq at all\n");
        FastqReader r(tf.path);
        FastqRecord rec;
        CHECK_THROWS_AS(r.next(rec), MigecError);
    }
    {
        TempFile tf;
        write_text(tf.path, "@r1\nACGT\nX\nIIII\n");  // missing '+'
        FastqReader r(tf.path);
        FastqRecord rec;
        CHECK_THROWS_AS(r.next(rec), MigecError);
    }
    {
        TempFile tf;
        write_text(tf.path, "@r1\nACGT\n+\nII\n");  // quality shorter than sequence
        FastqReader r(tf.path);
        FastqRecord rec;
        CHECK_THROWS_AS(r.next(rec), MigecError);
    }
    {
        TempFile tf;
        write_text(tf.path, "@r1\nACGT\n+\n");  // truncated mid-record
        FastqReader r(tf.path);
        FastqRecord rec;
        CHECK_THROWS_AS(r.next(rec), MigecError);
    }
}

TEST_CASE("a missing final newline is tolerated") {
    TempFile tf;
    write_text(tf.path, "@r1\nACGT\n+\nIIII");
    FastqReader r(tf.path);
    FastqRecord rec;
    REQUIRE(r.next(rec));
    CHECK(rec.qual == "IIII");
    CHECK_FALSE(r.next(rec));
}

TEST_CASE("writer round trips through gzip") {
    TempFile tf(".fq.gz");
    {
        FastqWriter w(tf.path);
        w.write("r1", "CB:Z:AAAA\tRX:Z:CCCC", "ACGTACGT", "IIIIIIII");
        w.write("r2", "", "GG", "##");
        w.close();
    }
    FastqReader r(tf.path);
    FastqRecord rec;
    REQUIRE(r.next(rec));
    CHECK(rec.name == "r1");
    CHECK(rec.comment == "CB:Z:AAAA\tRX:Z:CCCC");  // TAB-separated tags survive intact
    CHECK(rec.seq == "ACGTACGT");
    REQUIRE(r.next(rec));
    CHECK(rec.name == "r2");
    CHECK(rec.comment.empty());
    CHECK_FALSE(r.next(rec));
}
