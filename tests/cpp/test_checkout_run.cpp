// The whole-file driver: paired input, strand normalisation, and the property that matters most
// about the thread count -- that it changes nothing but the wall clock.

#include "doctest.h"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "migec/checkout.hpp"
#include "migec/fastq.hpp"
#include "migec/types.hpp"

using namespace migec;

namespace {

constexpr const char* kS1 = "aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN";
constexpr const char* kS2 = "aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN";

std::string temp_path(const std::string& suffix) {
    char tmpl[] = "/tmp/migec_run_XXXXXX";
    const int fd = ::mkstemp(tmpl);
    REQUIRE(fd >= 0);
    ::close(fd);
    std::string p(tmpl);
    std::remove(p.c_str());
    return p + suffix;
}

std::string tagged_read(const char* tag, const std::string& umi, const std::string& payload) {
    return std::string("aa") + tag + "cagtggtatcaacgcagagt" + umi.substr(0, 4) + "t" +
           umi.substr(4, 4) + "t" + umi.substr(8, 4) + payload;
}

std::string random_umi(std::mt19937& rng) {
    const char* B = "ACGT";
    std::string s;
    for (int i = 0; i < 12; ++i) s.push_back(B[rng() % 4]);
    return s;
}

std::string random_payload(std::mt19937& rng, size_t n) {
    const char* B = "ACGT";
    std::string s;
    for (size_t i = 0; i < n; ++i) s.push_back(B[rng() % 4]);
    return s;
}

// Writes a single-end FASTQ of `n` reads over both sample tags, and returns its path.
std::string write_reads(size_t n, uint32_t seed) {
    const std::string path = temp_path(".fq");
    FastqWriter w(path);
    std::mt19937 rng(seed);
    for (size_t i = 0; i < n; ++i) {
        const std::string seq =
            tagged_read(i % 2 ? "AGA" : "ACT", random_umi(rng), random_payload(rng, 60));
        w.write("r" + std::to_string(i), "", seq, std::string(seq.size(), 'I'));
    }
    w.close();
    return path;
}

std::string slurp(const std::string& path) {
    FastqReader r(path);
    FastqRecord rec;
    std::string out;
    while (r.next(rec)) {
        out += rec.name;
        out += '\t';
        out += rec.comment;
        out += '\t';
        out += rec.seq;
        out += '\n';
    }
    return out;
}

PatternSet two_samples() {
    PatternSet set;
    set.add("S1", kS1);
    set.add("S2", kS2);
    return set;
}

}  // namespace

TEST_CASE("thread count changes the wall clock and nothing else") {
    // If demultiplexing depended on -t, no two runs could be compared. Chunks are matched in
    // parallel but appended in input order, so this has to hold exactly, not approximately.
    const std::string in = write_reads(5000, 42);
    const PatternSet set = two_samples();

    std::vector<std::string> digests;
    for (int threads : {1, 3, 8}) {
        CheckoutRequest req;
        req.r1 = in;
        req.out_prefix = temp_path("_t" + std::to_string(threads) + ".");
        req.threads = threads;
        req.chunk_reads = 97;  // deliberately not a divisor of 5000
        CheckoutStats st = run_checkout(set, CheckoutParams{}, req);

        CHECK(st.counters.total == 5000);
        CHECK(st.counters.assigned == 5000);
        CHECK(st.threads == threads);
        digests.push_back(slurp(req.out_prefix + "S1.fq.gz") +
                          slurp(req.out_prefix + "S2.fq.gz"));
    }
    CHECK(digests[0] == digests[1]);
    CHECK(digests[0] == digests[2]);
    std::remove(in.c_str());
}

TEST_CASE("a run reports its own throughput and peak memory") {
    const std::string in = write_reads(500, 7);
    CheckoutRequest req;
    req.r1 = in;
    req.out_prefix = temp_path(".");
    CheckoutStats st = run_checkout(two_samples(), CheckoutParams{}, req);

    CHECK(st.wall_seconds > 0.0);
    CHECK(st.reads_per_second > 0.0);
    CHECK(st.peak_rss_bytes > 0);
    // The UMI counters are the allocation that scales with the library, so they are reported
    // separately from the process total.
    CHECK(st.umi_memory_bytes > 0);
    CHECK(st.umi_memory_bytes < st.peak_rss_bytes);
    std::remove(in.c_str());
}

TEST_CASE("paired input: the mate is carried through and tagged") {
    const std::string r1 = temp_path("_1.fq");
    const std::string r2 = temp_path("_2.fq");
    {
        FastqWriter w1(r1), w2(r2);
        const std::string seq = tagged_read("ACT", "AAAACCCCGGGG", "GACTCAGGGTTTCCAGGCCACAACTGCA");
        const std::string mate = "TTTTGGGGCCCCAAAATTTTGGGGCCCC";
        w1.write("p0", "", seq, std::string(seq.size(), 'I'));
        w2.write("p0", "", mate, std::string(mate.size(), 'I'));
        w1.close();
        w2.close();
    }
    CheckoutRequest req;
    req.r1 = r1;
    req.r2 = r2;
    req.out_prefix = temp_path(".");
    CheckoutStats st = run_checkout(two_samples(), CheckoutParams{}, req);
    CHECK(st.counters.assigned == 1);
    CHECK(st.counters.normalised == 0);

    FastqReader out1(req.out_prefix + "S1_R1.fq.gz");
    FastqReader out2(req.out_prefix + "S1_R2.fq.gz");
    FastqRecord a, b;
    REQUIRE(out1.next(a));
    REQUIRE(out2.next(b));
    CHECK(std::string(a.seq) == "GACTCAGGGTTTCCAGGCCACAACTGCA");
    // The mate is passed through whole -- there is no second tag to trim it by yet -- but it must
    // carry the barcode, or nothing downstream can group the pair.
    CHECK(std::string(b.seq) == "TTTTGGGGCCCCAAAATTTTGGGGCCCC");
    CHECK(std::string(b.comment).find("RX:Z:AAAACCCCGGGG") != std::string::npos);
    std::remove(r1.c_str());
    std::remove(r2.c_str());
}

TEST_CASE("a pair with the tag in R2 is swapped, not discarded") {
    // Amplicon libraries are sequenced in both orientations. A MIG holding both orientations of
    // one molecule silently loses half its reads at consensus, so checkout has to normalise here.
    const std::string r1 = temp_path("_1.fq");
    const std::string r2 = temp_path("_2.fq");
    const std::string tagged = tagged_read("ACT", "AAAACCCCGGGG", "GACTCAGGGTTTCCAGGCCACAACTGCA");
    {
        FastqWriter w1(r1), w2(r2);
        const std::string mate = "TTTTGGGGCCCCAAAATTTTGGGGCCCC";
        w1.write("p0", "", mate, std::string(mate.size(), 'I'));   // untagged mate first
        w2.write("p0", "", tagged, std::string(tagged.size(), 'I'));
        w1.close();
        w2.close();
    }
    CheckoutRequest req;
    req.r1 = r1;
    req.r2 = r2;
    req.out_prefix = temp_path(".");
    CheckoutStats st = run_checkout(two_samples(), CheckoutParams{}, req);
    CHECK(st.counters.assigned == 1);
    CHECK(st.counters.normalised == 1);

    FastqReader out1(req.out_prefix + "S1_R1.fq.gz");
    FastqRecord a;
    REQUIRE(out1.next(a));
    CHECK(std::string(a.seq) == "GACTCAGGGTTTCCAGGCCACAACTGCA");  // R1 always carries the tag
    std::remove(r1.c_str());
    std::remove(r2.c_str());
}

TEST_CASE("a single-end read carrying the tag on the other strand is flipped") {
    const std::string tagged = tagged_read("ACT", "AAAACCCCGGGG", "GACTCAGGGTTTCCAGGCCACAACTGCA");
    std::string rc = tagged, rq(tagged.size(), 'I');
    reverse_complement(rc, rq);

    const std::string in = temp_path(".fq");
    {
        FastqWriter w(in);
        w.write("r0", "", rc, rq);
        w.close();
    }
    CheckoutRequest req;
    req.r1 = in;
    req.out_prefix = temp_path(".");
    CheckoutStats st = run_checkout(two_samples(), CheckoutParams{}, req);
    CHECK(st.counters.assigned == 1);
    CHECK(st.counters.normalised == 1);

    FastqReader out(req.out_prefix + "S1.fq.gz");
    FastqRecord a;
    REQUIRE(out.next(a));
    CHECK(std::string(a.seq) == "GACTCAGGGTTTCCAGGCCACAACTGCA");
    std::remove(in.c_str());
}

TEST_CASE("truncated R2 is an error, not a silent half-run") {
    const std::string r1 = temp_path("_1.fq");
    const std::string r2 = temp_path("_2.fq");
    {
        FastqWriter w1(r1), w2(r2);
        const std::string s = tagged_read("ACT", "AAAACCCCGGGG", "GACTCAGGGTTTCCAGG");
        w1.write("p0", "", s, std::string(s.size(), 'I'));
        w1.write("p1", "", s, std::string(s.size(), 'I'));
        w2.write("p0", "", s, std::string(s.size(), 'I'));
        w1.close();
        w2.close();
    }
    CheckoutRequest req;
    req.r1 = r1;
    req.r2 = r2;
    req.out_prefix = temp_path(".");
    CHECK_THROWS_AS(run_checkout(two_samples(), CheckoutParams{}, req), MigecError);
    std::remove(r1.c_str());
    std::remove(r2.c_str());
}

TEST_CASE("two rows with one sample id are one sample, not two handles on one file") {
    // A MIGEC barcode table declares a sample sequenced with more than one tag as several rows
    // sharing the id. Opening a file per row means two FILE* on one path, whose interleaved
    // writes are not even a valid gzip stream -- and the summary would still report success.
    const std::string in = write_reads(2000, 11);
    PatternSet set;
    set.add("S1", kS1);
    set.add("S1", kS2);  // same sample, second tag

    CheckoutRequest req;
    req.r1 = in;
    req.out_prefix = temp_path(".");
    req.chunk_reads = 64;
    CheckoutStats st = run_checkout(set, CheckoutParams{}, req);

    REQUIRE(st.sample_ids.size() == 1);
    CHECK(st.sample_ids[0] == "S1");
    CHECK(st.sample_reads[0] == 2000);
    // Both tags' reads land in one file, readable, with none lost.
    size_t n = 0;
    FastqReader r(req.out_prefix + "S1.fq.gz");
    FastqRecord rec;
    while (r.next(rec)) ++n;
    CHECK(n == 2000);
    // ...and in one UMI counter, so the collision statistics describe the sample rather than a tag.
    CHECK(st.umi_counts.size() == 1);
    CHECK(st.umi_counts[0].total() == 2000);
    std::remove(in.c_str());
}

TEST_CASE("a sample that got no reads is still a readable empty file") {
    // Zero bytes is not a gzip stream: `gzip -t` and `zcat` reject it. Skipping the per-chunk
    // empty members is right, but the file still needs one.
    const std::string in = write_reads(200, 12);
    PatternSet set;
    set.add("S1", kS1);
    set.add("S3", "aaTTTcagtggtatcaacgcagagtNNNNtNNNNtNNNN");  // a tag no read carries

    CheckoutRequest req;
    req.r1 = in;
    req.out_prefix = temp_path(".");
    CheckoutStats st = run_checkout(set, CheckoutParams{}, req);
    CHECK(st.sample_reads[1] == 0);

    FastqReader r(req.out_prefix + "S3.fq.gz");  // throws if the file is not a gzip stream
    FastqRecord rec;
    CHECK_FALSE(r.next(rec));
    std::remove(in.c_str());
}

TEST_CASE("a UMI longer than the packed representation is rejected at the pattern, not the read") {
    // Otherwise it throws inside pack_barcode on a worker thread, where an escaping exception
    // calls std::terminate: an abort with no message rather than an error naming the bad row.
    CHECK_THROWS_AS(BarcodePattern::compile(std::string(33, 'N') + "ACGT"), MigecError);
    CHECK_NOTHROW(BarcodePattern::compile(std::string(32, 'N') + "ACGT"));
}
