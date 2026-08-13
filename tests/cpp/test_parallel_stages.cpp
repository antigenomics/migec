// The parallel paths in refine and assemble, and the helper both of them run on.
//
// What is asserted here is the property that makes a thread count safe to change: the OUTPUT does
// not depend on it. Not "looks the same", not "the counts agree" -- the bytes. A demultiplexer or
// an assembler whose result shifts with -t produces results that cannot be compared between runs,
// and the failure is invisible in every summary.
//
// This file is also the target for the thread sanitizer:
//
//     cmake -S . -B build-tsan -DMIGEC_TESTS=ON -DCMAKE_BUILD_TYPE=RelWithDebInfo \
//           -DCMAKE_CXX_FLAGS="-fsanitize=thread -g"
//     cmake --build build-tsan -j && ./build-tsan/migec_tests -ts=parallel

#include "doctest.h"

#include <atomic>
#include <cstdio>
#include <fstream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include "migec/assemble.hpp"
#include "migec/parallel.hpp"
#include "migec/refine.hpp"
#include "migec/types.hpp"

using namespace migec;

namespace {

std::string temp_dir_path() {
    char tmpl[] = "/tmp/migec_par_XXXXXX";
    REQUIRE(::mkdtemp(tmpl) != nullptr);
    return std::string(tmpl);
}

// A tagged FASTQ in checkout's own output format: RX/QX/BC in the comment, payload in the read.
// Deliberately messy -- barcode errors, uneven MIG sizes, a few molecules deep enough to split --
// so the parallel paths run over groups of every shape rather than one.
std::string write_corpus(const std::string& dir, int molecules, unsigned seed) {
    const std::string path = dir + "/reads.fq";
    std::mt19937 rng(seed);
    std::ofstream out(path);
    const char* bases = "ACGT";
    uint64_t read_id = 0;
    for (int m = 0; m < molecules; ++m) {
        std::string umi;
        for (int i = 0; i < 12; ++i) umi += bases[rng() % 4];
        std::string payload;
        for (int i = 0; i < 60; ++i) payload += bases[rng() % 4];
        const int depth = 1 + static_cast<int>(rng() % 9);
        for (int r = 0; r < depth; ++r) {
            std::string seq = payload;
            if (rng() % 5 == 0) seq[rng() % seq.size()] = bases[rng() % 4];
            std::string tag = umi;
            if (rng() % 7 == 0) tag[rng() % tag.size()] = bases[rng() % 4];  // a barcode error
            out << "@r" << read_id++ << " RX:Z:" << tag << "\tQX:Z:" << std::string(12, 'I')
                << "\tBC:Z:S1\n"
                << seq << "\n+\n" << std::string(seq.size(), 'I') << "\n";
        }
    }
    return path;
}

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    std::ostringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

}  // namespace

TEST_SUITE("parallel") {

TEST_CASE("parallel_for runs every item exactly once") {
    for (int threads : {1, 2, 4, 9}) {
        const size_t n = 1000;
        std::vector<int> seen(n, 0);
        std::atomic<int> concurrent{0}, peak{0};
        parallel_for(n, threads, [&](size_t i, int worker) {
            const int now = concurrent.fetch_add(1) + 1;
            int was = peak.load();
            while (now > was && !peak.compare_exchange_weak(was, now)) {
            }
            CHECK(worker >= 0);
            CHECK(worker < threads);
            seen[i] += 1;
            concurrent.fetch_sub(1);
        });
        for (size_t i = 0; i < n; ++i) CHECK(seen[i] == 1);
        CHECK(peak.load() <= threads);
    }
}

TEST_CASE("an exception in a worker reaches the caller, and the lowest index wins") {
    // Whichever item throws FIRST in index order is the one reported, however the threads
    // happened to be scheduled -- otherwise the same broken input gives a different message on
    // every run and nothing is reproducible.
    for (int attempt = 0; attempt < 20; ++attempt) {
        bool threw = false;
        try {
            parallel_for(64, 8, [](size_t i, int) {
                if (i == 7 || i == 40) throw MigecError("boom " + std::to_string(i));
            });
        } catch (const MigecError& e) {
            threw = true;
            CHECK(std::string(e.what()) == "boom 7");
        }
        CHECK(threw);
    }
}

TEST_CASE("parallel_for on zero items does nothing and does not hang") {
    int calls = 0;
    parallel_for(0, 8, [&](size_t, int) { ++calls; });
    CHECK(calls == 0);
}

TEST_CASE("refine output is byte-identical at every thread count") {
    const std::string dir = temp_dir_path();
    const std::string reads = write_corpus(dir, 400, 11);

    std::string reference;
    RefineStats reference_stats;
    for (int threads : {1, 2, 3, 8}) {
        RefineRequest req;
        req.input = reads;
        req.output_dir = dir + "/refine" + std::to_string(threads);
        req.sample_id = "S1";
        req.correction.threads = threads;
        const RefineStats st = refine(req);
        const std::string bytes = slurp(req.output_dir + "/S1.fq.gz");
        REQUIRE(!bytes.empty());
        if (reference.empty()) {
            reference = bytes;
            reference_stats = st;
            continue;
        }
        CHECK(bytes == reference);
        // ...and the decisions themselves, not just the file: a merge that moved would change the
        // molecule count without necessarily changing the byte count.
        CHECK(st.merged == reference_stats.merged);
        CHECK(st.merged_reads == reference_stats.merged_reads);
        CHECK(st.molecules == reference_stats.molecules);
        CHECK(st.barcodes == reference_stats.barcodes);
        CHECK(st.estimated_error == doctest::Approx(reference_stats.estimated_error));
    }
}

TEST_CASE("assemble output is byte-identical at every thread count") {
    const std::string dir = temp_dir_path();
    const std::string reads = write_corpus(dir, 500, 23);

    std::string reference_fastq, reference_table;
    AssembleStats reference_stats;
    for (int threads : {1, 2, 5, 16}) {
        AssembleRequest req;
        req.input = reads;
        req.output_dir = dir + "/asm" + std::to_string(threads);
        req.sample_id = "S1";
        req.threads = threads;
        const AssembleStats st = assemble(req);
        const std::string fastq = slurp(req.output_dir + "/S1.consensus.fq.gz");
        const std::string table = slurp(req.output_dir + "/S1.mig.tsv");
        REQUIRE(!fastq.empty());
        REQUIRE(!table.empty());
        if (reference_fastq.empty()) {
            reference_fastq = fastq;
            reference_table = table;
            reference_stats = st;
            continue;
        }
        CHECK(fastq == reference_fastq);
        CHECK(table == reference_table);
        CHECK(st.groups == reference_stats.groups);
        CHECK(st.molecules == reference_stats.molecules);
        CHECK(st.groups_split == reference_stats.groups_split);
        CHECK(st.mean_quality == doctest::Approx(reference_stats.mean_quality));
    }
}

TEST_CASE("the bucket count does not depend on the thread count") {
    // The reason the bytes above match: if -t chose how finely the input was cut, it would also
    // choose the gzip member boundaries, and two runs would differ byte-wise while holding
    // identical records.
    const std::string dir = temp_dir_path();
    const std::string reads = write_corpus(dir, 100, 5);
    int reference_buckets = 0;
    for (int threads : {1, 4, 32}) {
        AssembleRequest req;
        req.input = reads;
        req.output_dir = dir + "/b" + std::to_string(threads);
        req.sample_id = "S1";
        req.threads = threads;
        const AssembleStats st = assemble(req);
        if (!reference_buckets) reference_buckets = st.buckets;
        CHECK(st.buckets == reference_buckets);
        CHECK(st.buckets >= (1 << kMinBucketBits));
    }
}

TEST_CASE("assembling the same input twice gives the same bytes") {
    // Determinism across runs, not just across thread counts: nothing may depend on an address, a
    // hash seed, or the order two temp files happened to be written in.
    const std::string dir = temp_dir_path();
    const std::string reads = write_corpus(dir, 300, 77);
    std::string first;
    for (int run = 0; run < 3; ++run) {
        AssembleRequest req;
        req.input = reads;
        req.output_dir = dir + "/run" + std::to_string(run);
        req.sample_id = "S1";
        const std::string bytes = (assemble(req), slurp(req.output_dir + "/S1.consensus.fq.gz"));
        if (first.empty()) first = bytes;
        CHECK(bytes == first);
    }
}

TEST_CASE("the intake limit stops at the read it says") {
    const std::string dir = temp_dir_path();
    const std::string reads = write_corpus(dir, 400, 31);

    AssembleRequest req;
    req.input = reads;
    req.output_dir = dir + "/limited";
    req.sample_id = "S1";
    req.limit.reads = 100;
    const AssembleStats st = assemble(req);
    CHECK(st.reads == 100);
    CHECK(st.limited);

    AssembleRequest whole;
    whole.input = reads;
    whole.output_dir = dir + "/whole";
    whole.sample_id = "S1";
    const AssembleStats all = assemble(whole);
    CHECK(all.reads > st.reads);
    CHECK_FALSE(all.limited);
}

TEST_CASE("the barcode limit counts barcodes, not reads") {
    const std::string dir = temp_dir_path();
    const std::string reads = write_corpus(dir, 400, 13);

    AssembleRequest req;
    req.input = reads;
    req.output_dir = dir + "/umilimit";
    req.sample_id = "S1";
    req.limit.umis = 50;
    const AssembleStats st = assemble(req);
    CHECK(st.limited);
    CHECK(st.groups <= 50);
    // Each of those barcodes brought its reads with it, so more reads than barcodes were taken.
    CHECK(st.reads > st.groups);
}

}  // TEST_SUITE
