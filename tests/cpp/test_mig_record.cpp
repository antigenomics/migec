// The .mig format is the contract every stage depends on. These tests are the reason it can be
// changed safely: a round trip must be exact, and a truncated file must fail loudly rather than
// look like a short one.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <vector>

#include "migec/mig_record.hpp"
#include "migec/types.hpp"

using namespace migec;

namespace {

struct TempFile {
    std::string path;
    TempFile() {
        char tmpl[] = "/tmp/migec_test_XXXXXX";
        int fd = mkstemp(tmpl);
        if (fd >= 0) ::close(fd);
        path = tmpl;
    }
    ~TempFile() { std::remove(path.c_str()); }
};

struct Synth {
    std::vector<std::string> seq1, qual1, seq2, qual2;
    std::vector<MigRecord> recs;
};

// Deliberately includes zero-length reads: an empty mate is what a single-end file looks like,
// and an empty *record* is what a fully quality-trimmed read looks like.
Synth make_records(size_t n, uint32_t seed) {
    std::mt19937 rng(seed);
    Synth s;
    s.seq1.reserve(n); s.qual1.reserve(n); s.seq2.reserve(n); s.qual2.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        size_t l1 = (i % 97 == 0) ? 0 : 30 + rng() % 120;
        size_t l2 = (i % 89 == 0) ? 0 : 30 + rng() % 120;
        std::string a(l1, 'A'), qa(l1, 'I'), b(l2, 'C'), qb(l2, '5');
        for (size_t j = 0; j < l1; ++j) {
            a[j] = "ACGTN"[rng() % 5];
            qa[j] = char_from_phred(static_cast<uint8_t>(rng() % 42));
        }
        for (size_t j = 0; j < l2; ++j) {
            b[j] = "ACGT"[rng() % 4];
            qb[j] = char_from_phred(static_cast<uint8_t>(rng() % 42));
        }
        s.seq1.push_back(std::move(a)); s.qual1.push_back(std::move(qa));
        s.seq2.push_back(std::move(b)); s.qual2.push_back(std::move(qb));
    }
    s.recs.resize(n);
    for (size_t i = 0; i < n; ++i) {
        MigRecord& r = s.recs[i];
        r.cell = (static_cast<uint64_t>(rng()) << 32) ^ rng();
        r.umi = (static_cast<uint64_t>(rng()) << 32) ^ rng();
        r.src_index = i;
        r.flags = static_cast<uint16_t>(rng() % 256);
        r.umi_minq = static_cast<uint8_t>(rng() % 42);
        r.cell_minq = static_cast<uint8_t>(rng() % 42);
        r.seq1 = s.seq1[i]; r.qual1 = s.qual1[i];
        r.seq2 = s.seq2[i]; r.qual2 = s.qual2[i];
    }
    return s;
}

MigHeader make_header() {
    MigHeader h;
    h.umi_len = 12;
    h.cell_len = 16;
    h.bucket_index = 3;
    h.bucket_bits = 4;
    h.paired = true;
    h.sample_id = "S1";
    h.provenance = R"({"cmd":"migec checkout","version":"2.0.0.dev0"})";
    h.quality_calibration = {0.5f, 0.4f, 0.3f};
    return h;
}

}  // namespace

TEST_CASE("mig round trip is exact") {
    TempFile tf;
    const size_t n = 5000;
    Synth s = make_records(n, 42);
    const MigHeader h = make_header();

    {
        MigWriter w(tf.path, h, 64 << 10);  // small blocks, so the test crosses many boundaries
        for (const auto& r : s.recs) w.write(r);
        w.close();
        CHECK(w.records_written() == n);
    }

    MigReader r(tf.path);
    CHECK(r.header().umi_len == h.umi_len);
    CHECK(r.header().cell_len == h.cell_len);
    CHECK(r.header().bucket_index == h.bucket_index);
    CHECK(r.header().bucket_bits == h.bucket_bits);
    CHECK(r.header().paired == h.paired);
    CHECK(r.header().sample_id == h.sample_id);
    CHECK(r.header().provenance == h.provenance);
    REQUIRE(r.header().quality_calibration.size() == 3);
    CHECK(r.header().quality_calibration[1] == doctest::Approx(0.4f));

    MigRecord got;
    size_t i = 0;
    while (r.next(got)) {
        REQUIRE(i < n);
        const MigRecord& want = s.recs[i];
        CHECK(got.cell == want.cell);
        CHECK(got.umi == want.umi);
        CHECK(got.src_index == want.src_index);
        CHECK(got.flags == want.flags);
        CHECK(got.umi_minq == want.umi_minq);
        CHECK(got.cell_minq == want.cell_minq);
        CHECK(got.seq1 == want.seq1);
        CHECK(got.qual1 == want.qual1);
        CHECK(got.seq2 == want.seq2);
        CHECK(got.qual2 == want.qual2);
        ++i;
    }
    CHECK(i == n);
    CHECK(r.records_declared() == n);
}

TEST_CASE("empty file round trips") {
    TempFile tf;
    {
        MigWriter w(tf.path, make_header());
        w.close();
    }
    MigReader r(tf.path);
    MigRecord rec;
    CHECK_FALSE(r.next(rec));
    CHECK(r.records_declared() == 0);
}

TEST_CASE("truncation at every byte is an error, never a short read") {
    TempFile full;
    Synth s = make_records(400, 7);
    {
        MigWriter w(full.path, make_header(), 8 << 10);
        for (const auto& r : s.recs) w.write(r);
        w.close();
    }
    std::string blob;
    {
        std::FILE* f = std::fopen(full.path.c_str(), "rb");
        REQUIRE(f != nullptr);
        char buf[4096];
        size_t got;
        while ((got = std::fread(buf, 1, sizeof(buf), f)) > 0) blob.append(buf, got);
        std::fclose(f);
    }
    REQUIRE(blob.size() > 100);

    // Step over the file rather than testing all N cuts: enough boundaries to hit headers,
    // payloads and the footer without making the suite slow.
    const size_t step = blob.size() / 40 + 1;
    for (size_t cut = 1; cut < blob.size(); cut += step) {
        TempFile part;
        {
            std::FILE* f = std::fopen(part.path.c_str(), "wb");
            REQUIRE(f != nullptr);
            std::fwrite(blob.data(), 1, cut, f);
            std::fclose(f);
        }
        bool threw = false;
        size_t seen = 0;
        try {
            MigReader r(part.path);
            MigRecord rec;
            while (r.next(rec)) ++seen;
        } catch (const MigecError&) {
            threw = true;
        }
        // Truncating anywhere before the footer must be detected. The only way a truncated file
        // may read to completion is if the cut landed after the whole footer, which cannot
        // happen for cut < size.
        INFO("cut at " << cut << " of " << blob.size() << ", records read " << seen);
        CHECK(threw);
    }
}

TEST_CASE("a corrupt block body is caught by the CRC") {
    TempFile tf;
    Synth s = make_records(200, 11);
    {
        MigWriter w(tf.path, make_header(), 8 << 10);
        for (const auto& r : s.recs) w.write(r);
        w.close();
    }
    // Flip a byte late in the file, past the header and inside a compressed payload.
    std::FILE* f = std::fopen(tf.path.c_str(), "r+b");
    REQUIRE(f != nullptr);
    std::fseek(f, 0, SEEK_END);
    long size = std::ftell(f);
    std::fseek(f, size / 2, SEEK_SET);
    int c = std::fgetc(f);
    std::fseek(f, size / 2, SEEK_SET);
    std::fputc(c ^ 0xFF, f);
    std::fclose(f);

    bool threw = false;
    try {
        MigReader r(tf.path);
        MigRecord rec;
        while (r.next(rec)) {}
    } catch (const MigecError&) {
        threw = true;
    }
    CHECK(threw);
}

TEST_CASE("mismatched sequence and quality lengths are rejected on write") {
    TempFile tf;
    MigWriter w(tf.path, make_header());
    MigRecord r;
    r.seq1 = "ACGT";
    r.qual1 = "II";
    CHECK_THROWS_AS(w.write(r), MigecError);
}
