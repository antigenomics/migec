#include "doctest.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <random>
#include <string>

#include "migec/types.hpp"
#include "migec/umi_stats.hpp"

using namespace migec;

namespace {

uint64_t umi(const std::string& s) { return pack_barcode(s); }

}  // namespace

TEST_CASE("coverage histogram bins by powers of two") {
    UmiCounts c(4);
    c.add(umi("AAAA"), 1);
    c.add(umi("AAAC"), 1);
    c.add(umi("AAAG"), 3);   // bin 1: [2,4)
    c.add(umi("AAAT"), 8);   // bin 3: [8,16)
    c.add(umi("AACA"), 100); // bin 6: [64,128)

    CoverageHistogram h = c.histogram();
    CHECK(h.units[0] == 2);
    CHECK(h.reads[0] == 2);
    CHECK(h.units[1] == 1);
    CHECK(h.reads[1] == 3);
    CHECK(h.units[3] == 1);
    CHECK(h.units[6] == 1);
    CHECK(h.total_reads() == 113);
    CHECK(h.total_units() == 5);
    CHECK(h.mean_reads_per_umi() == doctest::Approx(113.0 / 5.0));

    // 108 of 113 reads sit in MIGs of >= 5 reads.
    CHECK(h.reads_in_migs_at_least(5) == doctest::Approx(108.0 / 113.0));
    CHECK(h.over_sequenced());
}

TEST_CASE("a uniform UMI has full entropy and its real length") {
    // Every 4^3 = 64 barcode of length 3, once each: perfectly uniform by construction.
    UmiCounts c(3);
    const char* B = "ACGT";
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k) {
                std::string s = {B[i], B[j], B[k]};
                c.add(umi(s), 1);
            }

    UmiComposition comp = c.composition(false);
    CHECK(comp.length == 3);
    for (int j = 0; j < 3; ++j) {
        CHECK(comp.entropy(j) == doctest::Approx(2.0));
        CHECK(comp.information(j) == doctest::Approx(0.0));
        CHECK(comp.collision(j) == doctest::Approx(0.25));
    }
    CHECK(comp.total_entropy() == doctest::Approx(6.0));
    CHECK(comp.total_information() == doctest::Approx(0.0));
    CHECK(comp.effective_length() == doctest::Approx(3.0));
    CHECK(comp.effective_space() == doctest::Approx(64.0));
}

TEST_CASE("a skewed UMI loses effective length, and Shannon overstates the space") {
    // Position 0 is fixed to A, the rest uniform: one base of the barcode is dead.
    UmiCounts c(3);
    const char* B = "ACGT";
    for (int j = 0; j < 4; ++j)
        for (int k = 0; k < 4; ++k) {
            std::string s = {'A', B[j], B[k]};
            c.add(umi(s), 1);
        }
    UmiComposition comp = c.composition(false);
    CHECK(comp.entropy(0) == doctest::Approx(0.0));
    CHECK(comp.collision(0) == doctest::Approx(1.0));
    CHECK(comp.effective_length() == doctest::Approx(2.0));   // 3 nominal, 2 usable
    CHECK(comp.effective_space() == doctest::Approx(16.0));

    // Now a mildly skewed position: Shannon-derived space (2^H) must EXCEED the collision-derived
    // space, which is the whole reason the collision form is the one used for decisions.
    UmiCounts s(1);
    s.add(umi("A"), 70);
    s.add(umi("C"), 10);
    s.add(umi("G"), 10);
    s.add(umi("T"), 10);
    UmiComposition sc = s.composition(true);
    const double shannon_space = std::pow(2.0, sc.total_entropy());
    CHECK(shannon_space > sc.effective_space());
}

TEST_CASE("expected collisions follows the birthday bound") {
    UmiCounts c(6);
    const char* B = "ACGT";
    // Uniform-ish sample over 6 positions.
    std::mt19937 rng(7);
    for (int i = 0; i < 4000; ++i) {
        std::string s;
        for (int j = 0; j < 6; ++j) s.push_back(B[rng() % 4]);
        c.add(umi(s), 1);
    }
    UmiComposition comp = c.composition(false);
    // 4^6 = 4096 barcodes; with 100 molecules E[colliding pairs] = 100^2/2 / 4096 ~ 1.22
    CHECK(comp.expected_collisions(100.0) == doctest::Approx(1.22).epsilon(0.15));
}

namespace {

// A realistic background library: a decision about one barcode depends on how many other
// molecules there are and how their sizes are distributed, so testing correction on three
// barcodes tests nothing.
UmiCounts background_library(int n_umis, int umi_len, uint32_t seed) {
    UmiCounts c(umi_len);
    std::mt19937 rng(seed);
    std::lognormal_distribution<double> size(2.0, 0.8);  // median ~7 reads, long tail
    const char* B = "ACGT";
    for (int i = 0; i < n_umis; ++i) {
        std::string s;
        for (int j = 0; j < umi_len; ++j) s.push_back(B[rng() % 4]);
        c.add(umi(s), static_cast<uint32_t>(std::max(1.0, size(rng))));
    }
    return c;
}

// CorrectionResult is indexed in parallel with counts.entries() rather than keyed by barcode --
// 12 bytes per UMI instead of two hash maps' worth. These translate back to barcodes for the
// assertions.
bool merged_away(const UmiCounts& c, const CorrectionResult& r, const std::string& s) {
    const size_t i = index_of(c, umi(s));
    return i != static_cast<size_t>(-1) && r.root[i] != i;
}
uint64_t root_of(const UmiCounts& c, const CorrectionResult& r, const std::string& s) {
    return c.entries()[r.root[index_of(c, umi(s))]].key;
}
uint32_t reads_of(const UmiCounts& c, const CorrectionResult& r, const std::string& s) {
    return r.corrected[index_of(c, umi(s))];
}

}  // namespace

TEST_CASE("a small child of a large parent is merged") {
    // The unambiguous case: a barcode with a couple of reads sitting one substitution from a
    // deeply covered molecule is what a UMI sequencing error looks like.
    UmiCounts c = background_library(2000, 12, 5);
    c.add(umi("ACGTACGTACGT"), 10000);
    c.add(umi("ACGTACGTACGA"), 2);

    CorrectionParams p;
    p.sequencing_error = 1e-3;
    CorrectionResult r = correct_umis(c, p);

    REQUIRE(merged_away(c, r, "ACGTACGTACGA"));
    CHECK(root_of(c, r, "ACGTACGTACGA") == umi("ACGTACGTACGT"));
    CHECK(reads_of(c, r, "ACGTACGTACGT") == 10002);
    CHECK(reads_of(c, r, "ACGTACGTACGA") == 0);
}

TEST_CASE("a neighbour of comparable size is not merged") {
    // Two molecules of similar abundance one substitution apart are far more likely to be two real
    // molecules than a parent and its error child -- no error turns 10000 reads into 9000.
    UmiCounts c = background_library(2000, 12, 6);
    c.add(umi("ACGTACGTACGT"), 10000);
    c.add(umi("ACGTACGTACGA"), 9000);

    CorrectionParams p;
    p.sequencing_error = 1e-3;
    CorrectionResult r = correct_umis(c, p);
    CHECK(!merged_away(c, r, "ACGTACGTACGA"));
}

TEST_CASE("an isolated low-coverage UMI keeps its reads") {
    // The explicit requirement: a molecule seen 3-5 times with NO plausible parent is information,
    // not noise. It must survive correction untouched.
    UmiCounts c = background_library(2000, 12, 7);
    c.add(umi("TTTTTTTTTTTT"), 4);  // 12 substitutions from anything else, by construction

    CorrectionResult r = correct_umis(c);
    CHECK(!merged_away(c, r, "TTTTTTTTTTTT"));
    CHECK(reads_of(c, r, "TTTTTTTTTTTT") == 4);
}

TEST_CASE("a child is never merged into a smaller or equal barcode") {
    UmiCounts c = background_library(500, 12, 8);
    c.add(umi("ACGTACGTACGT"), 50);
    c.add(umi("ACGTACGTACGA"), 50);
    CorrectionResult r = correct_umis(c);
    CHECK(!merged_away(c, r, "ACGTACGTACGT"));
    CHECK(!merged_away(c, r, "ACGTACGTACGA"));
}

TEST_CASE("chains resolve to a root, never to a cycle") {
    UmiCounts c(8);
    c.add(umi("ACGTACGT"), 100000);
    c.add(umi("ACGTACGA"), 300);   // child of the above
    c.add(umi("ACGTACAA"), 5);     // child of the child

    CorrectionParams p;
    p.sequencing_error = 1e-2;
    CorrectionResult r = correct_umis(c, p);

    // Whatever merges, every surviving barcode must be a root and the reads must be conserved.
    uint64_t total = 0;
    for (size_t i = 0; i < r.corrected.size(); ++i) {
        if (r.corrected[i] == 0) continue;
        CHECK(r.root[i] == i);  // every surviving barcode is its own root
        total += r.corrected[i];
    }
    CHECK(total == 100305);
}

TEST_CASE("reads are conserved by correction, always") {
    std::mt19937 rng(11);
    UmiCounts c(6);
    const char* B = "ACGT";
    uint64_t expected = 0;
    for (int i = 0; i < 2000; ++i) {
        std::string s;
        for (int j = 0; j < 6; ++j) s.push_back(B[rng() % 4]);
        const uint32_t n = 1 + rng() % 50;
        c.add(umi(s), n);
        expected += n;
    }
    CorrectionResult r = correct_umis(c);
    uint64_t total = 0;
    for (uint32_t n : r.corrected) total += n;
    CHECK(total == c.total());
    CHECK(c.total() <= expected);  // duplicate draws merged at add() time
}

TEST_CASE("the error rate estimate recovers an injected rate") {
    // Build parents, then add a 1-substitution child to a known fraction of them.
    std::mt19937 rng(3);
    UmiCounts c(10);
    const char* B = "ACGT";
    const int n_parents = 3000;
    const double eps = 2e-3;
    const int parent_size = 200;  // expected children per parent = 3L(1-e^{-c eps}) ~ 9.9

    std::vector<std::string> parents;
    for (int i = 0; i < n_parents; ++i) {
        std::string s;
        for (int j = 0; j < 10; ++j) s.push_back(B[rng() % 4]);
        parents.push_back(s);
        c.add(umi(s), parent_size);
    }
    std::binomial_distribution<int> nchild(30, 1.0 - std::exp(-parent_size * eps));
    for (const auto& s : parents) {
        const int k = nchild(rng);
        for (int t = 0; t < k; ++t) {
            std::string ch = s;
            const int pos = static_cast<int>(rng() % 10);
            char nb = B[rng() % 4];
            while (nb == ch[static_cast<size_t>(pos)]) nb = B[rng() % 4];
            ch[static_cast<size_t>(pos)] = nb;
            c.add(umi(ch), 1);
        }
    }

    UmiComposition comp = c.composition(false);
    const double est = estimate_umi_error(c, comp);
    INFO("estimated " << est << " vs injected " << eps);
    CHECK(est > eps / 3.0);
    CHECK(est < eps * 3.0);
}

TEST_CASE("a half-full barcode space is flagged and collision-corrected") {
    // 128 distinct barcodes drawn at random from the 256 of length 4. MIGEC disabled correction
    // outright in this regime; here it keeps running -- the collision prior makes it
    // self-limiting -- and the molecule count is corrected upward for the collisions that no
    // method can see.
    UmiCounts c(4);
    const char* B = "ACGT";
    std::mt19937 rng(21);
    while (c.distinct() < 128) {
        std::string s;
        for (int j = 0; j < 4; ++j) s.push_back(B[rng() % 4]);
        c.add(umi(s), 10);
    }
    CorrectionResult r = correct_umis(c);
    CHECK(r.saturated);
    CHECK(r.merged == 0);  // every barcode sits at equal depth; nothing looks like a child
    // 256 * -ln(1 - 128/256) = 177
    CHECK(r.molecules_corrected > static_cast<double>(r.molecules_observed));
    CHECK(r.molecules_corrected == doctest::Approx(177.0).epsilon(0.15));
}

TEST_CASE("an over-full barcode space declines to estimate rather than reporting zero collisions") {
    // Occupancy above 90% makes the space estimate collapse onto the observed count, which would
    // report "no collisions" for the most collided library possible. Refuse instead.
    UmiCounts c(4);
    const char* B = "ACGT";
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                for (int l = 0; l < 4; ++l) c.add(umi({B[i], B[j], B[k], B[l]}), 10);
    CorrectionResult r = correct_umis(c);
    CHECK(r.saturated);
    CHECK(r.molecules_corrected == doctest::Approx(static_cast<double>(r.molecules_observed)));
}

TEST_CASE("the flush path gives the same answer as one big buffer") {
    // Counts are accumulated into a bounded buffer that is sorted and folded into a sorted array
    // whenever it fills. A run whose buffer flushes a hundred times must be indistinguishable from
    // one that never flushes -- otherwise the answer depends on how much RAM was available.
    std::mt19937 rng(31);
    std::vector<uint64_t> draws;
    const char* B = "ACGT";
    for (int i = 0; i < 20000; ++i) {
        std::string s;
        for (int j = 0; j < 8; ++j) s.push_back(B[rng() % 4]);
        draws.push_back(umi(s));
    }

    UmiCounts big(8, 1u << 20);   // never flushes until the end
    UmiCounts small(8, 64);       // flushes ~300 times
    for (uint64_t d : draws) {
        big.add(d);
        small.add(d);
    }

    REQUIRE(big.distinct() == small.distinct());
    CHECK(big.total() == small.total());
    CHECK(big.total() == 20000);
    const std::vector<UmiCounts::Entry>& a = big.entries();
    const std::vector<UmiCounts::Entry>& b = small.entries();
    for (size_t i = 0; i < a.size(); ++i) {
        CHECK(a[i].key == b[i].key);
        CHECK(a[i].count == b[i].count);
    }
    // Sorted, which is what the range partition and the neighbourhood search both rely on.
    for (size_t i = 1; i < a.size(); ++i) CHECK(a[i - 1].key < a[i].key);
}

TEST_CASE("merging two counters is the same as counting once") {
    UmiCounts x(6), y(6), both(6);
    std::mt19937 rng(17);
    const char* B = "ACGT";
    for (int i = 0; i < 3000; ++i) {
        std::string s;
        for (int j = 0; j < 6; ++j) s.push_back(B[rng() % 4]);
        (i % 2 ? y : x).add(umi(s));
        both.add(umi(s));
    }
    x.merge(y);
    REQUIRE(x.distinct() == both.distinct());
    CHECK(x.total() == both.total());
    for (size_t i = 0; i < both.entries().size(); ++i) {
        CHECK(x.entries()[i].key == both.entries()[i].key);
        CHECK(x.entries()[i].count == both.entries()[i].count);
    }
}

TEST_CASE("the counter costs a small constant per distinct UMI") {
    // The number this class exists for. A hash map of the same contents runs to ~48 bytes per
    // entry once nodes and the bucket array are counted, which is the difference between a run
    // fitting in memory and not.
    UmiCounts c(12, 4096);
    std::mt19937 rng(23);
    const char* B = "ACGT";
    for (int i = 0; i < 100000; ++i) {
        std::string s;
        for (int j = 0; j < 12; ++j) s.push_back(B[rng() % 4]);
        c.add(umi(s));
    }
    const double per_umi = static_cast<double>(c.memory_bytes()) / static_cast<double>(c.distinct());
    CHECK(per_umi < 32.0);
}

namespace {

// A diverse background, because everything in correction is measured against the composition of
// the library it is in: a background of near-identical barcodes collapses the effective space and
// the whole library reads as saturated.
std::vector<std::string> background(size_t n, int len, uint32_t seed) {
    std::mt19937 rng(seed);
    std::vector<std::string> out;
    while (out.size() < n) {
        std::string s(static_cast<size_t>(len), 'A');
        for (char& c : s) c = "ACGT"[rng() & 3u];
        out.push_back(s);
    }
    return out;
}

// Fills `counts` and returns an evidence block whose payloads are per-barcode unique unless the
// barcode appears in `shared`, which all get the same sequence.
BarcodeEvidence make_evidence(const UmiCounts& counts, int len, int width,
                             const std::vector<std::string>& shared) {
    BarcodeEvidence ev;
    ev.payload_width = width;
    ev.payload.assign(counts.entries().size() * static_cast<size_t>(width), 'A');
    for (size_t i = 0; i < counts.entries().size(); ++i) {
        const std::string b = unpack_barcode(counts.entries()[i].key, len);
        const bool is_shared =
            std::find(shared.begin(), shared.end(), b) != shared.end();
        // The packed key fills from the HIGH bits, so the low 32 are zero for a 12 nt barcode --
        // seeding from them would give every payload the same sequence.
        const uint64_t k = counts.entries()[i].key;
        std::mt19937 rng(is_shared ? 42u : static_cast<uint32_t>((k >> 32) ^ k));
        for (int j = 0; j < width; ++j) {
            ev.payload[i * static_cast<size_t>(width) + static_cast<size_t>(j)] =
                "ACGT"[rng() & 3u];
        }
    }
    return ev;
}

}  // namespace

TEST_CASE("payload evidence merges two singletons the count ratio cannot") {
    // The shallow case: a barcode and its error child, one read each. No count asymmetry exists,
    // so the count-ratio rule can never fire -- but both reads are reads of the same molecule and
    // say so.
    UmiCounts c(12);
    c.add(umi("ACGTACGTACGT"), 1);
    c.add(umi("ACGTACGTACGA"), 1);  // distance 1
    for (const std::string& b : background(3000, 12, 7)) c.add(umi(b), 1);

    CorrectionParams p;
    CHECK(correct_umis(c, p).merged == 0);  // no evidence without the reads

    const BarcodeEvidence ev =
        make_evidence(c, 12, 60, {"ACGTACGTACGT", "ACGTACGTACGA"});
    const CorrectionResult r = correct_umis(c, p, ev);
    CHECK(r.merged == 1);
    CHECK(r.merged_by_payload == 1);
    CHECK(r.payload_clonality < 0.1);
}

TEST_CASE("a disagreeing payload refuses a merge the counts would have made") {
    // 100 reads against 1 at distance 1 is exactly what the count ratio calls a child. It is not
    // one: the reads are of a different molecule, and that has to win.
    UmiCounts c(12);
    c.add(umi("ACGTACGTACGT"), 100);
    c.add(umi("ACGTACGTACGA"), 1);
    for (const std::string& b : background(3000, 12, 8)) c.add(umi(b), 1);

    CorrectionParams p;
    CHECK(correct_umis(c, p).merged >= 1);  // counts alone merge it

    const BarcodeEvidence ev = make_evidence(c, 12, 60, {});  // every payload its own
    const CorrectionResult r = correct_umis(c, p, ev);
    const size_t child = index_of(c, umi("ACGTACGTACGA"));
    CHECK(r.root[child] == child);  // not merged into anything
}

TEST_CASE("a clonal library gets no help from payload agreement, and says so") {
    // Every molecule carries the same sequence, so agreement is worth nothing -- and the measured
    // clonality is what says that, rather than the evidence being silently over-trusted.
    UmiCounts c(12);
    c.add(umi("ACGTACGTACGT"), 1);
    c.add(umi("ACGTACGTACGA"), 1);
    std::vector<std::string> all{"ACGTACGTACGT", "ACGTACGTACGA"};
    for (const std::string& b : background(3000, 12, 9)) {
        c.add(umi(b), 1);
        all.push_back(b);
    }
    const BarcodeEvidence ev = make_evidence(c, 12, 60, all);  // one clone, all identical
    const CorrectionResult r = correct_umis(c, CorrectionParams{}, ev);
    CHECK(r.payload_clonality > 0.9);
    CHECK(r.merged_by_payload == 0);
}

TEST_CASE("barcode base quality sharpens the error prior") {
    // The same pair, once with a high-quality mismatching base and once with a low-quality one.
    // Only the second is a plausible miscall, and the posterior has to see the difference.
    auto merged_with = [](float err_at_last) {
        UmiCounts c(12);
        c.add(umi("ACGTACGTACGT"), 20);
        c.add(umi("ACGTACGTACGA"), 1);
        for (const std::string& b : background(3000, 12, 10)) c.add(umi(b), 1);
        BarcodeEvidence ev;
        ev.position_error.assign(c.entries().size() * 12, 1e-6f);
        for (size_t i = 0; i < c.entries().size(); ++i) {
            if (unpack_barcode(c.entries()[i].key, 12) == "ACGTACGTACGA") {
                ev.position_error[i * 12 + 11] = err_at_last;
            }
        }
        CorrectionParams p;
        p.sequencing_error = 1e-6;   // fixed, so the only difference is the per-base term
        // ...and no polymerase component, which is the hypothesis that does NOT depend on the
        // reported quality: an early-PCR child carries a high Phred in every read, so leaving it
        // in would merge the pair on count asymmetry alone and hide what is being tested.
        p.polymerase_error = 0.0;
        const CorrectionResult r = correct_umis(c, p, ev);
        const size_t child = index_of(c, umi("ACGTACGTACGA"));
        return r.root[child] != child;
    };
    CHECK_FALSE(merged_with(1e-6f));  // a confident base is not a miscall
    CHECK(merged_with(0.3f));         // Q5 at exactly the base that differs
}

// --- spilling: the range partition that bounds the counters (roadmap item 1) -----------------

TEST_CASE("a spilled counter answers exactly as a resident one") {
    // The whole claim: partitioning changes the memory, never the numbers. Both counters see the
    // same adds in the same order; one is forced to spill many times over.
    std::mt19937_64 rng(12345);
    std::vector<std::pair<uint64_t, uint32_t>> adds;
    for (int i = 0; i < 20000; ++i) {
        std::string s;
        for (int j = 0; j < 8; ++j) s += "ACGT"[rng() % 4];
        adds.emplace_back(umi(s), static_cast<uint32_t>(1 + rng() % 5));
    }

    UmiCounts resident(8);
    UmiCounts spilling(8);
    const std::string dir =
        (std::filesystem::temp_directory_path() / "migec_spill_test").string();
    std::filesystem::remove_all(dir);
    // A tiny budget so it spills repeatedly: a key ends up written to its bucket many times over,
    // which is the case reduction-on-read has to get right.
    spilling.enable_spill(dir, 4096, 4);
    for (const auto& a : adds) {
        resident.add(a.first, a.second);
        spilling.add(a.first, a.second);
    }

    CHECK(spilling.spilled());
    CHECK(spilling.total() == resident.total());
    CHECK(spilling.distinct() == resident.distinct());

    const CoverageHistogram hr = resident.histogram();
    const CoverageHistogram hs = spilling.histogram();
    CHECK(hs.reads == hr.reads);
    CHECK(hs.units == hr.units);
    CHECK(hs.mean_reads_per_umi() == doctest::Approx(hr.mean_reads_per_umi()));

    const UmiComposition cr = resident.composition(false);
    const UmiComposition cs = spilling.composition(false);
    REQUIRE(cs.length == cr.length);
    for (int j = 0; j < cr.length; ++j) {
        for (int b = 0; b < 4; ++b) {
            CHECK(cs.freq[static_cast<size_t>(j)][static_cast<size_t>(b)] ==
                  doctest::Approx(cr.freq[static_cast<size_t>(j)][static_cast<size_t>(b)]));
        }
    }
    CHECK(cs.effective_length() == doctest::Approx(cr.effective_length()));

    // for_each must still deliver ascending keys, each exactly once: the range partition is
    // ordered by construction and everything downstream assumes it.
    std::vector<uint64_t> keys;
    spilling.for_each([&keys](const UmiCounts::Entry& e) { keys.push_back(e.key); });
    CHECK(keys.size() == resident.distinct());
    CHECK(std::is_sorted(keys.begin(), keys.end()));
    CHECK(std::adjacent_find(keys.begin(), keys.end()) == keys.end());

    std::filesystem::remove_all(dir);
}

TEST_CASE("a spilled counter refuses the accessors that need every entry") {
    // Never: silently returning the resident remainder would be a wrong answer wearing the shape
    // of a right one -- `find()` would report "absent" for a barcode that is merely spilled.
    const std::string dir =
        (std::filesystem::temp_directory_path() / "migec_spill_refuse").string();
    std::filesystem::remove_all(dir);
    UmiCounts c(8);
    c.enable_spill(dir, 64, 3);
    // Past the internal flush threshold: a spill happens on flush, so a handful of adds sit in the
    // buffer and never reach the partition at all.
    std::mt19937_64 rng(7);
    for (int i = 0; i < 20000; ++i) {
        std::string s;
        for (int j = 0; j < 8; ++j) s += "ACGT"[rng() % 4];
        c.add(umi(s), 1);
    }
    c.distinct();  // forces the flush, and with it the spill
    REQUIRE(c.spilled());
    CHECK_THROWS_AS(c.entries(), MigecError);
    CHECK_THROWS_AS(c.find(umi("AAAAAAAA")), MigecError);
    // ...but the streaming ones still work, which is the point.
    CHECK(c.distinct() > 0);
    CHECK(c.histogram().units[0] > 0);
    std::filesystem::remove_all(dir);
}

TEST_CASE("enable_spill validates its arguments where they are attributable") {
    UmiCounts c(4);
    CHECK_THROWS_AS(c.enable_spill("/tmp/migec_bad", 1024, 0), MigecError);
    CHECK_THROWS_AS(c.enable_spill("/tmp/migec_bad", 1024, 21), MigecError);
    CHECK_THROWS_AS(c.enable_spill("/tmp/migec_bad", 0, 4), MigecError);
    // The rotated second correction pass needs a second prefix's worth of barcode, so a partition
    // wider than half the barcode is refused where the number is attributable -- not later, as a
    // merge count that is quietly short.
    CHECK_THROWS_AS(c.enable_spill("/tmp/migec_bad", 1024, 5), MigecError);
    UmiCounts wide(12);
    wide.enable_spill((std::filesystem::temp_directory_path() / "migec_bits_ok").string(), 1024, 8);
}

// --- bucketed correction: the rotated second pass (roadmap item 2) ---------------------------

TEST_CASE("rotate_barcode is a rotation of the barcode, not of the word") {
    const uint64_t k = pack_barcode("ACGTACGA");
    CHECK(unpack_barcode(rotate_barcode(k, 8, 3), 8) == "TACGAACG");
    CHECK(rotate_barcode(k, 8, 0) == k);
    CHECK(rotate_barcode(k, 8, 8) == k);
    CHECK(rotate_barcode(rotate_barcode(k, 8, 3), 8, 5) == k);
    // A full-width barcode has no unused tail to shift into, which is the case that overflows a
    // mask written for the general one.
    std::string full;
    for (int i = 0; i < 32; ++i) full += "ACGT"[i % 4];
    const uint64_t f = pack_barcode(full);
    CHECK(unpack_barcode(rotate_barcode(f, 32, 1), 32) == full.substr(1) + full.substr(0, 1));
}

TEST_CASE("a spilled counter corrects to the same answer as a resident one") {
    // The claim roadmap items 1 and 2 are one item for: a range partition splits a barcode from
    // its neighbour whenever the substitution lands in the partitioned prefix, so a bucketed
    // correction that took one pass would silently stop correcting exactly those children. Here
    // every position carries children, including the two the 8-bit partition cuts on.
    std::mt19937_64 rng(99);
    const int L = 12;
    std::vector<std::pair<uint64_t, uint32_t>> adds;
    for (int i = 0; i < 4000; ++i) {
        std::string s;
        for (int j = 0; j < L; ++j) s += "ACGT"[rng() % 4];
        const uint64_t parent = umi(s);
        adds.emplace_back(parent, 30);
        // One child per parent, at a position that cycles over the whole barcode: a quarter of
        // them land in the partitioned prefix and are invisible to the first pass.
        const size_t j = static_cast<size_t>(i % L);
        std::string child = s;
        const std::string others = std::string("ACGT").erase(
            std::string("ACGT").find(s[j]), 1);
        child[j] = others[rng() % 3];
        adds.emplace_back(umi(child), 1);
    }

    UmiCounts resident(L);
    UmiCounts spilling(L);
    const std::string dir =
        (std::filesystem::temp_directory_path() / "migec_spill_correct").string();
    std::filesystem::remove_all(dir);
    spilling.enable_spill(dir, 1 << 14, 8);  // 8 bits = a 4-base prefix, spilled many times over
    for (const auto& a : adds) {
        resident.add(a.first, a.second);
        spilling.add(a.first, a.second);
    }
    REQUIRE(spilling.spilled());
    REQUIRE(spilling.distinct() == resident.distinct());

    CorrectionParams p;
    p.threads = 1;
    const CorrectionResult r = correct_umis(resident, p);
    const CorrectionResult s = correct_umis(spilling, p);

    REQUIRE(r.merged > 3000);  // the resident answer is the thing being reproduced
    CHECK(s.merged == r.merged);
    CHECK(s.merged_reads == r.merged_reads);
    CHECK(s.molecules_observed == r.molecules_observed);
    CHECK(s.molecules_corrected == doctest::Approx(r.molecules_corrected).epsilon(0.01));
    CHECK(s.estimated_error == doctest::Approx(r.estimated_error).epsilon(0.05));
    // Never: the scalars are the answer, and the per-entry arrays are absent rather than wrong --
    // they are indexed against entries(), which is the array being bounded.
    CHECK(s.root.empty());
    CHECK(s.corrected.empty());
    std::filesystem::remove_all(dir);
}

TEST_CASE("a bucket left by a dead run is not read back as part of the next one") {
    // Never: a spill appends. A stale bucket is a well-formed bucket, so its counts would be added
    // to this library and nothing would say so.
    const std::string dir =
        (std::filesystem::temp_directory_path() / "migec_spill_stale").string();
    std::filesystem::remove_all(dir);
    std::filesystem::create_directories(dir);
    {
        std::ofstream junk(dir + "/umi_0.bin", std::ios::binary);
        const UmiCounts::Entry ghost{umi("AAAAAAAA"), 424242};
        junk.write(reinterpret_cast<const char*>(&ghost), sizeof(ghost));
    }
    UmiCounts c(8);
    c.enable_spill(dir, 4096, 4);
    std::mt19937_64 rng(11);
    for (int i = 0; i < 20000; ++i) {
        std::string s;
        for (int j = 0; j < 8; ++j) s += "ACGT"[rng() % 4];
        c.add(umi(s), 1);
    }
    REQUIRE(c.spilled());
    uint64_t reads = 0;
    c.for_each([&reads](const UmiCounts::Entry& e) { reads += e.count; });
    CHECK(reads == 20000);  // 424242 would be hard to miss
    std::filesystem::remove_all(dir);
}

TEST_CASE("bucketed correction refuses evidence it cannot index") {
    const std::string dir =
        (std::filesystem::temp_directory_path() / "migec_spill_evidence").string();
    std::filesystem::remove_all(dir);
    UmiCounts c(8);
    c.enable_spill(dir, 64, 4);
    std::mt19937_64 rng(5);
    for (int i = 0; i < 20000; ++i) {
        std::string s;
        for (int j = 0; j < 8; ++j) s += "ACGT"[rng() % 4];
        c.add(umi(s), 1);
    }
    REQUIRE(c.distinct() > 0);
    REQUIRE(c.spilled());
    BarcodeEvidence ev;
    ev.payload_width = 4;
    ev.payload.assign(4 * c.distinct(), 'A');
    CHECK_THROWS_AS(correct_umis(c, CorrectionParams{}, ev), MigecError);
    std::filesystem::remove_all(dir);
}

TEST_CASE("a truncated spill file is refused, not read past") {
    // Never: a spill file is a whole number of entries. A run killed mid-write leaves a partial
    // one, and sizing the buffer at `bytes / 16` while reading `bytes` overruns the heap by up to
    // 15 bytes -- silently, because the read itself succeeds.
    const std::string dir =
        (std::filesystem::temp_directory_path() / "migec_spill_truncated").string();
    std::filesystem::remove_all(dir);
    UmiCounts c(8);
    c.enable_spill(dir, 4096, 4);
    std::mt19937_64 rng(3);
    for (int i = 0; i < 20000; ++i) {
        std::string s;
        for (int j = 0; j < 8; ++j) s += "ACGT"[rng() % 4];
        c.add(umi(s), 1);
    }
    REQUIRE(c.spilled());
    // Chop three bytes off the first non-empty bucket, as a killed process would.
    for (const auto& e : std::filesystem::directory_iterator(dir)) {
        if (std::filesystem::file_size(e.path()) >= 16) {
            std::filesystem::resize_file(e.path(), std::filesystem::file_size(e.path()) - 3);
            break;
        }
    }
    CHECK_THROWS_AS(c.distinct(), MigecError);
    std::filesystem::remove_all(dir);
}

TEST_CASE("a sample id that would escape the output directory is refused") {
    // Never: an id becomes a FILE NAME in every stage, and in assemble and refine it comes from the
    // DATA -- the BC:Z: tag of the first read. `../escape` walks out of the output directory.
    CHECK_THROWS_AS(validate_sample_id("../escape", "test"), MigecError);
    CHECK_THROWS_AS(validate_sample_id("a/b", "test"), MigecError);
    CHECK_THROWS_AS(validate_sample_id("-rf", "test"), MigecError);
    CHECK_THROWS_AS(validate_sample_id("", "test"), MigecError);
    CHECK_THROWS_AS(validate_sample_id(std::string("bad\x01name"), "test"), MigecError);
    validate_sample_id("S1", "test");          // ordinary
    validate_sample_id("A2-i129", "test");     // a donor id with a dash inside
    validate_sample_id("S1.rep2", "test");     // a period is fine; the bucket suffix is 3 digits
}
