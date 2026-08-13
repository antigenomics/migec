#include "doctest.h"

#include <string>

#include "migec/pattern.hpp"
#include "migec/types.hpp"

using namespace migec;

namespace {

// The real MIGEC barcode table from misc/barcodes.txt: two fuzzy bases, a 3 nt degenerate sample
// tag, the SMART adapter in lowercase, and a 12 nt UMI split by lowercase spacers.
constexpr const char* kS1 = "aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN";
constexpr const char* kS2 = "aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN";

std::string qual_of(const std::string& s, char q = 'I') { return std::string(s.size(), q); }

}  // namespace

TEST_CASE("pattern compiles the MIGEC dialect") {
    BarcodePattern p = BarcodePattern::compile(kS1);
    CHECK(p.size() == 39);
    CHECK(p.umi_length() == 12);       // NNNN t NNNN t NNNN
    CHECK(p.scored_positions() == 27);  // everything but the 12 UMI positions

    CHECK_THROWS_AS(BarcodePattern::compile(""), MigecError);
    CHECK_THROWS_AS(BarcodePattern::compile("ACGTZ"), MigecError);
    CHECK_THROWS_AS(BarcodePattern::compile("NNNN"), MigecError);  // nothing to score
}

TEST_CASE("pattern extracts a non-contiguous UMI at the right offset") {
    BarcodePattern p = BarcodePattern::compile(kS1);
    //                      aa ACT cagtggtatcaacgcagagt NNNN t NNNN t NNNN
    const std::string tag = "aaACTcagtggtatcaacgcagagtAAAAtCCCCtGGGG";
    const std::string read = "GGGCC" + tag + "ACGTACGTACGTACGT";
    std::string q = qual_of(read);

    PatternMatch m = p.match(read, q);
    REQUIRE(m.found);
    CHECK(m.offset == 5);
    CHECK(m.umi == "AAAACCCCGGGG");  // the three runs concatenate into one 12 nt UMI
    CHECK(m.umi_qual == std::string(12, 'I'));
    CHECK(m.payload_begin == 5 + 39);
    // Trimming leaves exactly the payload.
    CHECK(read.substr(static_cast<size_t>(m.payload_begin)) == "ACGTACGTACGTACGT");
}

TEST_CASE("a mismatch on a bad base is cheap, on a good base it is fatal") {
    BarcodePattern p = BarcodePattern::compile("ACGTACGTACGTACGT");
    const std::string good = "ACGTACGTACGTACGT";
    std::string bad = good;
    bad[0] = 'T';  // one substitution in an uppercase position

    // High quality everywhere: the mismatch costs ~9.55 bits against a ~2 bits/base match.
    PatternMatch hi = p.match(bad, qual_of(bad, 'I'));
    // Low quality at the mismatching base: nearly free.
    std::string q = qual_of(bad, 'I');
    q[0] = '#';  // Q2
    PatternMatch lo = p.match(bad, q);

    REQUIRE(hi.found);
    REQUIRE(lo.found);
    CHECK(lo.score > hi.score);
    CHECK(hi.score < p.match(good, qual_of(good)).score);
}

TEST_CASE("quality is read at the match offset, not the read start") {
    // MIGEC v1 indexed quality from position 0 regardless of where the tag matched. Put the tag
    // late in the read, make the *leading* bases terrible and the tag's own bases perfect: a
    // correctly indexed scorer is unaffected.
    BarcodePattern p = BarcodePattern::compile("ACGTACGTACGT");
    const std::string read = "TTTTTTTTTTTTTTTTTTTTACGTACGTACGT";
    std::string q(read.size(), 'I');
    for (int i = 0; i < 20; ++i) q[static_cast<size_t>(i)] = '#';  // Q2 before the tag

    PatternMatch m = p.match(read, q);
    REQUIRE(m.found);
    CHECK(m.offset == 20);
    // Same score as if the read were only the tag at full quality.
    PatternMatch bare = p.match("ACGTACGTACGT", std::string(12, 'I'));
    CHECK(m.score == doctest::Approx(bare.score).epsilon(1e-9));
}

TEST_CASE("degenerate IUPAC positions match their whole set") {
    BarcodePattern p = BarcodePattern::compile("ARYTACGTACGT");  // R = A|G, Y = C|T
    for (const std::string s : {"AACTACGTACGT", "AGCTACGTACGT", "AATTACGTACGT", "AGTTACGTACGT"}) {
        PatternMatch m = p.match(s, qual_of(s));
        INFO(s);
        CHECK(m.found);
    }
    // A base outside the set is a mismatch, and a degenerate match is worth less than an exact one
    // -- two possible bases is one bit less information.
    PatternMatch exact = p.match("AACTACGTACGT", qual_of("AACTACGTACGT"));
    BarcodePattern q = BarcodePattern::compile("AACTACGTACGT");
    PatternMatch full = q.match("AACTACGTACGT", qual_of("AACTACGTACGT"));
    CHECK(full.score > exact.score);
}

TEST_CASE("lowercase is a half-weight region, not a wildcard") {
    // MIGEC treated lowercase as matching anything, throwing the evidence away. Here a lowercase
    // mismatch still costs, just half as much.
    BarcodePattern upper = BarcodePattern::compile("ACGTACGTACGT");
    BarcodePattern lower = BarcodePattern::compile("acgtACGTACGT");
    std::string bad = "TCGTACGTACGT";  // mismatch at position 0

    const double du = upper.match("ACGTACGTACGT", qual_of(bad)).score -
                      upper.match(bad, qual_of(bad)).score;
    const double dl = lower.match("ACGTACGTACGT", qual_of(bad)).score -
                      lower.match(bad, qual_of(bad)).score;
    CHECK(dl < du);          // lowercase penalises less
    CHECK(dl > 0.0);         // but it still penalises
}

TEST_CASE("a read with no tag is rejected") {
    BarcodePattern p = BarcodePattern::compile(kS1);
    const std::string read(120, 'A');
    CHECK_FALSE(p.match(read, qual_of(read)).found);

    const std::string randomish =
        "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
    CHECK_FALSE(p.match(randomish, qual_of(randomish)).found);
}

TEST_CASE("max_offset anchors positional chemistries") {
    // 10x style: cell barcode then UMI at a fixed position. Anchoring means a chance match later
    // in the read cannot steal the placement.
    BarcodePattern p = BarcodePattern::compile("ACGTACGTACGTNNNNNNNNNNNN");
    const std::string read = "ACGTACGTACGTAAAACCCCGGGGTTTTTTTTACGTACGTACGTTTTTTTTTTTTT";
    MatchParams anchored;
    anchored.max_offset = 0;
    PatternMatch m = p.match(read, qual_of(read), anchored);
    REQUIRE(m.found);
    CHECK(m.offset == 0);
    CHECK(m.umi == "AAAACCCCGGGG");
}

TEST_CASE("pattern set assigns the right sample and flags ambiguity") {
    PatternSet set;
    set.add("S1", kS1);
    set.add("S2", kS2);

    const std::string r1 = "aaACTcagtggtatcaacgcagagtAAAAtCCCCtGGGGACGTACGT";
    const std::string r2 = "aaAGAcagtggtatcaacgcagagtTTTTtCCCCtGGGGACGTACGT";

    auto a1 = set.assign(r1, qual_of(r1));
    REQUIRE(a1.sample == 0);
    CHECK(set.samples()[static_cast<size_t>(a1.sample)] == "S1");
    CHECK(a1.match.umi == "AAAACCCCGGGG");

    auto a2 = set.assign(r2, qual_of(r2));
    REQUIRE(a2.sample == 1);
    CHECK(a2.match.umi == "TTTTCCCCGGGG");

    // A tag that is equidistant from both samples must not be assigned to either. S1 is ACT and
    // S2 is AGA; "AGT" is one substitution from each.
    const std::string amb = "aaAGTcagtggtatcaacgcagagtAAAAtCCCCtGGGGACGTACGT";
    auto a3 = set.assign(amb, qual_of(amb));
    CHECK(a3.sample == -1);
    CHECK(a3.ambiguous);
}

TEST_CASE("an unrelated read is assigned to no sample and is not ambiguous") {
    PatternSet set;
    set.add("S1", kS1);
    set.add("S2", kS2);
    const std::string read(80, 'G');
    auto a = set.assign(read, qual_of(read));
    CHECK(a.sample == -1);
    CHECK_FALSE(a.ambiguous);  // no match at all is a different outcome from a tie
}

TEST_CASE("the quality calibration table overrides nominal phred") {
    // Two-colour instruments emit few distinct Q values and the nominal error is optimistic.
    // Feeding a pessimistic table must make a high-Q mismatch cheaper.
    BarcodePattern p = BarcodePattern::compile("ACGTACGTACGT");
    const std::string bad = "TCGTACGTACGT";
    std::string q(bad.size(), 'I');  // Q40 nominal -> 1e-4

    MatchParams nominal;
    MatchParams calibrated;
    calibrated.quality_calibration.assign(61, 1e-2);  // measured: really 1%

    const double s_nom = p.match(bad, q, nominal).score;
    const double s_cal = p.match(bad, q, calibrated).score;
    CHECK(s_cal > s_nom);
}

TEST_CASE("a near-tie placement is ambiguous, not silently resolved by the prune") {
    // The offset scan abandons an offset once it cannot reach the bar. Setting that bar at the
    // incumbent best drops any runner-up that lands *within* min_margin of it -- exactly the
    // offsets the margin exists to detect. The margin then comes back as best-minus-nothing and
    // the read is reported as an unambiguous match at whichever placement came first.
    const BarcodePattern p = BarcodePattern::compile("ACGTACGTAC");
    MatchParams mp;
    mp.min_margin = 5.0;

    // The tag, then the tag again with one mismatch on a Q2 base: worth -0.60 bits instead of
    // +2.00, so the placements are 2.6 bits apart -- inside the 5-bit margin.
    std::string qual(20, 'I');
    qual[15] = '#';  // Q2
    const PatternMatch tie = p.match("ACGTACGTAC" "ACGTAAGTAC", qual, mp);
    CHECK_FALSE(tie.found);

    // A single placement still matches, and still reports the full score as its margin.
    const PatternMatch clear = p.match("ACGTACGTAC" "TTTTTTTTTT", std::string(20, 'I'), mp);
    CHECK(clear.found);
    CHECK(clear.offset == 0);
    CHECK(clear.margin > 5.0);
}

TEST_CASE("the acceptance bar is charged for the offsets actually scanned") {
    // MAGERI's dual-end handle: five bases, 10 bits. Enough against one offset and not against
    // sixty -- so billing an anchored scan for a scan it never performs refuses every read of a
    // design that is perfectly well determined.
    const BarcodePattern p = BarcodePattern::compile("NNNNNNNNNNNNTGACT");
    const std::string seq = "ACGTACGTACGTTGACT" + std::string(60, 'A');
    const std::string qual(seq.size(), char_from_phred(35));

    MatchParams free_scan;
    CHECK_FALSE(p.match(seq, qual, free_scan).found);

    MatchParams anchored;
    anchored.max_offset = 0;
    const PatternMatch m = p.match(seq, qual, anchored);
    CHECK(m.found);
    CHECK(m.umi == "ACGTACGTACGT");
    CHECK(p.default_min_score(seq.size(), 1, 0.01, 0) <
          p.default_min_score(seq.size(), 1, 0.01, -1));
}

TEST_CASE("a slave pattern extends the UMI rather than starting a new one") {
    PatternSet set;
    set.add("S1", "NNNNNNNNNNNNTGACT", "AGTCANNNNNNNNNNNN");
    CHECK(set.has_slave(0));
    CHECK(set.pattern(0).umi_length() == 12);
    CHECK(set.slave(0).umi_length() == 12);
    CHECK(set.umi_length(0) == 24);

    PatternSet plain;
    plain.add("S2", "NNNNNNNNNNNNTGACT");
    CHECK_FALSE(plain.has_slave(0));
    CHECK(plain.umi_length(0) == 12);
}
