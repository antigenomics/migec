#include "doctest.h"

#include <string>

#include "migec/checkout.hpp"
#include "migec/types.hpp"

using namespace migec;

namespace {

constexpr const char* kS1 = "aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN";
constexpr const char* kS2 = "aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN";
const std::string kPayload = "GACTCAGGGTTTCCAGGCCACAACTGCA";

std::string read_for(const char* tag_upper, const std::string& umi_a, const std::string& umi_b,
                     const std::string& umi_c) {
    return std::string("aa") + tag_upper + "cagtggtatcaacgcagagt" + umi_a + "t" + umi_b + "t" +
           umi_c + kPayload;
}

}  // namespace

TEST_CASE("checkout assigns, extracts the UMI and trims the tag away") {
    PatternSet set;
    set.add("S1", kS1);
    set.add("S2", kS2);
    Checkout co(set, CheckoutParams{});

    const std::string read = read_for("ACT", "AAAA", "CCCC", "GGGG");
    const std::string qual(read.size(), 'I');

    CheckoutRead r = co.process(read, qual);
    REQUIRE(r.ok);
    CHECK(r.sample == 0);
    CHECK(r.umi == "AAAACCCCGGGG");
    // Trimming leaves exactly the biological payload -- no adapter, no sample tag, no UMI.
    CHECK(std::string(r.seq) == kPayload);
    CHECK(r.qual.size() == r.seq.size());
    CHECK(co.counters().assigned == 1);
    CHECK(co.counters().per_sample[0] == 1);
}

TEST_CASE("trim modes") {
    PatternSet set;
    set.add("S1", kS1);
    const std::string read = read_for("ACT", "AAAA", "CCCC", "GGGG");
    const std::string qual(read.size(), 'I');

    CheckoutParams none;
    none.trim = TrimMode::kNone;
    Checkout keep(set, none);
    CheckoutRead r = keep.process(read, qual);
    REQUIRE(r.ok);
    CHECK(std::string(r.seq) == read);  // untouched; the UMI is still in the header only
}

TEST_CASE("the header tag string is SAM-conformant") {
    const std::string tags = Checkout::header_tags("AAAACCCCGGGG", "IIIIIIIIIIII", "S1");
    CHECK(tags == "RX:Z:AAAACCCCGGGG\tQX:Z:IIIIIIIIIIII\tBC:Z:S1");
    // TABs, not spaces: bwa -C copies the comment verbatim into the SAM record.
    CHECK(tags.find(' ') == std::string::npos);

    CHECK(Checkout::header_tags("ACGT", "", "") == "RX:Z:ACGT");
    CHECK(Checkout::header_tags("", "", "S3") == "BC:Z:S3");
    CHECK(Checkout::header_tags("", "", "").empty());
}

TEST_CASE("an unmatched read is counted, not assigned") {
    PatternSet set;
    set.add("S1", kS1);
    Checkout co(set, CheckoutParams{});
    const std::string read = "GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA";
    CheckoutRead r = co.process(read, std::string(read.size(), 'I'));
    CHECK_FALSE(r.ok);
    CHECK(co.counters().unmatched == 1);
    CHECK(co.counters().assigned == 0);
    CHECK(co.counters().total == 1);
}

TEST_CASE("an ambiguous sample tag is counted separately from an unmatched read") {
    // The diagnostic distinction that matters: "your barcodes are too close together" is a
    // different problem from "your pattern is wrong", and one number cannot say both.
    PatternSet set;
    set.add("S1", kS1);
    set.add("S2", kS2);
    Checkout co(set, CheckoutParams{});
    const std::string read = read_for("AGT", "AAAA", "CCCC", "GGGG");  // 1 sub from each tag
    CheckoutRead r = co.process(read, std::string(read.size(), 'I'));
    CHECK_FALSE(r.ok);
    CHECK(co.counters().ambiguous == 1);
    CHECK(co.counters().unmatched == 0);
}

TEST_CASE("a low-quality UMI is kept by default and dropped only on request") {
    PatternSet set;
    set.add("S1", kS1);
    const std::string read = read_for("ACT", "AAAA", "CCCC", "GGGG");
    std::string qual(read.size(), 'I');
    // Wreck one UMI base. Its position: 2 + 3 + 20 = 25 is the first UMI base.
    qual[25] = '#';  // Q2

    // Default: the read survives. Discarding a molecule because one UMI base read badly throws
    // away sequence the correction step can usually recover.
    Checkout lenient(set, CheckoutParams{});
    CHECK(lenient.process(read, qual).ok);

    CheckoutParams strict;
    strict.min_umi_quality = 15;  // MIGEC's old default
    Checkout picky(set, strict);
    CheckoutRead r = picky.process(read, qual);
    CHECK_FALSE(r.ok);
    CHECK(picky.counters().bad_umi == 1);
}

TEST_CASE("a read that is all tag and no payload is dropped") {
    PatternSet set;
    set.add("S1", kS1);
    CheckoutParams p;
    p.min_payload = 10;
    Checkout co(set, p);
    const std::string read = std::string("aaACTcagtggtatcaacgcagagtAAAAtCCCCtGGGG") + "ACGT";
    CheckoutRead r = co.process(read, std::string(read.size(), 'I'));
    CHECK_FALSE(r.ok);
    CHECK(co.counters().short_payload == 1);
}

TEST_CASE("counters account for every read exactly once") {
    PatternSet set;
    set.add("S1", kS1);
    set.add("S2", kS2);
    Checkout co(set, CheckoutParams{});

    co.process(read_for("ACT", "AAAA", "CCCC", "GGGG"), std::string(67, 'I'));
    co.process(read_for("AGA", "TTTT", "CCCC", "GGGG"), std::string(67, 'I'));
    co.process(read_for("AGT", "AAAA", "CCCC", "GGGG"), std::string(67, 'I'));
    const std::string junk(60, 'G');
    co.process(junk, std::string(junk.size(), 'I'));

    const CheckoutCounters& c = co.counters();
    CHECK(c.total == 4);
    CHECK(c.assigned + c.unmatched + c.ambiguous + c.short_payload + c.bad_umi == c.total);
    CHECK(c.assigned == 2);
}

TEST_CASE("quality calibration fits the slope and drops variable positions") {
    QualityCalibration c;
    // Two quality levels, each with an error rate exactly twice the nominal, over four positions.
    for (size_t p = 0; p < 4; ++p) {
        c.by_position[p][30][0] = 1'000'000 - 2000;   // Q30 nominal 1e-3, observed 2e-3
        c.by_position[p][30][1] = 2000;
        c.by_position[p][20][0] = 1'000'000 - 20000;  // Q20 nominal 1e-2, observed 2e-2
        c.by_position[p][20][1] = 20000;
    }
    // ...and one position that is simply variable: 30% mismatch at every quality.
    c.by_position[4][30][0] = 700'000;
    c.by_position[4][30][1] = 300'000;

    c.fit();
    CHECK(c.fitted);
    CHECK(c.positions_dropped == 1);
    CHECK(c.slope == doctest::Approx(2.0).epsilon(0.02));
    CHECK(c.quality_independent < 1e-4);
    // error() applies the slope and NOT the intercept: the intercept belongs to the standard
    // being measured against, which is a synthesised oligo.
    CHECK(c.error(30) == doctest::Approx(2e-3).epsilon(0.02));
}

TEST_CASE("a calibration with one quality level declines to fit") {
    QualityCalibration c;
    c.by_position[0][38][0] = 10'000'000;
    c.by_position[0][38][1] = 1000;
    c.fit();
    CHECK_FALSE(c.fitted);
    // ...and falls back to the nominal rate rather than to a line through one point.
    CHECK(c.error(38) == doctest::Approx(phred_error(38)));
}

TEST_CASE("an intercept the fit puts below zero is reported as zero, not as a negative rate") {
    QualityCalibration c;
    c.by_position[0][30][0] = 1'000'000; c.by_position[0][30][1] = 100;      // 1e-4, below nominal
    c.by_position[0][10][0] = 1'000'000; c.by_position[0][10][1] = 200'000;  // 2e-1
    c.fit();
    CHECK(c.fitted);
    CHECK(c.quality_independent >= 0.0);
    CHECK(c.slope >= 0.0);
}
