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
