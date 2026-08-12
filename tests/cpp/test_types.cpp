#include "doctest.h"

#include <string>

#include "migec/types.hpp"

using namespace migec;

TEST_CASE("barcode packing preserves order") {
    // The range partition and the sort both assume packed order == lexicographic order.
    CHECK(pack_barcode("AAAA") < pack_barcode("AAAC"));
    CHECK(pack_barcode("AAAC") < pack_barcode("AAAG"));
    CHECK(pack_barcode("AAAG") < pack_barcode("AAAT"));
    CHECK(pack_barcode("ACGT") < pack_barcode("CAAA"));
    CHECK(pack_barcode("TTTT") > pack_barcode("TTTG"));
}

TEST_CASE("barcode packing round trips and flags N") {
    bool has_n = false;
    uint64_t p = pack_barcode("ACGTACGTACGT", &has_n);
    CHECK_FALSE(has_n);
    CHECK(unpack_barcode(p, 12) == "ACGTACGTACGT");

    p = pack_barcode("ACGNACGT", &has_n);
    CHECK(has_n);
    CHECK(unpack_barcode(p, 8) == "ACGAACGT");  // N is stored as A, ambiguity is in the flag

    CHECK(unpack_barcode(pack_barcode(""), 0).empty());
    CHECK_THROWS_AS(pack_barcode(std::string(33, 'A')), MigecError);
}

TEST_CASE("bucket_of uses the top bits and is monotone") {
    CHECK(bucket_of(pack_barcode("AAAAAAAAAAAA"), 4) == 0);
    CHECK(bucket_of(pack_barcode("TTTTTTTTTTTT"), 4) == 15);
    CHECK(bucket_of(pack_barcode("ACGT"), 0) == 0);  // 0 bits == one bucket
    // Monotone: a lexicographically larger barcode never lands in an earlier bucket.
    uint32_t prev = 0;
    for (char c : std::string("ACGT")) {
        uint32_t b = bucket_of(pack_barcode(std::string(1, c) + "AAAAAAAAAAA"), 2);
        CHECK(b >= prev);
        prev = b;
    }
}

TEST_CASE("reverse_complement reverses quality too") {
    std::string seq = "ACGTN", qual = "12345";
    reverse_complement(seq, qual);
    CHECK(seq == "NACGT");
    CHECK(qual == "54321");

    std::string odd = "ACG", oq = "abc";
    reverse_complement(odd, oq);
    CHECK(odd == "CGT");
    CHECK(oq == "cba");

    std::string s2 = "ACGT", q2 = "AB";
    CHECK_THROWS_AS(reverse_complement(s2, q2), MigecError);
}

TEST_CASE("phred conversion is sanger and clamped") {
    CHECK(phred_from_char('!') == 0);
    CHECK(phred_from_char('I') == 40);
    CHECK(char_from_phred(40) == 'I');
    CHECK(phred_from_char(' ') == 0);          // below the offset clamps rather than wraps
    CHECK(phred_from_char('~') == kMaxPhred);  // above the cap clamps
    CHECK(phred_error(0) == doctest::Approx(1.0));
    CHECK(phred_error(30) == doctest::Approx(1e-3));
}

TEST_CASE("iupac masks match the standard") {
    CHECK(iupac_mask('A') == 0b0001);
    CHECK(iupac_mask('n') == 0b1111);
    CHECK(iupac_size(iupac_mask('N')) == 4);
    CHECK(iupac_size(iupac_mask('R')) == 2);
    CHECK((iupac_mask('R') >> base_code('A') & 1) == 1);  // R = A or G
    CHECK((iupac_mask('R') >> base_code('C') & 1) == 0);
    CHECK(iupac_mask('Z') == 0);  // not a IUPAC symbol
}
