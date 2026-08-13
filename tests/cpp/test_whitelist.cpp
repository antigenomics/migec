#include "doctest.h"

#include <cstdio>
#include <fstream>
#include <random>
#include <string>
#include <vector>

#include "migec/types.hpp"
#include "migec/whitelist.hpp"

using namespace migec;

namespace {

std::string write_list(const std::vector<std::string>& entries, const char* name) {
    std::string path = std::string(std::tmpnam(nullptr)) + name;
    std::ofstream fh(path);
    for (const std::string& e : entries) fh << e << "\n";
    return path;
}

std::vector<std::string> random_list(size_t n, int len, uint32_t seed) {
    std::mt19937 rng(seed);
    std::vector<std::string> out;
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        std::string s(static_cast<size_t>(len), 'A');
        for (char& c : s) c = "ACGT"[rng() & 3u];
        out.push_back(s);
    }
    return out;
}

}  // namespace

TEST_CASE("a whitelist of two lengths is refused rather than truncated") {
    const std::string path = write_list({"ACGTACGT", "ACGTAC"}, "_mixed");
    CHECK_THROWS_AS(Whitelist::load(path), MigecError);
    std::remove(path.c_str());
}

TEST_CASE("the 10x gem-group suffix is tolerated") {
    const std::string path = write_list({"ACGTACGTACGTACGT-1", "TTTTGGGGCCCCAAAA-1"}, "_suffix");
    const Whitelist w = Whitelist::load(path);
    CHECK(w.size() == 2);
    CHECK(w.length() == 16);
    CHECK(w.contains(pack_barcode("ACGTACGTACGTACGT")));
    std::remove(path.c_str());
}

TEST_CASE("a barcode one substitution off a whitelist entry is snapped to it") {
    std::vector<std::string> entries = random_list(5000, 16, 3);
    entries.push_back("ACGTACGTACGTACGT");
    const std::string path = write_list(entries, "_snap");
    const Whitelist w = Whitelist::load(path);
    std::vector<uint32_t> counts(w.size(), 10);

    WhitelistParams p;
    // Q20 at the base that differs: a plausible miscall.
    std::string qual(16, char_from_phred(30));
    qual[15] = char_from_phred(20);
    CHECK(w.correct("ACGTACGTACGTACGA", qual, counts, p, 1e-8) == "ACGTACGTACGTACGT");
    std::remove(path.c_str());
}

TEST_CASE("the background hypothesis stops everything being assigned to something") {
    // The failure this file exists to prevent. A barcode two substitutions from every entry has
    // no candidate at all; one that is one away must still lose when the background is likely and
    // the base that differs was called confidently.
    std::vector<std::string> entries = random_list(5000, 16, 4);
    entries.push_back("ACGTACGTACGTACGT");
    const std::string path = write_list(entries, "_bg");
    const Whitelist w = Whitelist::load(path);
    std::vector<uint32_t> counts(w.size(), 10);
    WhitelistParams p;

    const std::string confident(16, char_from_phred(40));
    // The priors are per barcode. Plenty of off-list barcodes about: the background wins even at
    // distance 1, which is the whole point -- an index-hopped or undeclared barcode must not be
    // absorbed into whichever entry happens to be nearest.
    CHECK(w.correct("ACGTACGTACGTACGA", confident, counts, p, 1e-4).empty());
    // Almost none: the same observation is now best explained as a miscall.
    CHECK(w.correct("ACGTACGTACGTACGA", confident, counts, p, 1e-12) == "ACGTACGTACGTACGT");
    std::remove(path.c_str());
}

TEST_CASE("an N is expanded, not discarded") {
    std::vector<std::string> entries = random_list(5000, 16, 5);
    entries.push_back("ACGTACGTACGTACGT");
    const std::string path = write_list(entries, "_n");
    const Whitelist w = Whitelist::load(path);
    std::vector<uint32_t> counts(w.size(), 10);
    WhitelistParams p;
    const std::string qual(16, char_from_phred(30));
    // The instrument declined to call the last base. It is consistent with all four, so the
    // whitelist entry is recoverable -- discarding the barcode would lose a real molecule.
    CHECK(w.correct("ACGTACGTACGTACGN", qual, counts, p, 1e-8) == "ACGTACGTACGTACGT");
    std::remove(path.c_str());
}

TEST_CASE("a used whitelist entry is a likelier parent than an unused one") {
    std::vector<std::string> entries = random_list(2000, 16, 6);
    entries.push_back("ACGTACGTACGTACGT");   // will be given reads
    entries.push_back("ACGTACGTACGTACGA");   // observed exactly; the neighbour of the above
    const std::string path = write_list(entries, "_prior");
    const Whitelist w = Whitelist::load(path);

    WhitelistParams p;
    const std::string qual(16, char_from_phred(25));
    std::vector<uint32_t> counts(w.size(), 0);
    counts[w.index_of(pack_barcode("ACGTACGTACGTACGT"))] = 100000;
    // Observing the neighbour of a heavily used entry, with a mediocre base where they differ.
    const std::string heavy = w.correct("ACGTACGTACGTACGA", qual, counts, p, 1e-7);

    std::vector<uint32_t> flat(w.size(), 0);
    const std::string unused = w.correct("ACGTACGTACGTACGA", qual, flat, p, 1e-7);
    CHECK(heavy == "ACGTACGTACGTACGT");
    CHECK(unused.empty());
    std::remove(path.c_str());
}


TEST_CASE("the background prior is measured, not assumed") {
    // Half the reads on barcodes that cannot be errors of anything on the list, spread over a
    // million distinct off-list barcodes, is 5e-7 for any one of them -- not 0.5.
    CHECK(Whitelist::measure_background(500, 1000, 1'000'000) == doctest::Approx(5e-7));
    CHECK(Whitelist::measure_background(0, 1000, 10) == doctest::Approx(0.0));
    CHECK(Whitelist::measure_background(10, 0, 10) == doctest::Approx(0.0));
}
