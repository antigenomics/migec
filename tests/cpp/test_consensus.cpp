#include "doctest.h"

#include <cmath>
#include <deque>
#include <random>
#include <string>
#include <vector>

#include "migec/consensus.hpp"
#include "migec/types.hpp"

using namespace migec;

namespace {

// The reads are views, so the strings have to outlive the group. A deque never moves what it
// already holds, and nothing here is big enough for the leak to matter.
std::vector<ConsensusRead> group(std::vector<std::string> seqs, char q = 'I') {
    static std::deque<std::string> arena;
    std::vector<ConsensusRead> out;
    for (const std::string& s : seqs) {
        arena.push_back(s);
        const std::string& seq = arena.back();
        arena.push_back(std::string(s.size(), q));
        out.push_back({seq, arena.back()});
    }
    return out;
}

}  // namespace

TEST_CASE("a clean group returns its own sequence") {
    ConsensusParams p;
    auto reads = group({"ACGTACGT", "ACGTACGT", "ACGTACGT"});
    auto out = assemble_group(reads, p);
    REQUIRE(out.size() == 1);
    CHECK(out[0].seq == "ACGTACGT");
    CHECK(out[0].reads == 3);
}

TEST_CASE("one erroneous read is outvoted") {
    ConsensusParams p;
    auto reads = group({"ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT"});
    auto out = assemble_group(reads, p);
    REQUIRE(out.size() == 1);
    CHECK(out[0].seq == "ACGTACGT");
}

TEST_CASE("emitted quality never exceeds the RT floor") {
    ConsensusParams p;
    p.rt_floor = 1e-4;  // -> Q40 cap
    std::vector<std::string> seqs(50, "ACGTACGTACGT");
    auto out = assemble_group(group(seqs, char_from_phred(40)), p);
    REQUIRE(out.size() == 1);
    for (char c : out[0].qual) CHECK(phred_from_char(c) <= 40);

    // An error before amplification is in every read, so more reads must not buy more quality.
    p.rt_floor = 1e-2;  // -> Q20 cap
    auto capped = assemble_group(group(seqs, char_from_phred(40)), p);
    for (char c : capped[0].qual) CHECK(phred_from_char(c) <= 20);
}

TEST_CASE("quality rises with coverage below the floor") {
    ConsensusParams p;
    p.rt_floor = 1e-12;  // effectively uncapped, so the consensus posterior is what is measured
    auto small = assemble_group(group({"AC", "AC"}, char_from_phred(20)), p);
    auto large = assemble_group(group(std::vector<std::string>(12, "AC"), char_from_phred(20)), p);
    CHECK(phred_from_char(large[0].qual[0]) > phred_from_char(small[0].qual[0]));
}

TEST_CASE("a position no read calls is N at quality zero") {
    ConsensusParams p;
    auto out = assemble_group(group({"ANGT", "ANGT", "ANGT"}), p);
    REQUIRE(out.size() == 1);
    CHECK(out[0].seq == "ANGT");
    CHECK(phred_from_char(out[0].qual[1]) == 0);
}

TEST_CASE("two molecules in one MIG are split on linkage, not on site count") {
    // 20 and 20. The threshold implies a minimum group size, because the strongest evidence a
    // pair of columns can carry is log10 C(n, n/2): a 50/50 split needs ~34 reads before it can
    // clear 8.68 at all. That is the price of a measured false-positive rate.
    ConsensusParams p;
    std::vector<std::string> seqs;
    for (int i = 0; i < 20; ++i) seqs.push_back("ACGTACGTACGTACGTACGT");
    for (int i = 0; i < 20; ++i) seqs.push_back("ACGTTCGTACGAACGTACGT");  // differs at 4 and 11
    auto out = assemble_group(group(seqs), p);
    REQUIRE(out.size() == 2);
    CHECK(out[0].seq == "ACGTACGTACGTACGTACGT");
    CHECK(out[1].seq == "ACGTTCGTACGAACGTACGT");
    CHECK(out[0].reads == 20);
    CHECK(out[1].reads == 20);
    CHECK(out[0].linkage > p.linkage_threshold);
}

TEST_CASE("scattered error at the same positions does not split") {
    // Two positions each carry 16 minor bases, as in the split case above, but they overlap on
    // only the 8 reads chance would give. Counting polymorphic sites cannot tell these apart;
    // linkage can.
    ConsensusParams p;
    std::vector<std::string> seqs;
    for (int i = 0; i < 40; ++i) seqs.push_back("ACGTACGTACGTACGTACGT");
    for (int i = 0; i < 16; ++i) seqs[static_cast<size_t>(i)][4] = 'T';
    for (int i = 8; i < 24; ++i) seqs[static_cast<size_t>(i)][11] = 'A';
    auto out = assemble_group(group(seqs), p);
    CHECK(out.size() == 1);
    CHECK(out[0].linkage < p.linkage_threshold);
}

TEST_CASE("a bad read carrying minors everywhere is not a subclone") {
    // Two reads are minor at every callable position. Their co-segregation is perfect, which is
    // exactly what a low-quality read looks like, and X3's null is what says so.
    ConsensusParams p;
    std::vector<std::string> seqs(20, "ACGTACGTACGTACGTACGT");
    for (size_t i = 0; i < 2; ++i) {
        for (size_t j = 0; j < 20; j += 4) seqs[i][j] = 'T';
    }
    auto out = assemble_group(group(seqs), p);
    CHECK(out.size() == 1);
}

TEST_CASE("a small group is never split") {
    ConsensusParams p;
    auto out = assemble_group(group({"ACGTACGT", "TCGTACGA", "ACGTACGT"}), p);
    CHECK(out.size() == 1);
}

TEST_CASE("consensus error is below 1e-5 at coverage 5") {
    // The M1 gate, in miniature: random reads at Q30 over a fixed template, five per molecule.
    ConsensusParams p;
    p.rt_floor = 1e-12;  // the gate is on the consensus, not on the chemistry
    std::mt19937 rng(7);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    const std::string truth = "ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT";
    const double e = 1e-3;
    size_t mismatches = 0, bases = 0;
    for (int m = 0; m < 400; ++m) {
        std::vector<std::string> seqs;
        for (int r = 0; r < 5; ++r) {
            std::string s = truth;
            for (char& c : s) {
                if (unit(rng) < e) {
                    const char* alt = "ACGT";
                    char pick = alt[static_cast<int>(unit(rng) * 4) & 3];
                    if (pick == c) pick = (c == 'A') ? 'C' : 'A';
                    c = pick;
                }
            }
            seqs.push_back(s);
        }
        auto out = assemble_group(group(seqs, char_from_phred(30)), p);
        REQUIRE(out.size() == 1);
        for (size_t j = 0; j < truth.size(); ++j) {
            ++bases;
            mismatches += out[0].seq[j] != truth[j];
        }
    }
    CHECK(static_cast<double>(mismatches) / static_cast<double>(bases) < 1e-5);
}

TEST_CASE("linkage_score is zero when there is nothing to test") {
    ConsensusParams p;
    CHECK(linkage_score(group({"ACGT", "ACGT", "ACGT", "ACGT", "ACGT", "ACGT"}), p) == 0.0);
}

TEST_CASE("tiled reads are placed into one contig") {
    ConsensusParams p;
    p.contig = true;
    std::mt19937 rng(11);
    std::string molecule;
    for (int i = 0; i < 300; ++i) molecule += "ACGT"[rng() & 3];
    std::vector<std::string> seqs;
    for (int start = 0; start + 90 <= 300; start += 30) seqs.push_back(molecule.substr(start, 90));
    auto out = assemble_group(group(seqs), p);
    REQUIRE(out.size() == 1);
    CHECK(out[0].seq == molecule);
    CHECK(out[0].components == 1);
}

TEST_CASE("a gap between two islands is never bridged") {
    // The two islands share a barcode and no sequence. One consensus over them would assert 100
    // nt that no read covers, which is the failure X1 says is mandatory to avoid.
    ConsensusParams p;
    p.contig = true;
    std::mt19937 rng(12);
    std::string molecule;
    for (int i = 0; i < 400; ++i) molecule += "ACGT"[rng() & 3];
    std::vector<std::string> seqs;
    for (int start : {0, 30, 60, 250, 280, 310}) seqs.push_back(molecule.substr(static_cast<size_t>(start), 90));
    auto out = assemble_group(group(seqs), p);
    REQUIRE(out.size() == 2);
    CHECK(out[0].components == 2);
    CHECK(out[0].seq.size() == 150);
    CHECK(out[1].seq.size() == 150);
    CHECK(out[0].seq == molecule.substr(0, 150));
    CHECK(out[1].seq == molecule.substr(250, 150));
}

TEST_CASE("placement is refused when the overlap is too short to be evidence") {
    ConsensusParams p;
    p.contig = true;
    p.min_overlap = 40;
    std::mt19937 rng(13);
    std::string molecule;
    for (int i = 0; i < 200; ++i) molecule += "ACGT"[rng() & 3];
    // 20 nt of overlap, under the 40 required.
    auto out = assemble_group(group({molecule.substr(0, 90), molecule.substr(70, 90)}), p);
    CHECK(out.size() == 2);
}

TEST_CASE("amplicon mode leaves every read at offset zero") {
    ConsensusParams p;  // contig defaults to false
    std::mt19937 rng(14);
    std::string molecule;
    for (int i = 0; i < 300; ++i) molecule += "ACGT"[rng() & 3];
    std::vector<std::string> seqs;
    for (int start = 0; start + 90 <= 300; start += 30) seqs.push_back(molecule.substr(start, 90));
    auto out = assemble_group(group(seqs), p);
    REQUIRE(out.size() == 1);
    CHECK(out[0].seq.size() == 90);
    CHECK(out[0].components == 1);
}

TEST_CASE("placing reads is deterministic whatever order they arrive in") {
    ConsensusParams p;
    p.contig = true;
    std::mt19937 rng(15);
    std::string molecule;
    for (int i = 0; i < 400; ++i) molecule += "ACGT"[rng() & 3];
    std::vector<std::string> forward, backward;
    for (int start : {0, 30, 60, 250, 280, 310}) {
        forward.push_back(molecule.substr(static_cast<size_t>(start), 90));
    }
    backward.assign(forward.rbegin(), forward.rend());
    auto a = assemble_group(group(forward), p);
    auto b = assemble_group(group(backward), p);
    REQUIRE(a.size() == b.size());
    for (size_t i = 0; i < a.size(); ++i) CHECK(a[i].seq == b[i].seq);
}

TEST_CASE("a contig column is judged on the reads that reach it") {
    // Position 300 is covered by two reads out of six. A minor allele there is a minority of two,
    // not of six, and using the group size would make every contig edge look polymorphic.
    ConsensusParams p;
    p.contig = true;
    std::mt19937 rng(16);
    std::string molecule;
    for (int i = 0; i < 400; ++i) molecule += "ACGT"[rng() & 3];
    std::vector<std::string> seqs;
    for (int start = 0; start + 90 <= 400; start += 60) {
        seqs.push_back(molecule.substr(static_cast<size_t>(start), 90));
    }
    auto out = assemble_group(group(seqs), p);
    REQUIRE(out.size() == 1);
    CHECK(out[0].seq.size() == molecule.size() - (molecule.size() - 90) % 60);
    CHECK(out[0].seq == molecule.substr(0, out[0].seq.size()));
}
