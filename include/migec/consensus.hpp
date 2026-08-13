// The consensus of one molecule's reads, and the quality that can honestly be claimed for it.
//
// This is the scientific claim of the whole pipeline, so the two numbers that bound it are
// measured rather than assumed, and both are recorded next to the constant they set:
//
//  * `rt_floor` -- an error made during reverse transcription or the first PCR cycle is in every
//    read of the molecule, and no consensus removes it. X2 measured it at 1.54e-4 on an HIV-1
//    Primer ID control (docs/quality_floor.rst), which caps every emitted quality at ~Q38. The
//    1e-6 that a first-pass design assumed is excluded by two orders of magnitude.
//
//  * `linkage_threshold` -- splitting a MIG into two consensuses is accepted on the strength of
//    co-segregation, not on a count of polymorphic sites. X3 measured the false-positive curve by
//    randomising the reads x positions minor-allele matrix while preserving BOTH margins
//    (docs/nulls.rst): the 1% point is at a Bonferroni'd -log10 p of 8.68, where the nominal
//    p < 0.01 the Poisson derivation gives calls 19x too many.
//
// Substitutions only, everywhere. Reads are already oriented and trimmed by checkout, so a group
// is ungapped and left-anchored; the group's width is the shortest read in it.

#ifndef MIGEC_CONSENSUS_HPP
#define MIGEC_CONSENSUS_HPP

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace migec {

struct ConsensusParams {
    // Measured in X2. Also the cap on emitted quality: -10 log10(rt_floor).
    double rt_floor = 1.0e-4;
    // Measured in X3. Raise it to split less.
    double linkage_threshold = 8.68;
    // A group smaller than this cannot show co-segregation worth testing.
    uint32_t min_split_reads = 6;
    // Positions considered for a split, strongest minor allele first.
    int max_split_columns = 8;
    // Reads needed on a minor allele before the position is callable at all.
    uint32_t min_minor_reads = 2;
    // Reads a subclone must keep to be emitted as its own consensus.
    uint32_t min_subclone_reads = 2;
    uint8_t min_quality = 2;
    // Reads sharing a barcode tile the molecule instead of starting at the same base, so they
    // have to be placed against each other and cut into overlap components first. X1 measured
    // this on 10x 3' GEX (docs/fragmented.rst): the co-terminal assumption is false for 92% of
    // groups, and 27.3% of them hold more than one component -- one consensus over those asserts
    // sequence across a gap no read covers.
    bool contig = false;
    // Exact seed length used to place two reads against each other, and the evidence needed.
    int seed_length = 15;
    int min_seed_votes = 3;
    int min_overlap = 20;
    // Measured error frequency per reported Phred, from the .mig header. Empty means use
    // 10^(-q/10) -- which on a 2-colour instrument is wrong by an order of magnitude, so prefer
    // the measured table wherever there is one.
    std::vector<float> calibration;
};

struct ConsensusRead {
    std::string_view seq;
    std::string_view qual;
    // Where this read starts in the group's coordinate frame. Zero for an amplicon, where
    // checkout's pattern has already anchored every read at the same base. Random-primed reads
    // tagged with one UMI tile the molecule instead, and `place_reads` works the offsets out.
    int offset = 0;
};

struct Consensus {
    std::string seq;
    std::string qual;
    uint32_t reads = 0;
    // Which overlap component of the group this is, and how many the group had. In amplicon mode
    // always 0 of 1. More than one means the barcode's reads did not all reach each other, so
    // these are separate contigs of the same molecule -- or, on a saturated barcode, of different
    // molecules, which is why contig mode reports the occupancy it was run at.
    uint32_t component = 0;
    uint32_t components = 1;
    // -log10 of the strongest co-segregation found in this group, after it was split. Reported so
    // that a borderline call is visible rather than being silently resolved by the threshold.
    double linkage = 0.0;
    // Mean posterior error over the emitted bases, before the floor is added. This is what the
    // consensus itself achieved; the floor is what the chemistry costs on top.
    double mean_error = 0.0;
};

// The strongest co-segregation of minor alleles over any pair of callable positions, as a
// -log10 p with a Bonferroni correction for the pairs tested. 0 when there is nothing to test.
double linkage_score(const std::vector<ConsensusRead>& reads, const ConsensusParams& params);

// Places every read against the others by exact seed matching and returns one vector of reads per
// overlap component, offsets normalised so each component starts at 0. Components are returned in
// descending size. A component is NEVER extended across a gap: two reads that share a barcode but
// no sequence are two contigs, not one.
std::vector<std::vector<ConsensusRead>> place_reads(const std::vector<ConsensusRead>& reads,
                                                    const ConsensusParams& params);

// One group of reads -> one consensus per molecule found in it. Never empty for a non-empty
// group: a group that fails to split comes back as a single consensus.
std::vector<Consensus> assemble_group(const std::vector<ConsensusRead>& reads,
                                      const ConsensusParams& params);

}  // namespace migec

#endif  // MIGEC_CONSENSUS_HPP
