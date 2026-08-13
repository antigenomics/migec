// The consensus of one molecule's reads, and the quality that can honestly be claimed for it.
//
// This is the scientific claim of the whole pipeline, so the two numbers that bound it are
// measured rather than assumed, and both are recorded next to the constant they set:
//
//  * `rt_floor` -- an error made during reverse transcription or the first PCR cycle is in every
//    read of the molecule, and no consensus removes it. X2 measured it at 1.54e-4 on an HIV-1
//    Primer ID control (docs/quality_floor.rst) and by 10x for their V(D)J RT, which caps every
//    emitted quality at Q40. A blanket 1e-6 is excluded by two orders of magnitude for an RT
//    protocol -- and is right for a DNA one, which is why the floor is named per chemistry.
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
    // The RT/first-cycle-PCR error that is in every read of a molecule and that no consensus
    // removes. Also the cap on emitted quality: -10 log10(rt_floor), so 1e-4 caps at Q40.
    //
    // Never: this is the ONE-MOLECULE floor, and every record migec emits is one molecule. 10x
    // states it exactly: "The estimated error rate for the V(D)J RT reaction is 1e-4 per base.
    // Therefore, assembled bases that are covered by a single UMI are assigned Q40, and bases
    // covered by at least two UMIs are assigned Q60." The Q60 branch is available only after
    // several molecules agree -- an RT error is common-mode within a molecule and independent
    // between them -- and combining molecules is arda's job, not this one. A per-molecule record
    // claiming Q50 is claiming the two-UMI confidence on one-UMI evidence.
    //
    // X2's own measurement agrees with 10x's number: 1.54e-4 [1.36e-4, 1.74e-4] on SRR1763769,
    // itself an upper bound (that library is 49.6% occupied on a 9 nt barcode, so collided MIGs
    // pool two templates and the disagreement is charged to the RT).
    //
    //     1e-4   default: V(D)J / RepSeq off an ordinary RT, and 10x's stated figure (caps at Q40)
    //     1e-5   ctDNA and exome -- 2-10x less, and cfDNA assays often have no RT at all, only a
    //            first cycle of a high-fidelity polymerase (caps at Q50)
    //     1e-6   a DNA-only workflow with a high-fidelity polymerase and few pre-amplification
    //            cycles, where the read-out is a low-frequency variant (caps at Q60)
    double rt_floor = 1.0e-4;
    // Counting mode: emit the group's most frequent EXACT sequence, with each base carrying the
    // best quality any read of that sequence reported for it. No column model, so no per-base
    // error correction and no sub-clustering -- what it buys is molecule counts, fast. A read that
    // disagrees anywhere is a different string and votes for nothing, which is why this is the
    // wrong mode for error suppression and the right one for a count.
    bool fast = false;
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
    // Reads carrying the emitted sequence exactly. Set by the fast path, where it is the whole of
    // the evidence: `support` of `reads` agreed, and the rest voted for some other string. Zero in
    // the full path, where the consensus is a per-column call that need match no read at all.
    uint32_t support = 0;
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
