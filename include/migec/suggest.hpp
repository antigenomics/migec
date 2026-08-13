// suggest: read the barcode layout off the data instead of off the protocol.
//
// A barcode pattern is three kinds of cycle interleaved, and the per-cycle base composition tells
// them apart without being told anything:
//
//   UMI       the synthesiser mixed all four bases, so each is near 1/4 and the cycle carries
//             ~2 bits. This is the "1/4 PWM trace" -- four flat lines at 25%.
//   constant  one base at ~100%: the adapter, the primer, or a sample tag. ~0 bits.
//   payload   biological sequence: uneven, but not 1/4 and not fixed. Somewhere in between, and
//             the giveaway is that it never settles.
//
// The boundary between them is the pattern, so `suggest` prints one that can be pasted into a
// barcode table. This is also the honest way to handle a protocol description that disagrees with
// what the sequencer produced -- the file wins.

#ifndef MIGEC_SUGGEST_HPP
#define MIGEC_SUGGEST_HPP

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace migec {

struct CycleStats {
    int cycle = 0;
    std::array<uint64_t, 5> counts{};  // A, C, G, T, other
    uint64_t phred_sum = 0;
    uint64_t n = 0;

    std::array<double, 4> frequencies() const;
    double entropy() const;         // Shannon, bits: ~2 for a UMI cycle, ~0 for a constant one
    double collision() const;       // sum_a p_a^2: 1/4 for a uniform cycle, 1 for a fixed one
    char consensus() const;
    double consensus_fraction() const;
    double mean_phred() const;
    // Distance from a flat 1/4 trace, as total variation: 0 when perfectly uniform, 0.75 when
    // fixed. This is the statistic that names a UMI cycle, and it is bounded and scale-free in a
    // way entropy is not -- entropy 1.9 bits could be 25/25/25/25 slightly perturbed or a genuine
    // three-base mixture, and those are not the same thing.
    double deviation_from_uniform() const;
};

enum class CycleKind { kUmi, kConstant, kVariable };

struct CycleProfile {
    std::vector<CycleStats> cycles;
    uint64_t reads = 0;
    size_t read_length = 0;
};

struct Segment {
    CycleKind kind = CycleKind::kVariable;
    int begin = 0;
    int end = 0;  // exclusive
    std::string consensus;   // the literal bases, for a constant run
    double mean_deviation = 0.0;
    int length() const { return end - begin; }
};

struct Suggestion {
    CycleProfile profile;
    std::vector<Segment> segments;
    std::string pattern;      // paste-ready, MIGEC dialect
    int umi_length = 0;
    int anchor_length = 0;    // constant bases the pattern scores
    // Reads whose first `pattern.size()` cycles the suggested pattern would accept. Not a
    // guarantee, but a pattern that does not match its own training data is not worth printing.
    double matched_fraction = 0.0;
    std::string note;
};

// Profile the first `n_cycles` of up to `max_reads` reads.
CycleProfile profile_cycles(const std::string& fastq_path, int n_cycles = 60,
                            uint64_t max_reads = 200000);

// Segment a profile and build a pattern from it.
//
// `umi_deviation` is how far from a flat 1/4 trace a cycle may sit and still be called UMI. The
// default is loose enough for real synthesiser bias -- oligo mixes are routinely 20/30/30/20 --
// and tight enough to exclude biological sequence, which is never that even across a run of
// positions.
Suggestion suggest_pattern(const CycleProfile& profile, double umi_deviation = 0.18,
                           double constant_fraction = 0.9, int min_run = 3);

}  // namespace migec

#endif  // MIGEC_SUGGEST_HPP
