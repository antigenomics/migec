// UMI statistics: coverage histogram, base composition and entropy, and count correction.
//
// These are the numbers that tell you whether a library is usable before you spend an hour
// assembling it, and they are also the inputs to the correction decision.
//
// A note on which entropy. Per-position Shannon entropy H_j is what a sequence logo draws, and it
// is the right thing to *display*. It is the wrong thing to compute a collision rate from: the
// probability that two independent molecules draw the same UMI is the Renyi entropy of order 2,
//
//     P_collision = sum_u p_u^2 = prod_j m_j       where m_j = sum_a p_j(a)^2
//
// and since H_2 <= H_1, using Shannon overestimates the usable barcode space and *underestimates*
// collisions -- the direction that silently merges distinct molecules. Both are reported; only the
// second is used for any decision.

#ifndef MIGEC_UMI_STATS_HPP
#define MIGEC_UMI_STATS_HPP

#include <array>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace migec {

// Reads per UMI, binned by powers of two: bin i covers [2^i, 2^(i+1)). MIGEC used 17 such bins and
// the convention is worth keeping because published histograms are drawn on it.
struct CoverageHistogram {
    static constexpr int kBins = 17;
    std::vector<uint64_t> reads;  // reads falling in MIGs of this size
    std::vector<uint64_t> units;  // distinct UMIs of this size

    CoverageHistogram() : reads(kBins, 0), units(kBins, 0) {}

    uint64_t total_reads() const;
    uint64_t total_units() const;
    double mean_reads_per_umi() const;
    // Fraction of reads sitting in MIGs of at least `min_size`. Monotone in depth and directly
    // interpretable, unlike MIGEC's "is there a peak" test which inverts on deeply sequenced
    // libraries.
    double reads_in_migs_at_least(uint32_t min_size) const;
    // True when the library is over-sequenced enough for consensus to mean anything.
    bool over_sequenced() const { return mean_reads_per_umi() >= 5.0; }
};

// Per-position base composition. Rows are positions, columns are A,C,G,T.
struct UmiComposition {
    int length = 0;
    std::vector<std::array<double, 4>> freq;  // normalised per position

    // Shannon entropy of position j, in bits. 2.0 for a uniform position.
    double entropy(int j) const;
    // Information content of position j: 2 - H_j. This is the logo's letter height.
    double information(int j) const;
    double total_entropy() const;      // sum_j H_j
    double total_information() const;  // 2L - sum_j H_j

    // Collision probability of position j: m_j = sum_a p_j(a)^2. 0.25 when uniform.
    double collision(int j) const;
    // Effective barcode length in bases: -sum_j log4(m_j). Equals the real length when uniform,
    // and is always <= it.
    double effective_length() const;
    // Usable barcode space: 1 / prod_j m_j.
    double effective_space() const;
    // Expected number of UMIs shared by two or more distinct molecules, for `n_molecules`:
    //     ~ n^2/2 * prod_j m_j     (valid while n * max_u p_u << 1)
    double expected_collisions(double n_molecules) const;
};

// Observed UMI counts. Keys are packed barcodes (see types.hpp).
class UmiCounts {
public:
    explicit UmiCounts(int umi_length) : length_(umi_length) {}

    void add(uint64_t packed, uint32_t n = 1) { counts_[packed] += n; }
    size_t distinct() const { return counts_.size(); }
    uint64_t total() const;
    int length() const { return length_; }
    const std::unordered_map<uint64_t, uint32_t>& map() const { return counts_; }

    CoverageHistogram histogram() const;
    // `weight_by_reads` draws the composition MIGEC calls `pwm.txt` (each UMI counted once per
    // read); false gives `pwm-units.txt` (each distinct UMI counted once). The unit-weighted one
    // is what you want for the collision arithmetic -- read weighting lets one huge MIG dictate
    // the composition.
    UmiComposition composition(bool weight_by_reads = false) const;

private:
    int length_;
    std::unordered_map<uint64_t, uint32_t> counts_;
};

struct CorrectionParams {
    // Per-base sequencing error in the UMI. Negative means "estimate from the data".
    double sequencing_error = -1.0;
    // Per-base polymerase error per PCR cycle, and the number of cycles. This is the component a
    // sequencing-only model misses: a substitution introduced in cycle 1-3 is present in ~50/25/12%
    // of the descendants and carries HIGH quality in every read, so a Poisson on the sequencing
    // rate assigns it essentially zero probability and it survives as a second molecule. That is
    // the dominant residual error in UMI counting.
    double polymerase_error = 1e-5;
    int pcr_cycles = 25;
    // Posterior above which a child is merged into its parent.
    double min_posterior = 0.95;
    // A child can never be larger than this fraction of its parent, whatever the posterior says.
    double max_child_fraction = 0.5;
};

struct CorrectionResult {
    // child packed UMI -> parent packed UMI. Chains are resolved, so every value is a root.
    std::unordered_map<uint64_t, uint64_t> parent;
    std::unordered_map<uint64_t, uint32_t> corrected;  // packed UMI -> corrected read count
    double estimated_error = 0.0;   // per-base UMI error actually used
    size_t merged = 0;              // number of distinct UMIs folded into a parent
    uint64_t merged_reads = 0;
    // Molecules observed after correction, and the collision-corrected estimate. Two molecules of
    // the same sequence sharing a UMI are undetectable, so the raw count is biased low:
    //     M_hat = S_eff * -ln(1 - M_obs/S_eff)
    size_t molecules_observed = 0;
    double molecules_corrected = 0.0;
    bool saturated = false;  // observed UMIs are a large fraction of the usable space
};

// Estimates the per-base UMI error rate from the excess of 1-mismatch neighbours over what
// independent draws would produce. The sibling term matters: two children of the same parent
// differing at the same position by different bases are themselves at distance 1, and counting
// them as parent-child pairs inflates the estimate by up to 2x.
//
//     E[D1](eps) = C(n,2) * P_collision
//                + 3L * sum_i (1 - exp(-c_i eps))         parent-child
//                + 3L * sum_i (1 - exp(-c_i eps))^2       sibling
double estimate_umi_error(const UmiCounts& counts, const UmiComposition& comp);

CorrectionResult correct_umis(const UmiCounts& counts, const CorrectionParams& params = {});

}  // namespace migec

#endif  // MIGEC_UMI_STATS_HPP
