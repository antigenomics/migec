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

#include <algorithm>
#include <array>
#include <cstdint>
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

// ---------------------------------------------------------------------------------------------
// How full the barcode space is, and what that costs. Every field here is arithmetic on the
// composition and the observed count -- no fitting -- but it is the arithmetic that decides
// whether a molecule count means anything, so it is reported rather than left to the reader.
//
// Molecules land in the barcode space independently, so occupancy per barcode is Poisson(lambda).
// What is observed is the number of *occupied* barcodes, and that is what pins lambda:
//
//     occupied = S (1 - e^-lambda)   =>   lambda = -ln(1 - observed/S),   molecules = S * lambda
//
// This is the birthday problem in its useful form. The naive "expected collided pairs ~ n^2/2S" is
// its small-lambda limit and is badly wrong once the space is half full, which is exactly when
// somebody wants the number.
struct BarcodeSpace {
    // The nominal space is 4^L over the L *captured* positions. Fixed letters written between the
    // N runs are not part of it: `NNNNtNNNNtNNNN` captures 12 bases, not 14, and its nominal space
    // is 4^12. The t's are scored pattern positions like any other constant.
    int length = 0;                 // L, the number of captured positions
    double nominal_space = 0.0;     // 4^L, if the synthesiser were perfect
    // ...which it is not. An oligo synthesiser mixing "N" does not deliver 25% of each base, so
    // the usable space is always smaller than 4^L. `effective_space` is what the observed base
    // composition actually supports, and `bias_loss` is how much was lost to that skew.
    double effective_space = 0.0;   // 1 / prod_j m_j  <=  4^L
    double effective_length = 0.0;  // -sum_j log4(m_j)  <=  L
    double bias_loss = 0.0;         // 1 - effective_space / nominal_space
    uint64_t observed = 0;          // distinct barcodes seen
    double occupancy = 0.0;         // observed / effective_space
    double lambda = 0.0;            // molecules per barcode
    double molecules = 0.0;         // S * lambda: the collision-corrected molecule count
    double hidden = 0.0;            // molecules - observed: molecules no barcode reports
    // P(a barcode holds more than one molecule | it holds at least one). This is the fraction of
    // MIGs that are really two or more molecules pooled, and it is the number that says whether a
    // consensus per barcode is meaningful at all.
    double p_multi = 0.0;
    // Beyond this occupancy the estimates above stop being estimates: S is inferred from the
    // observed barcodes, so as the space fills, `molecules` collapses onto `observed` and would
    // report "no collisions" for the most collided library possible.
    bool saturated = false;
};

BarcodeSpace barcode_space(const UmiComposition& comp, uint64_t observed_barcodes,
                           double saturation = 0.9);

// ---------------------------------------------------------------------------------------------
// What the barcode error rate *should* be, from first principles, next to what was estimated.
//
// Two independent processes put a substitution in a barcode, and they are predicted by things we
// already know rather than fitted:
//
//   sequencing -- the reported Phred, averaged as the mean of 10^(-Q/10) over the barcode bases.
//                 Not 10^(-mean Q/10): the low-Q tail dominates the error and averaging Q first
//                 hides it. On a 2-colour instrument the nominal value is itself unreliable, which
//                 is why the measured calibration table exists -- see MatchParams.
//   polymerase -- eps_pol per base per cycle, over the cycles that matter.
//
// Their sum is what the distance-1 estimator has to reproduce. When it does not, one of them is
// wrong, and the ratio says which way.
struct ErrorBudget {
    double from_phred = 0.0;     // mean 10^(-Q/10) over observed barcode bases
    double mean_phred = 0.0;     // the (misleading) arithmetic mean of Q, for contrast
    double from_polymerase = 0.0;  // eps_pol * cycles
    double predicted = 0.0;      // the two together
    double estimated = 0.0;      // what estimate_umi_error() found in the data
    double ratio = 0.0;          // estimated / predicted
    // Fraction of barcodes expected to carry at least one error: 1 - (1 - predicted)^L.
    double barcodes_with_error = 0.0;
    // The distance-1 estimator subtracts the independent-pair expectation from the observed pair
    // count. Once most of a barcode's 3L neighbours are themselves real barcodes, that is a small
    // difference of two large numbers and the estimate collapses towards zero. Measured on
    // simulated data: 0.9x of truth at 4% occupancy, 0.25x at 50%, 0.001x at 93%.
    double neighbour_occupancy = 0.0;  // fraction of the 3L shell expected to be occupied
    bool estimate_unreliable = false;
};

ErrorBudget error_budget(const UmiComposition& comp, const std::array<uint64_t, 61>& phred_counts,
                         double estimated, uint64_t observed_barcodes,
                         double polymerase_error = 1e-5, int pcr_cycles = 25);

// Observed UMI counts. Keys are packed barcodes (see types.hpp).
//
// Storage is a bounded append buffer that is periodically sorted and run-length reduced into a
// sorted (key, count) array, NOT a hash map. On a deeply sequenced run this is the single largest
// allocation in the process, and the difference is not academic: libstdc++'s
// unordered_map<uint64_t,uint32_t> costs ~48 bytes per distinct UMI once the node, the cached
// hash and the bucket array are counted, against 16 here. At 4e8 distinct UMIs -- an ordinary
// NovaSeq output at 5 reads per molecule -- that is 19 GB against 6.4 GB.
//
// Sorted order is not a side effect, it is the point: it is what the range partition and the
// neighbourhood search both want, and it makes the whole structure a flat scan instead of a
// pointer chase.
//
// The ceiling is real and is not fixed here. 16 bytes x (distinct UMIs) still does not fit a
// laptop for a full run, and the answer is the range partition -- process one bucket of the
// barcode space at a time, so this object only ever holds 1/2^bits of the library. That lands
// with `.mig` bucket output in M2; until then `memory_bytes()` is reported on every run and
// `checkout` warns when it goes past a threshold.
class UmiCounts {
public:
    struct Entry {
        uint64_t key;
        uint32_t count;
    };

    // `buffer_umis` is the CEILING on the unsorted append buffer, not its size. The buffer grows
    // with the data -- half the distinct count, so merges stay amortised O(1) per add -- and a
    // fixed-size buffer would otherwise cost that ceiling per sample whatever the sample holds,
    // which on a 96-plex sheet is gigabytes of empty space.
    explicit UmiCounts(int umi_length, size_t buffer_umis = 1u << 20)
        : length_(umi_length), buffer_limit_(buffer_umis ? buffer_umis : 1) {
        flush_at_ = std::min<size_t>(buffer_limit_, kMinBuffer);
    }

    void add(uint64_t packed, uint32_t n = 1) {
        buf_.push_back(Entry{packed, n});
        total_ += n;
        if (buf_.size() >= flush_at_) flush();
    }

    // Folds another counter in. Used to combine per-thread accumulators; the result does not
    // depend on the order the threads finished in.
    void merge(const UmiCounts& other);

    size_t distinct() const { flush(); return entries_.size(); }
    uint64_t total() const { return total_; }
    int length() const { return length_; }

    // Sorted by key, ascending. Flushes first, so it is O(n log n) on the first call after adds
    // and free afterwards.
    const std::vector<Entry>& entries() const { flush(); return entries_; }

    // Read count for a packed barcode, or nullptr. Binary search over the sorted array: the top
    // levels stay resident, which is why this beats a hash probe here despite the log factor.
    const uint32_t* find(uint64_t key) const;

    // Resident bytes actually held by this counter, buffer included.
    size_t memory_bytes() const;

    CoverageHistogram histogram() const;
    // `weight_by_reads` draws the composition MIGEC calls `pwm.txt` (each UMI counted once per
    // read); false gives `pwm-units.txt` (each distinct UMI counted once). The unit-weighted one
    // is what you want for the collision arithmetic -- read weighting lets one huge MIG dictate
    // the composition.
    UmiComposition composition(bool weight_by_reads = false) const;

private:
    // const because every accessor needs it and none of them change what the object *means*.
    void flush() const;

    static constexpr size_t kMinBuffer = 4096;

    int length_;
    size_t buffer_limit_;
    mutable size_t flush_at_ = kMinBuffer;
    uint64_t total_ = 0;
    mutable std::vector<Entry> buf_;
    mutable std::vector<Entry> entries_;
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
    // Note: This gate is what makes a singleton-vs-singleton merge impossible, so it is lifted when
    // payload evidence is available: two barcodes whose reads agree on the molecule are one
    // molecule whatever their counts are.
    double max_child_fraction = 0.5;
    // Per-base disagreement expected between two draft consensuses of the SAME molecule. Two
    // reads each carrying error e disagree at ~2e.
    double payload_error = 2e-3;
    // A payload pair is "the same sequence" below this mismatch fraction. Used only to measure
    // how often two *unrelated* barcodes carry the same sequence anyway -- clonality -- which is
    // what decides how much payload agreement is worth.
    double payload_same_fraction = 0.05;
    // Random barcode pairs sampled to measure that. 0 disables payload evidence.
    int payload_null_samples = 20000;
    // Threads for the neighbourhood scan; 0 means one per core.
    //
    // Never: this must not change a single merge. The scan is a pure function of the barcode table
    // and the evidence -- it reads no union-find state -- so it parallelises exactly, and the
    // decisions it produces are APPLIED serially afterwards in the same smallest-first order as
    // before. Splitting the apply as well would be a different algorithm: merges chain, and which
    // root a child lands on depends on which merges happened before it.
    int threads = 0;
};

// What a barcode's own reads say about it, beyond how many there are. Both fields are optional
// and both exist because the count ratio -- the only evidence the first version used -- carries
// nothing on a library sequenced at 1-3 reads per UMI, which is the common case rather than the
// exotic one. Measured in `scripts/correction_accuracy.py`: recall 0.80 at 3.1 reads/UMI, 0.02 at
// 1.1. A parent with 2 reads and a child with 1 is not an asymmetry, and two singletons are not
// one either.
//
// Both are laid out flat and indexed in parallel with UmiCounts::entries().
struct BarcodeEvidence {
    // Mean error probability at each barcode position, over that barcode's reads: entry i,
    // position j is at [i * length + j]. A sequencing miscall in the barcode carries a LOW Phred
    // at the base it changed, and an early-PCR child carries a high one in every read -- so this
    // separates the two mechanisms the global rate has to average over. Works at one read.
    std::vector<float> position_error;
    // Draft payload consensus per barcode, `payload_width` bases each, entry i at
    // [i * payload_width]. A barcode error child is a read of the PARENT'S molecule, so its
    // payload matches; an independent molecule one substitution away has its own. Works at one
    // read, and it is the only thing that does when both barcodes are singletons.
    std::vector<char> payload;
    int payload_width = 0;

    bool has_quality() const { return !position_error.empty(); }
    bool has_payload() const { return payload_width > 0 && !payload.empty(); }
};

struct CorrectionResult {
    // Both vectors are indexed in parallel with UmiCounts::entries(), which is what keeps this
    // O(12 bytes) per distinct UMI instead of two hash maps at ~48 each.
    //
    // `root[i]` is the index of the barcode entry i was folded into, or i itself when entry i is
    // a root. Chains are already resolved, so root[root[i]] == root[i].
    std::vector<uint32_t> root;
    // Read count after correction. A barcode that was merged away has 0 and its reads appear in
    // its root's count.
    std::vector<uint32_t> corrected;
    double estimated_error = 0.0;   // per-base UMI error actually used
    size_t merged = 0;              // number of distinct UMIs folded into a parent
    uint64_t merged_reads = 0;
    // Molecules observed after correction, and the collision-corrected estimate. Two molecules of
    // the same sequence sharing a UMI are undetectable, so the raw count is biased low:
    //     M_hat = S_eff * -ln(1 - M_obs/S_eff)
    size_t molecules_observed = 0;
    double molecules_corrected = 0.0;
    bool saturated = false;  // observed UMIs are a large fraction of the usable space
    // Fraction of random barcode pairs whose payloads agree anyway -- the library's clonality.
    // Payload agreement is worth log(1/clonality): decisive in a diverse repertoire, worth
    // nothing in a clonal one, and this says which this library is.
    double payload_clonality = 0.0;
    size_t merged_by_payload = 0;  // merges the count ratio alone would have refused
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

CorrectionResult correct_umis(const UmiCounts& counts, const CorrectionParams& params = {},
                              const BarcodeEvidence& evidence = {});

// Position of `key` in counts.entries(), or SIZE_MAX. For tests and for translating a barcode
// into the index space CorrectionResult uses.
size_t index_of(const UmiCounts& counts, uint64_t key);

}  // namespace migec

#endif  // MIGEC_UMI_STATS_HPP
