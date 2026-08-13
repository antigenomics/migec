// Whitelist correction: snap an observed barcode to the list of barcodes that were actually
// synthesised, when the evidence supports it.
//
// Never: The posterior needs a BACKGROUND hypothesis, and leaving it out is the failure this file
// exists to prevent. Without it the model says "the true barcode is one of these 737,000", so
// every observed barcode is assigned to *some* entry with posterior 1.0 -- an index-hopped read,
// an undeclared sample, free-floating ambient sequence, all of it silently absorbed into whichever
// whitelist entry happens to be nearest. With it, the competing explanation "this barcode is not
// on the list and was read correctly" is on the table and usually wins for those.
//
// The prior on that background is MEASURED, not assumed: barcodes at distance >= 2 from every
// whitelist entry cannot be single-substitution errors of anything on it, so the share of reads
// they carry is a lower bound on how much of the library is off-list.
//
// An `N` is expanded, never discarded. It is a base the instrument declined to call, so it is
// consistent with all four at e = 0.75 -- which is exactly what the per-position error term
// already encodes, and it means a barcode with one N is still correctable rather than thrown away.

#ifndef MIGEC_WHITELIST_HPP
#define MIGEC_WHITELIST_HPP

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace migec {

struct WhitelistParams {
    // Posterior a candidate must reach before the barcode is rewritten. High on purpose: a wrong
    // snap moves a molecule into another cell, and nothing downstream can tell.
    double min_posterior = 0.975;
    // Per-base error used when the observed barcode carries no quality string.
    double default_error = 1e-3;
    // Prior that the true barcode is NOT on the list. Negative means "measure it from the data".
    //
    // Note: This is a prior on THIS barcode, not on the library. The whitelist prior is a share of
    // (1 - background) spread over every entry, so a 737,000-entry list gives each candidate
    // ~1e-6; a background quoted as "1% of the library is off-list" would then be four orders of
    // magnitude larger than any candidate and would win every time. The comparable quantity is
    // the off-list read share DIVIDED BY the number of distinct off-list barcodes observed --
    // which is what `measure_background` computes.
    double background_prior = -1.0;
};

struct WhitelistStats {
    uint64_t barcodes = 0;          // distinct observed
    uint64_t exact = 0;             // already on the list
    uint64_t corrected = 0;         // snapped to a neighbour
    uint64_t off_list = 0;          // kept as they are: no candidate won
    uint64_t reads_corrected = 0;
    uint64_t far = 0;               // distance >= 2 from every entry -- what measures the prior
    double background_prior = 0.0;  // the value used, measured or given
};

// The list itself, packed and sorted for binary search. Loading is the only place a barcode
// length is fixed: every entry must be the same length, and a mismatch is an error rather than a
// silent truncation.
class Whitelist {
public:
    // One barcode per line; blank lines and `#` comments skipped. Accepts a plain or gzipped file,
    // and tolerates 10x's "-1" suffix.
    static Whitelist load(const std::string& path);

    bool empty() const { return keys_.empty(); }
    size_t size() const { return keys_.size(); }
    int length() const { return length_; }
    bool contains(uint64_t key) const;

    // The per-barcode background prior: the share of reads on barcodes that cannot be single
    // substitutions of anything on the list, spread over the distinct off-list barcodes seen.
    // Both halves are measured. `far_reads` and `far_barcodes` count barcodes at distance >= 2
    // from every entry; `off_barcodes` is every distinct barcode not exactly on the list.
    static double measure_background(uint64_t far_reads, uint64_t total_reads,
                                     uint64_t off_barcodes);

    // Best whitelist entry for `observed`, or empty when nothing beat the background. `counts`
    // is the observed read count of each whitelist entry, indexed as this class stores them, and
    // is what makes a heavily used barcode a likelier parent than an unused one.
    std::string correct(std::string_view observed, std::string_view qual,
                        const std::vector<uint32_t>& counts, const WhitelistParams& params,
                        double background_prior) const;

    size_t index_of(uint64_t key) const;

private:
    std::vector<uint64_t> keys_;  // sorted
    int length_ = 0;
};

}  // namespace migec

#endif  // MIGEC_WHITELIST_HPP
