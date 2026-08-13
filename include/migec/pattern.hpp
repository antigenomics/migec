// Barcode patterns: a degenerate adapter/sample tag with UMI positions, matched against a read.
//
// The grammar is MIGEC's, because the published barcode tables are written in it and we want them
// to keep working verbatim:
//
//     S1<TAB>aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
//
//   UPPERCASE   a scored position, matched exactly. IUPAC letters allowed (R = A|G, and so on),
//               which is what makes the tag "degenerate".
//   lowercase   a scored position with half weight -- the adapter region, where a mismatch is
//               expected and should not by itself reject the read.
//   N or n      a UMI position. Captured, never scored. Non-contiguous runs concatenate, so
//               NNNNtNNNNtNNNN yields one 12 nt UMI.
//   X or x      a CELL BARCODE position. Captured into a separate barcode, never scored. This is
//               the one extension to MIGEC's dialect, and `X` was chosen because it is not a
//               IUPAC symbol -- so no published MIGEC table can contain one and every existing
//               table keeps its exact meaning. `C` would have been the 10x convention and is a
//               base here, which would have silently reinterpreted real barcode tables.
//   .           a wildcard: neither scored nor captured.
//
// Acceptance is a quality-aware log-likelihood ratio, not a mismatch count. For a scored position
// with IUPAC set S of size m, observed base b and error probability e:
//
//     b in S:   log2( 4 * [ (1-e)/m + (m-1)e/(3m) ] )
//     b not in: w * log2( 4e/3 )                          w = 1.0 upper, 0.5 lower
//
// which is log2 P(base | the tag is here) / P(base | random sequence). At m=1 a match is worth
// +2.00 bits, a mismatch -9.55 bits at Q30 and -0.60 bits at Q2. So a mismatch on a bad base is
// nearly free and one on a good base is fatal -- which is what MIGEC's good/bad mismatch counting
// was reaching for, done continuously and without its two bugs (it indexed quality from the start
// of the read rather than the match offset, and a dangling `else` meant low-quality mismatches
// were never counted at all).

#ifndef MIGEC_PATTERN_HPP
#define MIGEC_PATTERN_HPP

#include <array>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

#include "migec/types.hpp"

namespace migec {

struct MatchParams {
    // Minimum score to accept, in bits. Default from the Bonferroni bound over the offsets
    // actually scanned; see default_min_score(). A negative value means "use the default".
    double min_score = -1.0;
    // The best offset must beat the runner-up by this much, else the placement is ambiguous.
    double min_margin = 5.0;
    // Search window. -1 searches the whole read; 0 anchors at position 0 (positional chemistries
    // like 10x); k allows the tag to start anywhere in [0, k].
    int max_offset = -1;
    // Nominal Phred is not the error rate on 2-colour instruments. If non-empty, this is indexed
    // by reported Phred and gives the measured error probability.
    std::vector<double> quality_calibration;
};

struct PatternMatch {
    bool found = false;
    int offset = -1;         // where the pattern starts in the read
    double score = 0.0;      // bits
    double margin = 0.0;     // bits over the runner-up offset
    int payload_begin = 0;   // first base after the matched pattern -- where trimming leaves you
    std::string umi;
    std::string umi_qual;
    std::string cell;        // empty unless the pattern has X positions
    std::string cell_qual;
};

class BarcodePattern {
public:
    // Throws MigecError on an empty pattern or one with no scored position.
    static BarcodePattern compile(std::string_view spec);

    PatternMatch match(std::string_view seq, std::string_view qual,
                       const MatchParams& params = MatchParams()) const;

    // Counts match/mismatch at this pattern's UNAMBIGUOUS scored positions, indexed by the
    // reported Phred. Those positions are known sequence -- the adapter and the sample tag -- so
    // a disagreement there is an instrument error and nothing else, which makes them the only
    // free calibration standard in the read. Degenerate (IUPAC) positions are skipped: a mismatch
    // against a 2-base set is not a miscall with probability 1.
    void calibrate(std::string_view seq, std::string_view qual, int offset,
                   std::vector<std::array<std::array<uint64_t, 2>, 61>>& by_position) const;

    size_t size() const { return mask_.size(); }
    int umi_length() const { return umi_length_; }
    int cell_length() const { return cell_length_; }
    int scored_positions() const { return scored_; }
    const std::string& spec() const { return spec_; }

    // Bonferroni over the offsets actually scanned, for a per-read false-match rate alpha:
    //     min_score = log2( n_offsets * n_patterns / alpha )
    // This is a starting point, not gospel -- reads are not i.i.d. uniform ACGT (shared primers,
    // composition bias), so calibrate against shuffled decoy patterns on real data.
    // `max_offset` must be the one the scan will actually use: anchoring the pattern at the read
    // start is what makes a short handle placeable at all, and charging it for offsets it never
    // tries would keep refusing it.
    double default_min_score(size_t read_length, size_t n_patterns = 1, double alpha = 0.01,
                             int max_offset = -1) const;

private:
    std::string spec_;
    std::vector<uint8_t> mask_;    // IUPAC mask per position, 0 for UMI/wildcard
    std::vector<float> weight_;    // 1.0 upper, 0.5 lower, 0 for unscored
    // 0 = not captured, 1 = UMI, 2 = cell barcode.
    std::vector<uint8_t> capture_;
    int umi_length_ = 0;
    int cell_length_ = 0;
    int scored_ = 0;
};

// One pattern per sample, as read from a barcode metadata table. Assignment picks the best-scoring
// sample and requires it to beat the runner-up sample by min_margin -- otherwise the read is
// ambiguous and is counted rather than arbitrarily assigned.
// One or two patterns per sample. The second -- MIGEC calls it the *slave* barcode, and it is
// column 3 of the published tables -- sits on the other mate, and its captured positions extend
// the UMI rather than starting a new one. That is how a 24 nt dual-end UMI is written as
// `NNNNNNNNNNNNtgact` and `agtcaNNNNNNNNNNNN`: twelve bases from each end of the molecule.
//
// Never: Both must match. A dual-end design that accepts a read on the master alone silently emits
// half-length UMIs alongside full-length ones, and every collision estimate downstream is then
// computed over two different barcode spaces at once.
class PatternSet {
public:
    void add(std::string sample_id, std::string_view spec, std::string_view slave = {});

    struct Assignment {
        int sample = -1;  // index into samples(); -1 = unassigned
        PatternMatch match;
        bool ambiguous = false;
    };

    Assignment assign(std::string_view seq, std::string_view qual,
                      const MatchParams& params = MatchParams()) const;

    size_t size() const { return patterns_.size(); }
    const std::vector<std::string>& samples() const { return samples_; }
    const BarcodePattern& pattern(size_t i) const { return patterns_[i]; }
    bool has_slave(size_t i) const { return slave_of_[i] >= 0; }
    const BarcodePattern& slave(size_t i) const { return slaves_[static_cast<size_t>(slave_of_[i])]; }
    // Total captured UMI length for a sample: master plus slave.
    int umi_length(size_t i) const {
        return patterns_[i].umi_length() + (has_slave(i) ? slave(i).umi_length() : 0);
    }

private:
    std::vector<std::string> samples_;
    std::vector<BarcodePattern> patterns_;
    std::vector<BarcodePattern> slaves_;
    std::vector<int> slave_of_;  // index into slaves_, or -1
};

}  // namespace migec

#endif  // MIGEC_PATTERN_HPP
