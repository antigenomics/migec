#include "migec/pattern.hpp"

#include <algorithm>
#include <cmath>

#include "migec/types.hpp"

namespace migec {
namespace {

// log2(4e/3) for a mismatch, and the matched-position term, both as functions of the error
// probability.
inline double mismatch_bits(double e) {
    // Clamp: e = 0 would give -inf and reject a read on one perfect-quality mismatch, which is
    // exactly the over-confidence the calibration table exists to prevent.
    if (e < 1e-6) e = 1e-6;
    return std::log2(4.0 * e / 3.0);
}

inline double match_bits(double e, int m) {
    // P(observe a member of S | tag is here) = (1-e)/m + (m-1)e/(3m):
    // either the true base is this one and was read correctly, or it was another member of S and
    // miscalled into this one.
    const double p = (1.0 - e) / m + (m - 1) * e / (3.0 * m);
    return std::log2(4.0 * p);
}

// Both of the above are log2 calls, which at ~25 ns each dominated everything: the scan touches a
// few thousand positions per read, so the transcendental was 90% of checkout's runtime. The score
// only ever depends on the reported Phred and the size of the IUPAC set, and both are small
// integers, so the whole thing tabulates into 1.2 kB.
struct ScoreTable {
    float mismatch[kMaxPhred + 1];
    float match[5][kMaxPhred + 1];  // indexed by |S| in 1..4; row 0 unused
};

ScoreTable build_table(const std::vector<double>& calibration) {
    ScoreTable t{};
    for (int q = 0; q <= kMaxPhred; ++q) {
        const double e = (!calibration.empty() && static_cast<size_t>(q) < calibration.size())
                             ? calibration[q]
                             : phred_error(static_cast<uint8_t>(q));
        t.mismatch[q] = static_cast<float>(mismatch_bits(e));
        for (int m = 1; m <= 4; ++m) t.match[m][q] = static_cast<float>(match_bits(e, m));
    }
    return t;
}

const ScoreTable& nominal_table() {
    static const ScoreTable t = build_table({});
    return t;
}

}  // namespace

BarcodePattern BarcodePattern::compile(std::string_view spec) {
    if (spec.empty()) throw MigecError("pattern: empty pattern");
    BarcodePattern p;
    p.spec_ = std::string(spec);
    p.mask_.reserve(spec.size());
    p.weight_.reserve(spec.size());
    p.capture_.reserve(spec.size());

    for (char c : spec) {
        if (c == 'N' || c == 'n') {
            p.mask_.push_back(0);
            p.weight_.push_back(0.0f);
            p.capture_.push_back(1);
            ++p.umi_length_;
            continue;
        }
        if (c == 'X' || c == 'x') {
            p.mask_.push_back(0);
            p.weight_.push_back(0.0f);
            p.capture_.push_back(2);
            ++p.cell_length_;
            continue;
        }
        if (c == '.') {
            p.mask_.push_back(0);
            p.weight_.push_back(0.0f);
            p.capture_.push_back(0);
            continue;
        }
        const uint8_t m = iupac_mask(c);
        if (m == 0) {
            throw MigecError(std::string("pattern: '") + c +
                             "' is not a IUPAC symbol, N, or '.' in \"" + std::string(spec) + "\"");
        }
        p.mask_.push_back(m);
        // Lowercase marks the fuzzy region (the adapter), where a mismatch is expected. It is a
        // weight, not a gate: MIGEC treated lowercase as "matches anything", which threw away the
        // evidence entirely.
        p.weight_.push_back(c >= 'a' && c <= 'z' ? 0.5f : 1.0f);
        p.capture_.push_back(0);
        ++p.scored_;
    }

    if (p.scored_ == 0) {
        throw MigecError("pattern: \"" + std::string(spec) +
                         "\" has no scored position, so it matches everywhere");
    }
    // Checked here rather than where the barcode is packed: pack_barcode runs on a worker thread
    // per read, and a pattern is compiled once, on the caller's thread, where the error is
    // attributable to the row of the barcode table that caused it.
    if (p.umi_length_ > kMaxBarcodeLen) {
        throw MigecError("pattern: \"" + std::string(spec) + "\" captures a " +
                         std::to_string(p.umi_length_) + " nt UMI; the packed representation holds " +
                         std::to_string(kMaxBarcodeLen));
    }
    if (p.cell_length_ > kMaxBarcodeLen) {
        throw MigecError("pattern: \"" + std::string(spec) + "\" captures a " +
                         std::to_string(p.cell_length_) +
                         " nt cell barcode; the packed representation holds " +
                         std::to_string(kMaxBarcodeLen));
    }
    return p;
}

double BarcodePattern::default_min_score(size_t read_length, size_t n_patterns,
                                         double alpha) const {
    const size_t n_offsets =
        read_length >= mask_.size() ? read_length - mask_.size() + 1 : 1;
    return std::log2(static_cast<double>(n_offsets) * static_cast<double>(n_patterns) / alpha);
}

PatternMatch BarcodePattern::match(std::string_view seq, std::string_view qual,
                                   const MatchParams& params) const {
    PatternMatch out;
    const size_t plen = mask_.size();
    if (seq.size() < plen) return out;
    if (!qual.empty() && qual.size() != seq.size()) {
        throw MigecError("pattern_match: sequence and quality lengths differ");
    }

    size_t last_offset = seq.size() - plen;
    if (params.max_offset >= 0) {
        last_offset = std::min<size_t>(last_offset, static_cast<size_t>(params.max_offset));
    }

    const double min_score = params.min_score >= 0.0
                                 ? params.min_score
                                 : default_min_score(seq.size());
    const auto& calib = params.quality_calibration;
    // ponytail: the calibrated path rebuilds its table per read (~300 log2, ~8 us). Nothing
    // supplies a calibration yet -- it arrives from the .mig header in M2 -- and when it does the
    // table gets built once at that boundary rather than here.
    ScoreTable local;
    const ScoreTable* tab = &nominal_table();
    if (!calib.empty()) {
        local = build_table(calib);
        tab = &local;
    }
    const bool have_qual = !qual.empty();

    double best = -1e300, second = -1e300;
    long best_off = -1;

    // An offset scoring below this can be neither the winner nor the runner-up that sets the
    // margin, so it can be abandoned the moment it is out of reach -- which is every wrong offset,
    // within a handful of positions.
    const double prune_floor = min_score - params.min_margin;
    // ...but the bar has to leave room for the runner-up. An offset that cannot reach the incumbent
    // best can still land within min_margin of it, which makes the placement ambiguous; pruning at
    // `best` drops it silently and the margin then comes back as `best - (-inf)`, so a read with
    // two near-equal placements is reported as an unambiguous match at whichever came first.
    const double margin_slack = params.min_margin > 0.0 ? params.min_margin : 0.0;

    // ponytail: plain O(offsets x pattern) scan, made cheap by the table above and the early exit
    // below rather than by a smarter algorithm. Upgrade path if a long read ever makes it matter:
    // shift-or over the longest uppercase run to generate candidate offsets, then score only those.
    for (size_t off = 0; off <= last_offset; ++off) {
        const double reachable = best - margin_slack;
        const double bar = reachable > prune_floor ? reachable : prune_floor;
        double s = 0.0;
        bool pruned = false;
        for (size_t i = 0; i < plen; ++i) {
            const uint8_t m = mask_[i];
            if (m == 0) continue;  // UMI or wildcard: never scored
            const uint8_t code = base_code(seq[off + i]);
            const uint8_t q = have_qual ? phred_from_char(qual[off + i]) : 30;
            if (code != kInvalidBase && ((m >> code) & 1u)) {
                s += tab->match[iupac_size(m)][q];
            } else {
                // An N in the read is not evidence against the tag; it is no evidence at all.
                // Scoring it as a mismatch at its (low) quality already handles that, and an
                // uncallable base carries Q2, worth -0.6 bits.
                s += weight_[i] * tab->mismatch[q];
            }
            // Once this offset cannot reach the bar even with every remaining position matching
            // perfectly, stop.
            if (s + 2.0 * static_cast<double>(plen - i - 1) < bar) {
                pruned = true;
                break;
            }
        }
        // A pruned offset carries a partial score, which must never be mistaken for a real one:
        // it is only known to be below the bar, and the bar can sit above the incumbent best.
        if (pruned) continue;
        if (s > best) {
            second = best;
            best = s;
            best_off = static_cast<long>(off);
        } else if (s > second) {
            second = s;
        }
    }

    if (best_off < 0 || best < min_score) return out;
    const double margin = second <= -1e299 ? best : best - second;
    if (margin < params.min_margin) return out;  // ambiguous placement

    out.found = true;
    out.offset = static_cast<int>(best_off);
    out.score = best;
    out.margin = margin;
    out.payload_begin = static_cast<int>(best_off + plen);
    out.umi.reserve(static_cast<size_t>(umi_length_));
    out.umi_qual.reserve(static_cast<size_t>(umi_length_));
    out.cell.reserve(static_cast<size_t>(cell_length_));
    out.cell_qual.reserve(static_cast<size_t>(cell_length_));
    for (size_t i = 0; i < plen; ++i) {
        const uint8_t what = capture_[i];
        if (!what) continue;
        const char base = seq[static_cast<size_t>(best_off) + i];
        const char q = qual.empty() ? char_from_phred(30) : qual[static_cast<size_t>(best_off) + i];
        if (what == 1) {
            out.umi.push_back(base);
            out.umi_qual.push_back(q);
        } else {
            out.cell.push_back(base);
            out.cell_qual.push_back(q);
        }
    }
    return out;
}

void BarcodePattern::calibrate(std::string_view seq, std::string_view qual, int offset,
                               std::vector<std::array<std::array<uint64_t, 2>, 61>>& by_position)
    const {
    if (offset < 0) return;
    const size_t begin = static_cast<size_t>(offset);
    if (begin + mask_.size() > seq.size()) return;
    for (size_t i = 0; i < mask_.size(); ++i) {
        const uint8_t m = mask_[i];
        // Unambiguous scored positions only. A degenerate position tolerates more than one base,
        // so a disagreement there is not a miscall with probability 1 and would bias the table.
        if (m == 0 || iupac_size(m) != 1) continue;
        const size_t at = begin + i;
        const uint8_t obs = base_code(seq[at]);
        if (obs == kInvalidBase) continue;
        if (i >= by_position.size()) break;
        const uint8_t q = at < qual.size() ? phred_from_char(qual[at]) : kMaxPhred;
        ++by_position[i][q > 60 ? 60 : q][((m >> obs) & 1u) ? 0 : 1];
    }
}

void PatternSet::add(std::string sample_id, std::string_view spec) {
    samples_.push_back(std::move(sample_id));
    patterns_.push_back(BarcodePattern::compile(spec));
}

PatternSet::Assignment PatternSet::assign(std::string_view seq, std::string_view qual,
                                          const MatchParams& params) const {
    Assignment out;
    // Within one sample the margin is about placement; between samples it is about identity.
    // Score every sample without the placement margin, then apply the margin once, across
    // samples -- otherwise two samples whose tags differ by one base reject each other's reads.
    MatchParams per_pattern = params;
    per_pattern.min_margin = 0.0;

    double best = -1e300, second = -1e300;
    for (size_t i = 0; i < patterns_.size(); ++i) {
        // Per pattern, because the Bonferroni bound counts the offsets *this* pattern is scanned
        // over, and patterns in one sheet need not be the same length.
        if (params.min_score < 0.0) {
            per_pattern.min_score = patterns_[i].default_min_score(seq.size(), patterns_.size());
        }
        PatternMatch m = patterns_[i].match(seq, qual, per_pattern);
        if (!m.found) continue;
        if (m.score > best) {
            second = best;
            best = m.score;
            out.sample = static_cast<int>(i);
            out.match = m;
        } else if (m.score > second) {
            second = m.score;
        }
    }
    if (out.sample < 0) return out;
    if (second > -1e299 && best - second < params.min_margin) {
        out.ambiguous = true;
        out.sample = -1;
        out.match = PatternMatch();
    }
    return out;
}

}  // namespace migec
