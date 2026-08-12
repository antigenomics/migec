#include "migec/pattern.hpp"

#include <algorithm>
#include <cmath>

#include "migec/types.hpp"

namespace migec {
namespace {

// log2(4e/3) for a mismatch, and the matched-position term, both as functions of the error
// probability. Tabulating by Phred is not enough once a calibration table is in play, so these
// stay as calls; they are a handful of instructions and the offset loop is the cost, not this.
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

}  // namespace

BarcodePattern BarcodePattern::compile(std::string_view spec) {
    if (spec.empty()) throw MigecError("pattern: empty pattern");
    BarcodePattern p;
    p.spec_ = std::string(spec);
    p.mask_.reserve(spec.size());
    p.weight_.reserve(spec.size());
    p.is_umi_.reserve(spec.size());

    for (char c : spec) {
        if (c == 'N' || c == 'n') {
            p.mask_.push_back(0);
            p.weight_.push_back(0.0f);
            p.is_umi_.push_back(1);
            ++p.umi_length_;
            continue;
        }
        if (c == '.') {
            p.mask_.push_back(0);
            p.weight_.push_back(0.0f);
            p.is_umi_.push_back(0);
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
        p.is_umi_.push_back(0);
        ++p.scored_;
    }

    if (p.scored_ == 0) {
        throw MigecError("pattern: \"" + std::string(spec) +
                         "\" has no scored position, so it matches everywhere");
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

    double best = -1e300, second = -1e300;
    long best_off = -1;

    // ponytail: plain O(offsets x pattern) scan. At 150 nt reads and a 40 nt pattern that is ~4400
    // scored positions per read, which is not the bottleneck next to gzip. Upgrade path if it ever
    // is: shift-or over the longest uppercase run to generate candidate offsets, then score only
    // those.
    for (size_t off = 0; off <= last_offset; ++off) {
        double s = 0.0;
        for (size_t i = 0; i < plen; ++i) {
            const uint8_t m = mask_[i];
            if (m == 0) continue;  // UMI or wildcard: never scored
            const char b = seq[off + i];
            const uint8_t code = base_code(b);
            uint8_t q = qual.empty() ? 30 : phred_from_char(qual[off + i]);
            double e = (!calib.empty() && q < calib.size()) ? calib[q] : phred_error(q);
            if (code != kInvalidBase && ((m >> code) & 1u)) {
                s += match_bits(e, iupac_size(m));
            } else {
                // An N in the read is not evidence against the tag; it is no evidence at all.
                // Scoring it as a mismatch at its (low) quality already handles that, and an
                // uncallable base carries Q2, worth -0.6 bits.
                s += weight_[i] * mismatch_bits(e);
            }
            // Once this offset cannot beat the incumbent even with every remaining position
            // matching perfectly, stop. Cheap and it cuts the common case (no match) short.
            if (s + 2.0 * static_cast<double>(plen - i - 1) < best) break;
        }
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
    for (size_t i = 0; i < plen; ++i) {
        if (!is_umi_[i]) continue;
        out.umi.push_back(seq[static_cast<size_t>(best_off) + i]);
        out.umi_qual.push_back(qual.empty() ? char_from_phred(30)
                                            : qual[static_cast<size_t>(best_off) + i]);
    }
    return out;
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
    if (per_pattern.min_score < 0.0 && !patterns_.empty()) {
        per_pattern.min_score = patterns_[0].default_min_score(seq.size(), patterns_.size());
    }

    double best = -1e300, second = -1e300;
    for (size_t i = 0; i < patterns_.size(); ++i) {
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
