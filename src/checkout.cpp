#include "migec/checkout.hpp"

#include <algorithm>

#include "migec/types.hpp"

namespace migec {

Checkout::Checkout(const PatternSet& patterns, CheckoutParams params)
    : patterns_(patterns), params_(std::move(params)) {
    counters_.per_sample.assign(patterns.size(), 0);
}

std::string Checkout::header_tags(const std::string& umi, const std::string& umi_qual,
                                  const std::string& sample) {
    // RX/QX are the SAM standard tags for a UMI and its qualities (fgbio, Picard and umi_tools all
    // read RX). BC is the sample barcode. Tabs between tags, one space before the first -- that is
    // what makes the comment survive `bwa mem -C` into a valid SAM record.
    std::string out;
    out.reserve(umi.size() * 2 + sample.size() + 24);
    if (!umi.empty()) {
        out += "RX:Z:";
        out += umi;
        if (!umi_qual.empty()) {
            out += "\tQX:Z:";
            out += umi_qual;
        }
    }
    if (!sample.empty()) {
        if (!out.empty()) out += "\t";
        out += "BC:Z:";
        out += sample;
    }
    return out;
}

CheckoutRead Checkout::process(std::string_view seq, std::string_view qual) {
    ++counters_.total;
    CheckoutRead out;

    PatternSet::Assignment a = patterns_.assign(seq, qual, params_.match);
    if (a.ambiguous) {
        ++counters_.ambiguous;
        return out;
    }
    if (a.sample < 0) {
        ++counters_.unmatched;
        return out;
    }

    const PatternMatch& m = a.match;
    if (params_.reject_umi_n &&
        m.umi.find_first_not_of("ACGTacgt") != std::string::npos) {
        ++counters_.bad_umi;
        return out;
    }
    if (params_.min_umi_quality > 0 && !m.umi_qual.empty()) {
        uint8_t worst = kMaxPhred;
        for (char c : m.umi_qual) worst = std::min(worst, phred_from_char(c));
        if (worst < params_.min_umi_quality) {
            ++counters_.bad_umi;
            return out;
        }
    }

    size_t begin = 0, end = seq.size();
    switch (params_.trim) {
        case TrimMode::kNone:
            break;
        case TrimMode::kPattern:
            // Everything 5' of the payload goes: adapter, sample tag, UMI. This is synthetic
            // sequence and leaving it in costs soft-clips at best and mismapping at worst.
            begin = static_cast<size_t>(m.payload_begin);
            break;
        case TrimMode::kPatternOnly:
            // Splice the pattern out, keeping the flank before it. Returning a view is impossible
            // here, so we keep the 3' side -- the flank is available to the caller via the match
            // offset if it is genuinely wanted.
            begin = static_cast<size_t>(m.payload_begin);
            break;
    }
    if (begin > end) begin = end;

    if (static_cast<int>(end - begin) < params_.min_payload) {
        ++counters_.short_payload;
        return out;
    }

    out.ok = true;
    out.sample = a.sample;
    out.umi = m.umi;
    out.umi_qual = m.umi_qual;
    out.seq = seq.substr(begin, end - begin);
    out.qual = qual.empty() ? std::string_view() : qual.substr(begin, end - begin);
    out.score = m.score;
    ++counters_.assigned;
    ++counters_.per_sample[static_cast<size_t>(a.sample)];
    return out;
}

}  // namespace migec
