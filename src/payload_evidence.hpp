// Shared first-draft observation semantics for correction, clonality and residual diagnostics.
#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include "migec/umi_stats.hpp"
#include "migec/types.hpp"
namespace migec::detail {
struct PayloadEmission {
    double delta, weight, same_mismatch, same_match, related_mismatch, related_match;
};

inline PayloadEmission payload_emission(int q1, int q2, const CorrectionParams& params) {
    double delta = params.payload_error;
    const double e1 = q1 >= 0 ? std::min(0.75, phred_error(static_cast<uint8_t>(q1))) : -1.0;
    const double e2 = q2 >= 0 ? std::min(0.75, phred_error(static_cast<uint8_t>(q2))) : -1.0;
    if (e1 == 0.75 || e2 == 0.75) {
        delta = 0.75;  // preserve exact zero information, including mixed/missing quality
    } else if (e1 >= 0.0 && e2 >= 0.0) {
        delta = std::max(delta, e1 + e2 - 4.0 * e1 * e2 / 3.0);
    }
    delta = std::clamp(delta, 1e-12, 0.75);
    const double weight = 1.0 - 4.0 * delta / 3.0;
    const double related = std::clamp(delta + params.payload_same_fraction * weight,
                                     1e-12, 1.0 - 1e-12);
    return {delta, weight, std::log(delta), std::log1p(-delta),
            std::log(related), std::log1p(-related)};
}

// Immutable model-parameter emissions, not cached query results. Profiling found
// per-candidate logarithms dominated scoring; default Phred pairs have a tiny table.
inline const auto& default_payload_emissions() {
    static const auto table = [] {
        constexpr size_t width = kMaxPhred + 2;
        std::array<PayloadEmission, width * width> values;
        const CorrectionParams params;
        for (size_t a = 0; a < width; ++a)
            for (size_t b = 0; b < width; ++b)
                values[a * width + b] = payload_emission(
                    a == width - 1 ? -1 : static_cast<int>(a),
                    b == width - 1 ? -1 : static_cast<int>(b), params);
        return values;
    }();
    return table;
}

struct PayloadComparison {
    int comparable = 0, mismatches = 0, raw_comparable = 0, raw_mismatches = 0;
    double excess = 0.0, information = 0.0;
    double log_same = 0.0, log_related = 0.0, log_random = 0.0;
    double latent_fraction() const {
        return information > 0.0 ? std::clamp(excess / information, 0.0, 1.0) : 0.0;
    }
};

// Payloads are first-read observations. The pair-disagreement floor is the
// existing payload_error, not two copies of that floor. Unknown Phred bytes use
// the documented historical floor, rather than being interpreted as Q93.
inline PayloadComparison compare_payload(const std::vector<char>& data, int width, bool quality,
                                  size_t a, size_t b, const CorrectionParams& params, int comparison_width = 0, bool quality_emissions = true) {
    PayloadComparison result;
    const size_t stride = static_cast<size_t>(width) * (quality ? 2 : 1);
    const int compared = comparison_width > 0 ? std::min(width, comparison_width) : width;
    if (!quality || !quality_emissions || params.cell_length > 0) {
        // Historical constant-emission path: aggregate counts before logarithms.
        const auto emission = payload_emission(-1, -1, params);
        for (int position = 0; position < compared; ++position) {
            const char left = data[a * stride + static_cast<size_t>(position)];
            const char right = data[b * stride + static_cast<size_t>(position)];
            if (!left || !right || left == 'N' || right == 'N') continue;
            ++result.comparable;
            result.mismatches += left != right;
        }
        const int matches = result.comparable - result.mismatches;
        result.raw_comparable = result.comparable;
        result.raw_mismatches = result.mismatches;
        result.excess = result.mismatches;
        result.information = result.comparable;
        result.log_same = result.mismatches * emission.same_mismatch + matches * emission.same_match;
        result.log_related = result.mismatches * emission.related_mismatch + matches * emission.related_match;
        result.log_random = result.mismatches * std::log(0.75) + matches * std::log(0.25);
        return result;
    }
    for (int position = 0; position < compared; ++position) {
        const size_t x = a * stride + static_cast<size_t>(position);
        const size_t y = b * stride + static_cast<size_t>(position);
        const char left = data[x], right = data[y];
        if (left == 0 || right == 0 || left == 'N' || right == 'N') continue;
        ++result.raw_comparable;
        result.raw_mismatches += left != right;
        constexpr size_t table_width = kMaxPhred + 2;
        size_t q1 = table_width - 1, q2 = table_width - 1;
        if (quality && params.cell_length == 0) {
            const unsigned char first = data[x + static_cast<size_t>(width)];
            const unsigned char second = data[y + static_cast<size_t>(width)];
            if (first >= 33 && first <= 126)
                q1 = std::min<size_t>(kMaxPhred, first - 33);
            if (second >= 33 && second <= 126)
                q2 = std::min<size_t>(kMaxPhred, second - 33);
        }
        const auto emission = params.payload_error == 2e-3 && params.payload_same_fraction == 0.05
            ? default_payload_emissions()[q1 * table_width + q2]
            : payload_emission(q1 == table_width - 1 ? -1 : static_cast<int>(q1),
                               q2 == table_width - 1 ? -1 : static_cast<int>(q2), params);
        if (emission.weight <= 0.0) continue;
        const bool mismatch = left != right;
        ++result.comparable;
        result.mismatches += mismatch;
        result.excess += static_cast<double>(mismatch) - emission.delta;
        result.information += emission.weight;
        result.log_same += mismatch ? emission.same_mismatch : emission.same_match;
        result.log_related += mismatch ? emission.related_mismatch : emission.related_match;
        result.log_random += mismatch ? std::log(0.75) : std::log(0.25);
    }
    return result;
}

}  // namespace migec::detail
