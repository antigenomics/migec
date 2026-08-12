#include "migec/umi_stats.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

#include "migec/types.hpp"

namespace migec {
namespace {

int log2_bin(uint32_t count) {
    int b = 0;
    while (count > 1 && b < CoverageHistogram::kBins - 1) {
        count >>= 1;
        ++b;
    }
    return b;
}

// Zero-truncated Poisson: P(X = k | X >= 1), the right likelihood for an error child, because a
// child with zero reads is not observed at all and must not carry probability mass.
double zt_poisson(uint32_t k, double lambda) {
    if (k == 0) return 0.0;
    if (lambda <= 0.0) return 0.0;
    double logp = -lambda + k * std::log(lambda);
    for (uint32_t i = 2; i <= k; ++i) logp -= std::log(static_cast<double>(i));
    const double denom = 1.0 - std::exp(-lambda);
    if (denom <= 0.0) return 0.0;
    return std::exp(logp) / denom;
}

// Probability mass that a polymerase error child of a parent with `c_par` reads is seen with
// exactly `c_child` reads.
//
// Under a branching process the child's share f of the family has density ~ 1/f^2 (Luria-Delbruck):
// an error entering at cycle k reaches a fraction ~ (1+e)^-k of the descendants, and that is a
// very heavy tail compared with a Poisson on the sequencing rate. This is the component a
// sequencing-only model misses, and it is why an error child can be several percent of its parent
// -- MIGEC merged below 10% and MAGERI below 1/20, both far above anything eps/3 predicts.
//
// Normalised over f in [1/(c_par+1), f_max], then transformed from f to a count via
// |df/dc| = c_par / (c_child + c_par)^2.
double ld_pmf(uint32_t c_child, uint32_t c_par, double f_max) {
    const double C = static_cast<double>(c_par);
    const double c = static_cast<double>(c_child);
    const double f = c / (c + C);
    const double f_min = 1.0 / (C + 1.0);
    if (f <= 0.0 || f > f_max || f_max <= f_min) return 0.0;
    const double norm = 1.0 / f_min - 1.0 / f_max;  // integral of 1/f^2 over [f_min, f_max]
    if (norm <= 0.0) return 0.0;
    const double pdf_f = (1.0 / (f * f)) / norm;
    const double dfdc = C / ((c + C) * (c + C));
    return pdf_f * dfdc;
}

}  // namespace

uint64_t CoverageHistogram::total_reads() const {
    return std::accumulate(reads.begin(), reads.end(), uint64_t{0});
}
uint64_t CoverageHistogram::total_units() const {
    return std::accumulate(units.begin(), units.end(), uint64_t{0});
}
double CoverageHistogram::mean_reads_per_umi() const {
    const uint64_t u = total_units();
    return u ? static_cast<double>(total_reads()) / static_cast<double>(u) : 0.0;
}
double CoverageHistogram::reads_in_migs_at_least(uint32_t min_size) const {
    const uint64_t tot = total_reads();
    if (!tot) return 0.0;
    uint64_t kept = 0;
    for (int b = 0; b < kBins; ++b) {
        if ((1u << b) >= min_size) kept += reads[static_cast<size_t>(b)];
    }
    return static_cast<double>(kept) / static_cast<double>(tot);
}

double UmiComposition::entropy(int j) const {
    double h = 0.0;
    for (double p : freq[static_cast<size_t>(j)]) {
        if (p > 0.0) h -= p * std::log2(p);
    }
    return h;
}
double UmiComposition::information(int j) const { return 2.0 - entropy(j); }
double UmiComposition::total_entropy() const {
    double h = 0.0;
    for (int j = 0; j < length; ++j) h += entropy(j);
    return h;
}
double UmiComposition::total_information() const { return 2.0 * length - total_entropy(); }
double UmiComposition::collision(int j) const {
    double m = 0.0;
    for (double p : freq[static_cast<size_t>(j)]) m += p * p;
    return m;
}
double UmiComposition::effective_length() const {
    double l = 0.0;
    for (int j = 0; j < length; ++j) l -= std::log(collision(j)) / std::log(4.0);
    return l;
}
double UmiComposition::effective_space() const {
    double prod = 1.0;
    for (int j = 0; j < length; ++j) prod *= collision(j);
    return prod > 0.0 ? 1.0 / prod : 0.0;
}
double UmiComposition::expected_collisions(double n_molecules) const {
    double prod = 1.0;
    for (int j = 0; j < length; ++j) prod *= collision(j);
    return 0.5 * n_molecules * n_molecules * prod;
}

uint64_t UmiCounts::total() const {
    uint64_t t = 0;
    for (const auto& kv : counts_) t += kv.second;
    return t;
}

CoverageHistogram UmiCounts::histogram() const {
    CoverageHistogram h;
    for (const auto& kv : counts_) {
        const int b = log2_bin(kv.second);
        h.reads[static_cast<size_t>(b)] += kv.second;
        h.units[static_cast<size_t>(b)] += 1;
    }
    return h;
}

UmiComposition UmiCounts::composition(bool weight_by_reads) const {
    UmiComposition c;
    c.length = length_;
    c.freq.assign(static_cast<size_t>(length_), {0.0, 0.0, 0.0, 0.0});
    double total = 0.0;
    for (const auto& kv : counts_) {
        const double w = weight_by_reads ? static_cast<double>(kv.second) : 1.0;
        for (int j = 0; j < length_; ++j) {
            const uint8_t code = static_cast<uint8_t>((kv.first >> (62 - 2 * j)) & 3u);
            c.freq[static_cast<size_t>(j)][code] += w;
        }
        total += w;
    }
    if (total > 0.0) {
        for (auto& row : c.freq) {
            for (double& v : row) v /= total;
        }
    }
    return c;
}

double estimate_umi_error(const UmiCounts& counts, const UmiComposition& comp) {
    const int L = counts.length();
    if (L <= 0 || counts.distinct() < 2) return 0.0;

    // Observed distinct-barcode pairs at Hamming distance 1, counted once each.
    uint64_t d1_obs = 0;
    const auto& m = counts.map();
    for (const auto& kv : m) {
        for (int j = 0; j < L; ++j) {
            const int shift = 62 - 2 * j;
            const uint64_t cur = (kv.first >> shift) & 3u;
            for (uint64_t b = 0; b < 4; ++b) {
                if (b == cur) continue;
                const uint64_t nb = (kv.first & ~(uint64_t{3} << shift)) | (b << shift);
                if (nb > kv.first && m.count(nb)) ++d1_obs;  // count each pair once
            }
        }
    }

    const double n = static_cast<double>(counts.distinct());
    double p_coll = 1.0;
    for (int j = 0; j < L; ++j) p_coll *= comp.collision(j);
    // Independent pairs that happen to sit at distance 1. Per position there are 3 alternatives,
    // so the distance-1 shell has probability ~ P_coll * 3L / (the per-position collision), which
    // for a near-uniform composition is 3L/4^L; use that form, guarded for the uniform case.
    const double p_d1 = p_coll * 3.0 * L;
    const double d1_ind = 0.5 * n * (n - 1.0) * p_d1;

    const double excess = static_cast<double>(d1_obs) - d1_ind;
    if (excess <= 0.0) return 0.0;

    // Bisect on log(eps) against the parent-child plus sibling expectation.
    auto expected = [&](double eps) {
        double parent_child = 0.0, sibling = 0.0;
        for (const auto& kv : m) {
            const double t = 1.0 - std::exp(-static_cast<double>(kv.second) * eps);
            parent_child += t;
            sibling += t * t;
        }
        return 3.0 * L * (parent_child + sibling);
    };

    double lo = 1e-8, hi = 0.2;
    if (expected(hi) < excess) return hi;
    if (expected(lo) > excess) return lo;
    for (int it = 0; it < 60; ++it) {
        const double mid = std::sqrt(lo * hi);
        if (expected(mid) < excess) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return std::sqrt(lo * hi);
}

CorrectionResult correct_umis(const UmiCounts& counts, const CorrectionParams& params) {
    CorrectionResult res;
    const int L = counts.length();
    const auto& m = counts.map();
    res.corrected.reserve(m.size());
    for (const auto& kv : m) res.corrected[kv.first] = kv.second;
    if (L <= 0 || m.size() < 2) {
        res.molecules_observed = m.size();
        res.molecules_corrected = static_cast<double>(m.size());
        return res;
    }

    const UmiComposition comp = counts.composition(false);
    double eps = params.sequencing_error;
    if (eps < 0.0) eps = estimate_umi_error(counts, comp);
    if (eps <= 0.0) eps = 1e-4;  // a floor, so correction still runs on a clean small library
    res.estimated_error = eps;

    double p_coll = 1.0;
    for (int j = 0; j < L; ++j) p_coll *= comp.collision(j);
    const double n = static_cast<double>(m.size());
    const double space = comp.effective_space();
    res.saturated = space > 0.0 && n > 0.05 * space;

    // Prior that a substitution at one specific position-and-base is polymerase-derived. Per base
    // per cycle it is eps_pol, over the cycles that matter, spread over the 3 alternative bases.
    const double rho_pol = std::min(
        0.9, params.polymerase_error * std::max(1, params.pcr_cycles) / 3.0 * 3.0 * L);

    // The independent hypothesis: some *other real molecule* happens to occupy this exact
    // neighbouring barcode. Its probability is (number of molecules) x (probability a molecule
    // draws that specific barcode) -- and p_coll is exactly that probability, since sum_u p_u^2 is
    // the chance two independent draws coincide.
    const double a_ind = n * p_coll;

    // ...and if it is a real molecule, its read count follows the library's own MIG size
    // distribution. Using the empirical distribution rather than a parametric one means the test
    // adapts to how deeply the library was sequenced without another tunable.
    std::unordered_map<uint32_t, double> size_pmf;
    for (const auto& kv : m) size_pmf[kv.second] += 1.0;
    for (auto& kv : size_pmf) kv.second /= n;
    const double size_floor = 1.0 / (n + 1.0);  // never claim a size is impossible

    // Order by count descending so a parent is always processed before its children, then walk
    // each barcode's 3L neighbourhood.
    std::vector<std::pair<uint64_t, uint32_t>> order(m.begin(), m.end());
    std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) {
        if (a.second != b.second) return a.second > b.second;
        return a.first < b.first;  // total order, so the result is reproducible
    });

    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const uint64_t child = it->first;
        const uint32_t c_child = it->second;
        uint64_t best_parent = 0;
        double best_post = 0.0;

        for (int j = 0; j < L; ++j) {
            const int shift = 62 - 2 * j;
            const uint64_t cur = (child >> shift) & 3u;
            for (uint64_t b = 0; b < 4; ++b) {
                if (b == cur) continue;
                const uint64_t cand = (child & ~(uint64_t{3} << shift)) | (b << shift);
                auto f = m.find(cand);
                if (f == m.end()) continue;
                const uint32_t c_par = f->second;
                if (c_par <= c_child) continue;  // a parent must be strictly larger
                if (static_cast<double>(c_child) >
                    params.max_child_fraction * static_cast<double>(c_par)) {
                    continue;
                }

                // Two ways to be an error child. Sequencing miscalls land on one specific
                // alternative base, so the rate per neighbour is eps/3, not eps.
                const double lam = static_cast<double>(c_par) * eps / 3.0;
                const double l_seq = zt_poisson(c_child, lam);
                const double l_pol = ld_pmf(c_child, c_par, params.max_child_fraction);
                const double l_err = (1.0 - rho_pol) * l_seq + rho_pol * l_pol;

                // ...against being a real molecule that happens to sit one substitution away and
                // to have this many reads.
                auto sp = size_pmf.find(c_child);
                const double p_size = sp == size_pmf.end() ? size_floor
                                                           : std::max(sp->second, size_floor);
                const double l_ind = std::max(a_ind * p_size, 1e-300);

                const double post = l_err / (l_err + l_ind);
                if (post > best_post) {
                    best_post = post;
                    best_parent = cand;
                }
            }
        }

        if (best_post >= params.min_posterior && best_parent != 0) {
            // Resolve the chain: the parent may itself already have been merged.
            uint64_t root = best_parent;
            for (int guard = 0; guard < 64; ++guard) {
                auto p = res.parent.find(root);
                if (p == res.parent.end()) break;
                root = p->second;
            }
            if (root == child) continue;  // never make a cycle
            res.parent[child] = root;
            res.corrected[root] += res.corrected[child];
            res.merged_reads += res.corrected[child];
            res.corrected.erase(child);
            ++res.merged;
        }
    }

    res.molecules_observed = res.corrected.size();
    // Two molecules drawing the same UMI are invisible to any method, so the observed count is
    // biased low. Inverting the Poisson occupancy recovers the estimate:
    //     M_hat = S_eff * -ln(1 - M_obs / S_eff)
    //
    // ...but only while the space is not nearly full. S_eff is estimated from the composition of
    // the observed barcodes, so as occupancy approaches 1 the estimate collapses onto M_obs and
    // the correction becomes circular -- it would report "no collisions" for the most collided
    // library possible. Past 90% occupancy we decline to estimate and say so via `saturated`.
    const double m_obs = static_cast<double>(res.molecules_observed);
    if (space > 0.0 && m_obs < 0.9 * space) {
        res.molecules_corrected = space * -std::log1p(-m_obs / space);
    } else {
        res.molecules_corrected = m_obs;  // not estimable; see CorrectionResult::saturated
        res.saturated = true;
    }
    return res;
}

}  // namespace migec
