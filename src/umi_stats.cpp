#include "migec/umi_stats.hpp"

#include <algorithm>
#include <cmath>
#include <array>
#include <numeric>
#include <unordered_map>
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

void UmiCounts::flush() const {
    if (buf_.empty()) return;
    std::sort(buf_.begin(), buf_.end(),
              [](const Entry& a, const Entry& b) { return a.key < b.key; });
    // Run-length reduce the buffer in place.
    size_t w = 0;
    for (size_t r = 0; r < buf_.size(); ++r) {
        if (w > 0 && buf_[w - 1].key == buf_[r].key) {
            buf_[w - 1].count += buf_[r].count;
        } else {
            buf_[w++] = buf_[r];
        }
    }
    buf_.resize(w);

    auto set_next_flush = [this] {
        flush_at_ = std::min(buffer_limit_, std::max(kMinBuffer, entries_.size() / 2));
    };

    if (entries_.empty()) {
        entries_.swap(buf_);
        buf_.clear();
        set_next_flush();
        return;
    }

    // Merge two sorted runs *backwards into the grown array* rather than into a fresh vector: at
    // this size the transient copy is the peak memory of the whole process.
    const size_t n = entries_.size();
    // reserve() before resize(): resize alone grows geometrically, so the array would sit at up to
    // twice the bytes it needs for the whole run. This is the largest allocation in the process.
    entries_.reserve(n + w);
    entries_.resize(n + w);
    size_t i = n, j = w, out = n + w;
    while (i > 0 && j > 0) {
        const Entry& a = entries_[i - 1];
        const Entry& b = buf_[j - 1];
        if (a.key == b.key) {
            entries_[--out] = Entry{a.key, a.count + b.count};
            --i;
            --j;
        } else if (a.key > b.key) {
            entries_[--out] = a;
            --i;
        } else {
            entries_[--out] = b;
            --j;
        }
    }
    while (j > 0) entries_[--out] = buf_[--j];
    // Equal keys collapsed, so the merged run can be shorter than the space reserved for it; the
    // survivors sit at the top and are shifted down.
    if (out > 0) {
        std::move(entries_.begin() + static_cast<long>(out), entries_.end(),
                  entries_.begin() + static_cast<long>(i));
        entries_.resize(i + (n + w - out));
    }
    buf_.clear();
    set_next_flush();
}

void UmiCounts::merge(const UmiCounts& other) {
    other.flush();
    for (const Entry& e : other.entries_) add(e.key, e.count);
}

const uint32_t* UmiCounts::find(uint64_t key) const {
    flush();
    auto it = std::lower_bound(entries_.begin(), entries_.end(), key,
                               [](const Entry& e, uint64_t k) { return e.key < k; });
    if (it == entries_.end() || it->key != key) return nullptr;
    return &it->count;
}

size_t UmiCounts::memory_bytes() const {
    return entries_.capacity() * sizeof(Entry) + buf_.capacity() * sizeof(Entry);
}

size_t index_of(const UmiCounts& counts, uint64_t key) {
    const std::vector<UmiCounts::Entry>& e = counts.entries();
    auto it = std::lower_bound(e.begin(), e.end(), key,
                               [](const UmiCounts::Entry& x, uint64_t k) { return x.key < k; });
    if (it == e.end() || it->key != key) return static_cast<size_t>(-1);
    return static_cast<size_t>(it - e.begin());
}

CoverageHistogram UmiCounts::histogram() const {
    CoverageHistogram h;
    for (const Entry& e : entries()) {
        const int b = log2_bin(e.count);
        h.reads[static_cast<size_t>(b)] += e.count;
        h.units[static_cast<size_t>(b)] += 1;
    }
    return h;
}

UmiComposition UmiCounts::composition(bool weight_by_reads) const {
    UmiComposition c;
    c.length = length_;
    c.freq.assign(static_cast<size_t>(length_), {0.0, 0.0, 0.0, 0.0});
    double total = 0.0;
    for (const Entry& e : entries()) {
        const double w = weight_by_reads ? static_cast<double>(e.count) : 1.0;
        for (int j = 0; j < length_; ++j) {
            const uint8_t code = static_cast<uint8_t>((e.key >> (62 - 2 * j)) & 3u);
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

    const std::vector<UmiCounts::Entry>& m = counts.entries();
    // Binary search over the entry array we already hold. Not UmiCounts::find, which flushes the
    // append buffer and could reallocate the very array `m` refers to.
    auto present = [&m](uint64_t key) {
        auto it = std::lower_bound(m.begin(), m.end(), key,
                                   [](const UmiCounts::Entry& e, uint64_t k) { return e.key < k; });
        return it != m.end() && it->key == key;
    };

    // Observed distinct-barcode pairs at Hamming distance 1, counted once each.
    uint64_t d1_obs = 0;
    for (const UmiCounts::Entry& e : m) {
        for (int j = 0; j < L; ++j) {
            const int shift = 62 - 2 * j;
            const uint64_t cur = (e.key >> shift) & 3u;
            for (uint64_t b = 0; b < 4; ++b) {
                if (b == cur) continue;
                const uint64_t nb = (e.key & ~(uint64_t{3} << shift)) | (b << shift);
                if (nb > e.key && present(nb)) ++d1_obs;  // count each pair once
            }
        }
    }

    const double n = static_cast<double>(counts.distinct());
    double p_coll = 1.0;
    for (int j = 0; j < L; ++j) p_coll *= comp.collision(j);
    // Independent pairs that happen to sit at distance 1: agree everywhere but position j, and
    // differ there. Summing over j,
    //     P_d1 = sum_j (prod_{k != j} m_k) * (1 - m_j) = P_coll * sum_j (1 - m_j)/m_j
    // which is 3L * P_coll only for a uniform composition. Since m_j > 1/4 whenever the
    // composition is skewed, the uniform form overstates the independent term, understates the
    // excess, and so *underestimates* the error rate -- the direction that leaves errors
    // uncorrected.
    double shell = 0.0;
    for (int j = 0; j < L; ++j) {
        const double mj = comp.collision(j);
        if (mj > 0.0) shell += (1.0 - mj) / mj;
    }
    const double d1_ind = 0.5 * n * (n - 1.0) * p_coll * shell;

    const double excess = static_cast<double>(d1_obs) - d1_ind;
    if (excess <= 0.0) return 0.0;

    // Bisect on log(eps) against the parent-child plus sibling expectation.
    //
    // For one parent with c reads and one specific neighbour (position j, base b), the chance some
    // read carries *that* error is 1 - (1 - eps/3)^c ~ 1 - exp(-c eps/3): eps/3, not eps, because a
    // miscall has to land on that one alternative base out of three. There are 3L such neighbours,
    // hence the 3L factor outside. Using eps in the exponent makes the expectation 3x too large at
    // small c*eps and so returns an eps 3x too small -- which it did, uniformly, at every
    // occupancy from 0.3% upwards.
    auto expected = [&](double eps) {
        double parent_child = 0.0, sibling = 0.0;
        for (const UmiCounts::Entry& e : m) {
            const double t = 1.0 - std::exp(-static_cast<double>(e.count) * eps / 3.0);
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

BarcodeSpace barcode_space(const UmiComposition& comp, uint64_t observed_barcodes,
                           double saturation) {
    BarcodeSpace b;
    b.length = comp.length;
    b.nominal_space = std::pow(4.0, static_cast<double>(comp.length));
    b.effective_space = comp.effective_space();
    b.effective_length = comp.effective_length();
    b.bias_loss = b.nominal_space > 0.0 ? 1.0 - b.effective_space / b.nominal_space : 0.0;
    b.observed = observed_barcodes;
    const double n = static_cast<double>(observed_barcodes);
    const double S = b.effective_space;
    if (S <= 0.0 || n <= 0.0) return b;

    b.occupancy = n / S;
    b.saturated = b.occupancy >= saturation;
    if (b.saturated) {
        // S is inferred from the observed barcodes, so at saturation the inversion collapses onto
        // the observed count and would report "no collisions" for the most collided library there
        // can be. Decline rather than mislead; the fields stay at their observed values.
        b.molecules = n;
        b.lambda = 0.0;
        return b;
    }
    b.lambda = -std::log1p(-b.occupancy);
    b.molecules = S * b.lambda;
    b.hidden = b.molecules - n;
    // P(k > 1 | k >= 1) for k ~ Poisson(lambda).
    const double e = std::exp(-b.lambda);
    const double occupied = 1.0 - e;
    b.p_multi = occupied > 0.0 ? (occupied - b.lambda * e) / occupied : 0.0;
    return b;
}

ErrorBudget error_budget(const UmiComposition& comp, const std::array<uint64_t, 61>& phred_counts,
                         double estimated, uint64_t observed_barcodes, double polymerase_error,
                         int pcr_cycles) {
    ErrorBudget b;
    uint64_t total = 0;
    double sum_e = 0.0, sum_q = 0.0;
    for (int q = 0; q <= kMaxPhred; ++q) {
        const uint64_t n = phred_counts[static_cast<size_t>(q)];
        if (!n) continue;
        total += n;
        sum_e += static_cast<double>(n) * phred_error(static_cast<uint8_t>(q));
        sum_q += static_cast<double>(n) * q;
    }
    if (total) {
        // The mean of 10^(-Q/10), not 10^(-mean Q/10). Averaging Q first hides the low-Q tail,
        // which is where nearly all of the error is.
        b.from_phred = sum_e / static_cast<double>(total);
        b.mean_phred = sum_q / static_cast<double>(total);
    }
    b.from_polymerase = polymerase_error * std::max(1, pcr_cycles);
    b.predicted = b.from_phred + b.from_polymerase;
    b.estimated = estimated;
    b.ratio = b.predicted > 0.0 ? estimated / b.predicted : 0.0;
    if (comp.length > 0 && b.predicted > 0.0 && b.predicted < 1.0) {
        b.barcodes_with_error = 1.0 - std::pow(1.0 - b.predicted, comp.length);
    }
    // The distance-1 shell around a barcode holds 3L neighbours. If a large share of those are
    // themselves real barcodes, the observed pair count is dominated by coincidence and the excess
    // the estimator reads is a small difference of two large numbers.
    const double S = comp.effective_space();
    if (S > 0.0 && comp.length > 0) {
        const double occ = static_cast<double>(observed_barcodes) / S;
        b.neighbour_occupancy = occ > 1.0 ? 1.0 : occ;
        b.estimate_unreliable = b.neighbour_occupancy > 0.05;
    }
    return b;
}

CorrectionResult correct_umis(const UmiCounts& counts, const CorrectionParams& params) {
    CorrectionResult res;
    const int L = counts.length();
    const std::vector<UmiCounts::Entry>& m = counts.entries();
    const size_t n_entries = m.size();

    res.root.resize(n_entries);
    res.corrected.resize(n_entries);
    for (size_t i = 0; i < n_entries; ++i) {
        res.root[i] = static_cast<uint32_t>(i);
        res.corrected[i] = m[i].count;
    }
    if (L <= 0 || n_entries < 2) {
        res.molecules_observed = n_entries;
        res.molecules_corrected = static_cast<double>(n_entries);
        return res;
    }

    // Index of a packed barcode in the sorted entry array. Binary search: the alternative would be
    // a side hash map, which is exactly the allocation this class exists to avoid.
    auto find_idx = [&m](uint64_t key) -> size_t {
        auto it = std::lower_bound(m.begin(), m.end(), key,
                                   [](const UmiCounts::Entry& e, uint64_t k) { return e.key < k; });
        if (it == m.end() || it->key != key) return static_cast<size_t>(-1);
        return static_cast<size_t>(it - m.begin());
    };

    const UmiComposition comp = counts.composition(false);
    double eps = params.sequencing_error;
    if (eps < 0.0) eps = estimate_umi_error(counts, comp);
    if (eps <= 0.0) eps = 1e-4;  // a floor, so correction still runs on a clean small library
    res.estimated_error = eps;

    double p_coll = 1.0;
    for (int j = 0; j < L; ++j) p_coll *= comp.collision(j);
    const double n = static_cast<double>(n_entries);
    const double space = comp.effective_space();
    res.saturated = space > 0.0 && n > 0.05 * space;

    // Prior that a neighbour one substitution away is polymerase-derived rather than a miscall:
    // eps_pol per base per cycle, over the cycles that matter, over the L barcode positions.
    const double rho_pol =
        std::min(0.9, params.polymerase_error * std::max(1, params.pcr_cycles) * L);

    // The independent hypothesis: some *other real molecule* happens to occupy this exact
    // neighbouring barcode. Its probability is (number of molecules) x (probability a molecule
    // draws that specific barcode) -- and p_coll is exactly that probability, since sum_u p_u^2 is
    // the chance two independent draws coincide.
    const double a_ind = n * p_coll;

    // ...and if it is a real molecule, its read count follows the library's own MIG size
    // distribution. Using the empirical distribution rather than a parametric one means the test
    // adapts to how deeply the library was sequenced without another tunable.
    std::unordered_map<uint32_t, double> size_pmf;  // keyed by MIG size, so it stays small
    for (const UmiCounts::Entry& e : m) size_pmf[e.count] += 1.0;
    for (auto& kv : size_pmf) kv.second /= n;
    const double size_floor = 1.0 / (n + 1.0);  // never claim a size is impossible

    // Order by count descending, then walk each barcode's 3L neighbourhood from the smallest MIG
    // upwards. Indices, not copies of the entries: 4 bytes each rather than 16.
    std::vector<uint32_t> order(n_entries);
    for (size_t i = 0; i < n_entries; ++i) order[i] = static_cast<uint32_t>(i);
    std::sort(order.begin(), order.end(), [&m](uint32_t a, uint32_t b) {
        if (m[a].count != m[b].count) return m[a].count > m[b].count;
        return m[a].key < m[b].key;  // total order, so the result is reproducible
    });

    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const size_t child_idx = *it;
        const uint64_t child = m[child_idx].key;
        const uint32_t c_child = m[child_idx].count;
        size_t best_parent = static_cast<size_t>(-1);
        double best_post = 0.0;

        for (int j = 0; j < L; ++j) {
            const int shift = 62 - 2 * j;
            const uint64_t cur = (child >> shift) & 3u;
            for (uint64_t b = 0; b < 4; ++b) {
                if (b == cur) continue;
                const uint64_t cand = (child & ~(uint64_t{3} << shift)) | (b << shift);
                const size_t cand_idx = find_idx(cand);
                if (cand_idx == static_cast<size_t>(-1)) continue;
                const uint32_t c_par = m[cand_idx].count;
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
                    best_parent = cand_idx;
                }
            }
        }

        if (best_post >= params.min_posterior && best_parent != static_cast<size_t>(-1)) {
            // Follow the parent to its current root -- it may itself already have been merged.
            uint32_t root = static_cast<uint32_t>(best_parent);
            for (int guard = 0; guard < 64 && res.root[root] != root; ++guard) root = res.root[root];
            if (root == child_idx) continue;  // never make a cycle
            res.root[child_idx] = root;
            res.corrected[root] += res.corrected[child_idx];
            res.merged_reads += res.corrected[child_idx];
            res.corrected[child_idx] = 0;
            ++res.merged;
        }
    }

    // Flatten. A child merged early can point at a parent that was itself merged later in the
    // walk, so the invariant "root[root[i]] == root[i]" only holds after this pass. The read
    // counts were already correct -- they follow the chain as it forms -- but a consumer reading
    // `root` directly would otherwise land on an intermediate.
    for (size_t i = 0; i < n_entries; ++i) {
        uint32_t r = res.root[i];
        for (int guard = 0; guard < 64 && res.root[r] != r; ++guard) r = res.root[r];
        res.root[i] = r;
    }

    res.molecules_observed = 0;
    for (uint32_t c : res.corrected) {
        if (c > 0) ++res.molecules_observed;
    }
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
