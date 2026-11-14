#include "migec/umi_stats.hpp"

#include <algorithm>
#include <cmath>
#include <array>
#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include "migec/parallel.hpp"
#include "migec/types.hpp"
#include "barcode_search.hpp"

namespace migec {
namespace {
using detail::barcode_neighbors;
using detail::BarcodeSearch;

int log2_bin(uint32_t count) {
    int b = 0;
    while (count > 1 && b < CoverageHistogram::kBins - 1) {
        count >>= 1;
        ++b;
    }
    return b;
}

// Poisson pmf, UNtruncated, as an expected count rather than a conditional probability.
//
// The zero-truncated form -- P(X = k | X >= 1) -- looks like the right likelihood for a child,
// since a child with zero reads is not observed. It is the wrong one here, because the quantity it
// is compared against (`a_ind * p_size`) is an expected *number* of neighbouring molecules, not a
// probability conditioned on one existing. Dividing by (1 - e^-lambda) cancels precisely the term
// that says whether an error child should exist at all -- and for a singleton child that term is
// the entire signal: ZT-Poisson(1, lambda) -> 1 for every small lambda, so the error rate, and
// with it the barcode's own base quality, stops mattering at exactly the coverage where nothing
// else is available either.
double poisson_pmf(uint32_t k, double lambda) {
    if (lambda <= 0.0) return 0.0;
    double logp = -lambda + k * std::log(lambda);
    for (uint32_t i = 2; i <= k; ++i) logp -= std::log(static_cast<double>(i));
    return std::exp(logp);
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
    for (int j = 0; j < length; ++j) {
        const double m = collision(j);
        // A sample that got no reads has no composition, so every m_j is 0 and the sum is +inf.
        // Infinity is not an effective length; it is the absence of one, and printing it into a
        // TSV column that everything downstream parses as a number is worse than saying zero.
        if (m <= 0.0) return 0.0;
        l -= std::log(m) / std::log(4.0);
    }
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

void UmiCounts::carry_evidence(int payload_width, bool quality) {
    if (!buf_.empty() || !entries_.empty()) {
        throw MigecError("UmiCounts::carry_evidence: the table already holds barcodes -- the "
                         "evidence has to be carried from the first add, or the barcodes added "
                         "before it would have none and nothing would say which");
    }
    ev_pw_ = std::max(0, payload_width);
    ev_err_ = quality ? length_ : 0;
}

void UmiCounts::push_evidence(const float* pos_err, std::string_view payload) {
    if (ev_averaged_) {
        throw MigecError("UmiCounts: a barcode was added after the table had been read. The "
                         "carried quality is accumulated as a sum and divided by the counts when "
                         "the table is handed over, so a later add would fold a raw sum into an "
                         "average");
    }
    if (ev_err_ > 0) {
        const size_t at = buf_err_.size();
        buf_err_.resize(at + static_cast<size_t>(ev_err_), 0.0f);
        if (pos_err) std::copy_n(pos_err, ev_err_, buf_err_.begin() + static_cast<long>(at));
    }
    if (ev_pw_ > 0) {
        const size_t at = buf_pay_.size();
        buf_pay_.resize(at + static_cast<size_t>(ev_pw_), 0);
        const size_t n = std::min(payload.size(), static_cast<size_t>(ev_pw_));
        std::copy_n(payload.begin(), n, buf_pay_.begin() + static_cast<long>(at));
    }
}

// The same flush, for a table that carries evidence. Two differences, both deliberate:
//
//   * the buffer is sorted through a PERMUTATION and stably, because the evidence lives in
//     parallel arrays and because the payload draft is the first read's -- an unstable sort would
//     make the draft depend on where the buffer boundaries happened to fall;
//   * the buffer is gathered into key order first, because the merge below reads it sequentially
//     and a permutation cannot be followed backwards through three arrays at once. The gather is
//     over the BUFFER, which is half the entry array at most.
//
// Never: it still merges BACKWARDS into the grown arrays. Merging forwards into fresh ones is
// three lines shorter and doubles the peak memory of exactly the table this machinery exists to
// bound -- measured at 1,788 MB against 441 MB on 1.5 M barcodes at a 256-base draft.
void UmiCounts::flush_evidence() const {
    const size_t w0 = buf_.size();
    const size_t E = static_cast<size_t>(ev_err_);
    const size_t P = static_cast<size_t>(ev_pw_);

    // Sort the buffer through a permutation and STABLY: the evidence lives in parallel arrays, and
    // the payload draft is the first read's, so an unstable sort would make the draft depend on
    // where the buffer boundaries happened to fall.
    {
        std::vector<uint32_t> ord(w0);
        for (size_t i = 0; i < w0; ++i) ord[i] = static_cast<uint32_t>(i);
        std::stable_sort(ord.begin(), ord.end(),
                         [this](uint32_t a, uint32_t b) { return buf_[a].key < buf_[b].key; });
        std::vector<Entry> sorted;
        std::vector<float> sorted_err;
        std::vector<char> sorted_pay;
        sorted.reserve(w0);
        sorted_err.reserve(buf_err_.size());
        sorted_pay.reserve(buf_pay_.size());
        for (uint32_t k : ord) {
            // Run-length reduce while gathering: the same barcode read twice in one buffer.
            if (!sorted.empty() && sorted.back().key == buf_[k].key) {
                sorted.back().count += buf_[k].count;
                for (size_t j = 0; j < E; ++j) {
                    sorted_err[(sorted.size() - 1) * E + j] += buf_err_[k * E + j];
                }
                if (P && sorted_pay[(sorted.size() - 1) * P] == 0) {
                    std::copy_n(buf_pay_.begin() + static_cast<long>(k * P), P,
                                sorted_pay.begin() + static_cast<long>((sorted.size() - 1) * P));
                }
                continue;
            }
            sorted.push_back(buf_[k]);
            if (E) sorted_err.insert(sorted_err.end(),
                                     buf_err_.begin() + static_cast<long>(k * E),
                                     buf_err_.begin() + static_cast<long>((k + 1) * E));
            if (P) sorted_pay.insert(sorted_pay.end(),
                                     buf_pay_.begin() + static_cast<long>(k * P),
                                     buf_pay_.begin() + static_cast<long>((k + 1) * P));
        }
        buf_.swap(sorted);
        buf_err_.swap(sorted_err);
        buf_pay_.swap(sorted_pay);
    }
    const size_t w = buf_.size();

    if (entries_.empty()) {
        entries_.swap(buf_);
        err_.swap(buf_err_);
        pay_.swap(buf_pay_);
        buf_.clear();
        buf_err_.clear();
        buf_pay_.clear();
        flush_at_ = std::min(buffer_cap(), std::max(kMinBuffer, entries_.size() / 2));
        if (spill_budget_ && entries_.size() * (sizeof(Entry) + ev_bytes()) > spill_budget_) {
            spill();
        }
        return;
    }

    const size_t n = entries_.size();
    // Note: geometric, unlike the count-only path's exact `reserve(n + w)`. The buffer is capped in
    // BYTES here, so on a wide draft it is small against the table and an exact reserve reallocates
    // (and copies) all three arrays on every flush. Measured at a 256-base draft over 1.5 M
    // barcodes: 5,857 MB peak against 1,812 MB. The slack is bounded by the spill budget, which is
    // what this path has and the other does not.
    auto grow = [](auto& v, size_t want) {
        if (v.capacity() < want) v.reserve(std::max(want, v.capacity() + v.capacity() / 2));
        v.resize(want);
    };
    grow(entries_, n + w);
    if (E) grow(err_, (n + w) * E);
    if (P) grow(pay_, (n + w) * P);

    // `out` only ever runs ahead of both read cursors, so a block written at `out` never clobbers
    // a source block still to be read.
    auto put = [&](size_t dst, bool from_entries, size_t src) {
        if (from_entries && dst == src) return;  // a run that did not move
        entries_[dst] = from_entries ? entries_[src] : buf_[src];
        if (E) {
            const std::vector<float>& s = from_entries ? err_ : buf_err_;
            std::copy_n(s.begin() + static_cast<long>(src * E), E,
                        err_.begin() + static_cast<long>(dst * E));
        }
        if (P) {
            const std::vector<char>& s = from_entries ? pay_ : buf_pay_;
            std::copy_n(s.begin() + static_cast<long>(src * P), P,
                        pay_.begin() + static_cast<long>(dst * P));
        }
    };

    size_t i = n, j = w, out = n + w;
    while (i > 0 && j > 0) {
        const Entry& a = entries_[i - 1];
        const Entry& b = buf_[j - 1];
        if (a.key == b.key) {
            --out;
            --i;
            --j;
            // The older run wins the payload draft, and that is what makes it the FIRST read's
            // however many times the buffer was flushed. Fold the newer one's quality in first,
            // then write the older one over the slot.
            const uint32_t total = a.count + b.count;
            if (E) {
                for (size_t q = 0; q < E; ++q) err_[i * E + q] += buf_err_[j * E + q];
            }
            if (P && pay_[i * P] == 0) {
                std::copy_n(buf_pay_.begin() + static_cast<long>(j * P), P,
                            pay_.begin() + static_cast<long>(i * P));
            }
            put(out, true, i);
            entries_[out].count = total;
        } else if (a.key > b.key) {
            put(--out, true, --i);
        } else {
            put(--out, false, --j);
        }
    }
    while (j > 0) put(--out, false, --j);
    // Equal keys collapsed, so the merged run can be shorter than the space reserved for it; the
    // survivors sit at the top and are shifted down.
    if (out > 0) {
        const size_t kept = n + w - out;
        for (size_t k = 0; k < kept; ++k) put(i + k, true, out + k);
        entries_.resize(i + kept);
        if (E) err_.resize((i + kept) * E);
        if (P) pay_.resize((i + kept) * P);
    }
    buf_.clear();
    buf_err_.clear();
    buf_pay_.clear();
    flush_at_ = std::min(buffer_cap(), std::max(kMinBuffer, entries_.size() / 2));
    if (spill_budget_ && entries_.size() * (sizeof(Entry) + ev_bytes()) > spill_budget_) spill();
}

void UmiCounts::flush() const {
    if (buf_.empty()) return;
    if (carries_evidence()) {
        flush_evidence();
        return;
    }
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
        flush_at_ = std::min(buffer_cap(), std::max(kMinBuffer, entries_.size() / 2));
    };

    // Never: the FIRST flush is the one that swaps the buffer in wholesale, and it has to check the
    // budget like any other. Returning early here meant a counter whose whole library arrives in
    // one buffer's worth never spilled at all -- the partition switched itself off on exactly the
    // small-and-simple case, and reported a resident answer that happened to be right.
    if (entries_.empty()) {
        entries_.swap(buf_);
        buf_.clear();
        set_next_flush();
        if (spill_budget_ && entries_.size() * sizeof(Entry) > spill_budget_) spill();
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
    if (spill_budget_ && entries_.size() * sizeof(Entry) > spill_budget_) spill();
}

void UmiCounts::enable_spill(const std::string& directory, size_t budget_bytes, int bits) {
    if (bits < 1 || bits > 20) {
        throw MigecError("UmiCounts::enable_spill: bits must be in 1..20, got " +
                         std::to_string(bits));
    }
    if (budget_bytes == 0) throw MigecError("UmiCounts::enable_spill: budget must be positive");
    // Correction runs a second pass on keys rotated by the width of the partitioned prefix, and
    // the two prefixes have to be disjoint or a substitution inside the first one is invisible to
    // both passes. The prefix is (bits + 1) / 2 bases, so it has to fit twice into the barcode.
    // Refused here, on the caller's thread, where the number is attributable -- not later, as a
    // merge count that is quietly short.
    if (2 * ((bits + 1) / 2) > length_) {
        throw MigecError("UmiCounts::enable_spill: " + std::to_string(bits) +
                         " partition bits need " + std::to_string(2 * ((bits + 1) / 2)) +
                         " barcode positions for the rotated correction pass, but the barcode is " +
                         std::to_string(length_) + " long");
    }
    spill_dir_ = directory;
    spill_budget_ = budget_bytes;
    spill_bits_ = bits;
    // Note: the directory is created by the first spill, not here. Most runs never reach the
    // budget, and a stage that mkdirs next to its output on every run has to explain itself.
}

void UmiCounts::require_resident(const char* what) const {
    if (!spill_paths_.empty()) {
        throw MigecError(std::string("UmiCounts: ") + what +
                         " needs every entry resident, but this counter has spilled to a range "
                         "partition -- that is the point of spilling. Use for_each(), which "
                         "streams one bucket at a time.");
    }
}

// The top `spill_bits_` bits of the key decide the bucket, so buckets are key-ordered ranges and
// concatenating them in index order is ascending key order.
void UmiCounts::spill() const {
    if (entries_.empty()) return;
    const size_t n_buckets = size_t{1} << spill_bits_;
    if (spill_paths_.empty()) {
        std::filesystem::create_directories(spill_dir_);
        spill_paths_.resize(n_buckets);
        if (ev_err_ > 0) spill_err_paths_.resize(n_buckets);
        if (ev_pw_ > 0) spill_pay_paths_.resize(n_buckets);
        std::error_code ec;
        auto claim = [&](std::vector<std::string>& into, size_t b, const char* stem) {
            into[b] = (std::filesystem::path(spill_dir_) /
                       (stem + std::to_string(b) + ".bin")).string();
            // Never: every spill APPENDS, so a bucket left behind by a run that died would be read
            // back as part of this library and its counts added to it -- silently, since a stale
            // bucket is a well-formed one. Only the files this counter is about to own are removed,
            // once, and nothing else in the directory is touched.
            std::filesystem::remove(into[b], ec);
        };
        for (size_t b = 0; b < n_buckets; ++b) {
            claim(spill_paths_, b, "umi_");
            // The evidence rides in its own file per bucket, on the same key and the same bit
            // count, so bucket b of one is bucket b of the other entry for entry. Separate files
            // rather than one interleaved record because a bucket is written as three contiguous
            // runs and interleaving would cost a write per entry.
            if (ev_err_ > 0) claim(spill_err_paths_, b, "err_");
            if (ev_pw_ > 0) claim(spill_pay_paths_, b, "pay_");
        }
    }
    // entries_ is already sorted, so each bucket is one contiguous run: find the boundaries and
    // append each run in one write rather than per entry.
    const int shift = 64 - spill_bits_;
    size_t i = 0;
    while (i < entries_.size()) {
        const uint64_t b = entries_[i].key >> shift;
        size_t j = i;
        while (j < entries_.size() && (entries_[j].key >> shift) == b) ++j;
        auto append_run = [&](const std::vector<std::string>& paths, const void* data,
                              size_t bytes) {
            const std::string& path = paths[static_cast<size_t>(b)];
            std::ofstream out(path, std::ios::binary | std::ios::app);
            if (!out) throw MigecError("UmiCounts: cannot write spill file " + path);
            out.write(static_cast<const char*>(data), static_cast<std::streamsize>(bytes));
            if (!out) throw MigecError("UmiCounts: short write to " + path);
        };
        append_run(spill_paths_, entries_.data() + i, (j - i) * sizeof(Entry));
        if (ev_err_ > 0) {
            append_run(spill_err_paths_, err_.data() + i * static_cast<size_t>(ev_err_),
                       (j - i) * static_cast<size_t>(ev_err_) * sizeof(float));
        }
        if (ev_pw_ > 0) {
            append_run(spill_pay_paths_, pay_.data() + i * static_cast<size_t>(ev_pw_),
                       (j - i) * static_cast<size_t>(ev_pw_));
        }
        i = j;
    }
    // Release the capacity, not just the size: shrinking to zero size while holding the array is
    // exactly the allocation being bounded.
    std::vector<Entry>().swap(entries_);
    std::vector<float>().swap(err_);
    std::vector<char>().swap(pay_);
    distinct_known_ = false;
}

void UmiCounts::for_each(const std::function<void(const Entry&)>& fn) const {
    for_each_bucket([&fn](const std::vector<Entry>& bucket) {
        for (const Entry& e : bucket) fn(e);
    });
}

void UmiCounts::for_each_bucket(const std::function<void(const std::vector<Entry>&)>& fn) const {
    flush();
    if (spill_paths_.empty()) {
        fn(entries_);
        return;
    }
    // Anything added since the last spill is still resident and belongs in its bucket.
    const int shift = 64 - spill_bits_;
    std::vector<Entry> bucket;
    for (size_t b = 0; b < spill_paths_.size(); ++b) {
        bucket.clear();
        std::ifstream in(spill_paths_[b], std::ios::binary);
        if (in) {
            in.seekg(0, std::ios::end);
            const std::streamoff bytes = in.tellg();
            in.seekg(0, std::ios::beg);
            if (bytes > 0) {
                // Never: a spill file is a whole number of entries. A run killed mid-write leaves
                // a partial one, and reading `bytes` into a buffer sized `bytes / 16` overruns the
                // heap by up to 15 bytes -- silently, since the read itself succeeds.
                if (static_cast<size_t>(bytes) % sizeof(Entry) != 0) {
                    throw MigecError("UmiCounts: spill file " + spill_paths_[b] + " holds " +
                                     std::to_string(bytes) +
                                     " bytes, which is not a whole number of entries -- it was "
                                     "truncated, most likely by a run that was killed");
                }
                bucket.resize(static_cast<size_t>(bytes) / sizeof(Entry));
                in.read(reinterpret_cast<char*>(bucket.data()), bytes);
                if (!in) throw MigecError("UmiCounts: short read from " + spill_paths_[b]);
            }
        }
        for (const Entry& e : entries_) {
            if ((e.key >> shift) == b) bucket.push_back(e);
        }
        if (bucket.empty()) continue;
        // A key can appear once per spill plus once resident, so reduction happens HERE rather
        // than at spill time -- which is what makes a spill O(1) per entry.
        std::sort(bucket.begin(), bucket.end(),
                  [](const Entry& a, const Entry& c) { return a.key < c.key; });
        size_t w = 0;
        for (size_t r = 0; r < bucket.size(); ++r) {
            if (w > 0 && bucket[w - 1].key == bucket[r].key) {
                bucket[w - 1].count += bucket[r].count;
            } else {
                bucket[w++] = bucket[r];
            }
        }
        bucket.resize(w);
        fn(bucket);
    }
}

namespace {

// A whole spill file, as a vector of T. Empty when the file is not there, which is what an empty
// bucket looks like.
template <typename T>
std::vector<T> read_spill(const std::string& path, const char* what) {
    std::vector<T> out;
    if (path.empty()) return out;
    std::ifstream in(path, std::ios::binary);
    if (!in) return out;
    in.seekg(0, std::ios::end);
    const std::streamoff bytes = in.tellg();
    in.seekg(0, std::ios::beg);
    if (bytes <= 0) return out;
    // Never: a spill file is a whole number of records. A run killed mid-write leaves a partial
    // one, and reading `bytes` into a buffer sized `bytes / sizeof(T)` overruns the heap --
    // silently, since the read itself succeeds.
    if (static_cast<size_t>(bytes) % sizeof(T) != 0) {
        throw MigecError(std::string("UmiCounts: spill file ") + path + " holds " +
                         std::to_string(bytes) + " bytes, which is not a whole number of " + what +
                         " -- it was truncated, most likely by a run that was killed");
    }
    out.resize(static_cast<size_t>(bytes) / sizeof(T));
    in.read(reinterpret_cast<char*>(out.data()), bytes);
    if (!in) throw MigecError(std::string("UmiCounts: short read from ") + path);
    return out;
}

}  // namespace

void UmiCounts::for_each_bucket(
    const std::function<void(const std::vector<Entry>&, const BarcodeEvidence&)>& fn) const {
    flush();
    const size_t E = static_cast<size_t>(ev_err_);
    const size_t P = static_cast<size_t>(ev_pw_);
    if (!carries_evidence()) {
        const BarcodeEvidence none;
        for_each_bucket([&](const std::vector<Entry>& bucket) { fn(bucket, none); });
        return;
    }

    // `position_error` accumulates as a SUM, because a sum is what survives a partition -- it adds
    // across spill generations and across the resident tail. The posterior wants the mean over the
    // barcode's reads, so the division happens here, once.
    auto average = [&](const std::vector<Entry>& bucket, std::vector<float>& err) {
        for (size_t i = 0; i < bucket.size(); ++i) {
            const float c = static_cast<float>(bucket[i].count);
            for (size_t j = 0; j < E; ++j) err[i * E + j] /= c;
        }
    };

    if (spill_paths_.empty()) {
        // Never: divided IN PLACE. Copying the evidence to divide it would momentarily double the
        // resident set of exactly the table this machinery exists to bound. Adding after the table
        // has been handed over is refused in `push_evidence` rather than silently folding a raw
        // sum into an average.
        if (E && !ev_averaged_) {
            average(entries_, err_);
            ev_averaged_ = true;
        }
        BarcodeEvidence ev;
        if (E) ev.position_error.assign(err_.begin(), err_.end());
        if (P) {
            ev.payload.assign(pay_.begin(), pay_.end());
            ev.payload_width = ev_pw_;
        }
        // A resident table hands the evidence over by value because the caller may hold it past
        // the call; it is one allocation on a path that runs once.
        fn(entries_, ev);
        return;
    }

    const int shift = 64 - spill_bits_;
    for (size_t b = 0; b < spill_paths_.size(); ++b) {
        std::vector<Entry> raw = read_spill<Entry>(spill_paths_[b], "entries");
        std::vector<float> raw_err =
            E ? read_spill<float>(spill_err_paths_[b], "quality floats") : std::vector<float>();
        std::vector<char> raw_pay =
            P ? read_spill<char>(spill_pay_paths_[b], "payload bases") : std::vector<char>();
        if (E && raw_err.size() != raw.size() * E) {
            throw MigecError("UmiCounts: spill bucket " + std::to_string(b) + " holds " +
                             std::to_string(raw.size()) + " entries but " +
                             std::to_string(raw_err.size()) +
                             " quality floats -- the two files disagree");
        }
        // Anything added since the last spill is still resident and belongs in its bucket. It goes
        // last, which is the newest run, so "older wins" still picks the first read's payload.
        for (size_t i = 0; i < entries_.size(); ++i) {
            if ((entries_[i].key >> shift) != b) continue;
            raw.push_back(entries_[i]);
            if (E) raw_err.insert(raw_err.end(), err_.begin() + static_cast<long>(i * E),
                                  err_.begin() + static_cast<long>((i + 1) * E));
            if (P) raw_pay.insert(raw_pay.end(), pay_.begin() + static_cast<long>(i * P),
                                  pay_.begin() + static_cast<long>((i + 1) * P));
        }
        if (raw.empty()) continue;

        // A key can appear once per spill plus once resident, so reduction happens HERE. Stable,
        // through a permutation, so the surviving payload draft is still the first read's.
        std::vector<uint32_t> ord(raw.size());
        for (size_t i = 0; i < ord.size(); ++i) ord[i] = static_cast<uint32_t>(i);
        std::stable_sort(ord.begin(), ord.end(),
                         [&raw](uint32_t x, uint32_t y) { return raw[x].key < raw[y].key; });

        std::vector<Entry> bucket;
        BarcodeEvidence ev;
        ev.payload_width = ev_pw_;
        bucket.reserve(raw.size());
        if (E) ev.position_error.reserve(raw.size() * E);
        if (P) ev.payload.reserve(raw.size() * P);
        for (uint32_t k : ord) {
            if (!bucket.empty() && bucket.back().key == raw[k].key) {
                bucket.back().count += raw[k].count;
                for (size_t j = 0; j < E; ++j) {
                    ev.position_error[(bucket.size() - 1) * E + j] += raw_err[k * E + j];
                }
                if (P && ev.payload[(bucket.size() - 1) * P] == 0) {
                    std::copy_n(raw_pay.begin() + static_cast<long>(k * P), P,
                                ev.payload.begin() + static_cast<long>((bucket.size() - 1) * P));
                }
                continue;
            }
            bucket.push_back(raw[k]);
            if (E) ev.position_error.insert(ev.position_error.end(),
                                            raw_err.begin() + static_cast<long>(k * E),
                                            raw_err.begin() + static_cast<long>((k + 1) * E));
            if (P) ev.payload.insert(ev.payload.end(),
                                     raw_pay.begin() + static_cast<long>(k * P),
                                     raw_pay.begin() + static_cast<long>((k + 1) * P));
        }
        if (E) average(bucket, ev.position_error);
        fn(bucket, ev);
    }
}

UmiCounts UmiCounts::rotated_copy(int r, const std::string& directory) const {
    UmiCounts out(length_, buffer_limit_);
    if (carries_evidence()) out.carry_evidence(ev_pw_, ev_err_ > 0);
    if (spill_budget_) out.enable_spill(directory, spill_budget_, spill_bits_);
    const size_t E = static_cast<size_t>(ev_err_);
    std::vector<float> sums(E);
    const size_t rotation = length_ > 0 ? static_cast<size_t>((r % length_ + length_) % length_) : 0;
    for_each_bucket([&](const std::vector<Entry>& bucket, const BarcodeEvidence& ev) {
        for (size_t i = 0; i < bucket.size(); ++i) {
            // The evidence arrives averaged and is stored as a sum, so it is multiplied back by
            // the count on the way in. Anything else divides by the count twice and the rotated
            // pass weighs a barcode quality that is a factor of its own depth too small.
            for (size_t j = 0; j < E; ++j) {
                sums[j] = ev.position_error[i * E + (j + rotation) % E] * static_cast<float>(bucket[i].count);
            }
            out.add(rotate_barcode(bucket[i].key, length_, r), bucket[i].count,
                    E ? sums.data() : nullptr,
                    ev.has_payload()
                        ? std::string_view(ev.payload.data() + i * static_cast<size_t>(ev_pw_),
                                           static_cast<size_t>(ev_pw_))
                        : std::string_view());
        }
    });
    return out;
}

size_t UmiCounts::distinct() const {
    flush();
    if (spill_paths_.empty()) return entries_.size();
    if (!distinct_known_) {
        uint64_t n = 0;
        for_each([&n](const Entry&) { ++n; });
        distinct_cache_ = n;
        distinct_known_ = true;
    }
    return static_cast<size_t>(distinct_cache_);
}

void UmiCounts::merge(const UmiCounts& other) {
    other.flush();
    other.require_resident("merge()");
    if (carries_evidence() || other.carries_evidence()) {
        throw MigecError("UmiCounts::merge: this counter carries per-barcode evidence, and merging "
                         "would have to reduce it too. Add the reads to one counter instead");
    }
    for (const Entry& e : other.entries_) add(e.key, e.count);
}

const uint32_t* UmiCounts::find(uint64_t key) const {
    flush();
    require_resident("find()");
    auto it = std::lower_bound(entries_.begin(), entries_.end(), key,
                               [](const Entry& e, uint64_t k) { return e.key < k; });
    if (it == entries_.end() || it->key != key) return nullptr;
    return &it->count;
}

size_t UmiCounts::memory_bytes() const {
    return entries_.capacity() * sizeof(Entry) + buf_.capacity() * sizeof(Entry) +
           (err_.capacity() + buf_err_.capacity()) * sizeof(float) +
           pay_.capacity() + buf_pay_.capacity();
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
    for_each([&h](const Entry& e) {
        const int b = log2_bin(e.count);
        h.reads[static_cast<size_t>(b)] += e.count;
        h.units[static_cast<size_t>(b)] += 1;
    });
    return h;
}

UmiComposition UmiCounts::composition(bool weight_by_reads) const {
    UmiComposition c;
    c.length = length_;
    c.freq.assign(static_cast<size_t>(length_), {0.0, 0.0, 0.0, 0.0});
    double total = 0.0;
    const int L = length_;
    for_each([&](const Entry& e) {
        const double w = weight_by_reads ? static_cast<double>(e.count) : 1.0;
        for (int j = 0; j < L; ++j) {
            const uint8_t code = static_cast<uint8_t>((e.key >> (62 - 2 * j)) & 3u);
            c.freq[static_cast<size_t>(j)][code] += w;
        }
        total += w;
    });
    if (total > 0.0) {
        for (auto& row : c.freq) {
            for (double& v : row) v /= total;
        }
    }
    return c;
}

namespace {

// Distinct MIG sizes and how many barcodes carry each, ascending. The expectation below sums over
// it rather than over the entries, which is what lets the sum be accumulated bucket by bucket --
// and it is kept in size order rather than in a hash order so the sum is the same everywhere.
using SizeHist = std::vector<std::pair<uint32_t, uint64_t>>;

SizeHist size_hist_of(const std::vector<UmiCounts::Entry>& m) {
    std::unordered_map<uint32_t, uint64_t> h;
    for (const UmiCounts::Entry& e : m) ++h[e.count];
    SizeHist out(h.begin(), h.end());
    std::sort(out.begin(), out.end());
    return out;
}

// Keep the established distance-one excess estimator distinct from the correction radius.
// Partition ownership remains position-based, and each unordered pair contributes once.
uint64_t d1_from_search(const std::vector<UmiCounts::Entry>& entries,
                        const BarcodeSearch& search, int j0, int j1, int threads) {
    uint64_t total = 0;
    search.consume(1, threads, [&](size_t, size_t, seqtree::CountBuffer&& block) {
        for (size_t row = 0; row < block.data.size(); row += 6) {
            const size_t query = block.data[row], reference = block.data[row + 1];
            if (block.data[row + 3] != 1 || reference <= query) continue;
            const uint64_t differences = entries[query].key ^ entries[reference].key;
            for (int j = j0; j < j1; ++j) {
                if (((differences >> (62 - 2 * j)) & 3u) != 0) { ++total; break; }
            }
        }
    });
    return total;
}

uint64_t d1_census(const std::vector<UmiCounts::Entry>& entries, int length,
                   int j0, int j1, int threads) {
    if (entries.size() < 2 || j1 <= j0) return 0;
    const BarcodeSearch search(entries, length);
    return d1_from_search(entries, search, j0, j1, threads);
}

// Inverts the distance-1 census into a per-base error rate. Split out of `estimate_umi_error` so
// the bucketed driver can feed it a census summed over buckets and both rotations: the census is
// the only part that needs the table, and everything here is a function of the library totals.
double solve_umi_error(uint64_t d1_obs, double n, const UmiComposition& comp, int L,
                       const std::vector<std::pair<uint32_t, uint64_t>>& sizes) {
    if (L <= 0 || n < 2.0) return 0.0;

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
        for (const auto& kv : sizes) {
            const double t = 1.0 - std::exp(-static_cast<double>(kv.first) * eps / 3.0);
            const double w = static_cast<double>(kv.second);
            parent_child += w * t;
            sibling += w * t * t;
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

}  // namespace

double estimate_umi_error(const UmiCounts& counts, const UmiComposition& comp, int threads) {
    const int L = counts.length();
    if (L <= 0 || counts.distinct() < 2) return 0.0;
    // Every position, because this table is the whole library: `entries()` refuses a spilled
    // counter, and a partition is censused by `correct_umis` in two passes instead.
    const std::vector<UmiCounts::Entry>& m = counts.entries();
    return solve_umi_error(d1_census(m, L, 0, L, threads), static_cast<double>(m.size()), comp, L,
                           size_hist_of(m));
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

namespace {

// log C(n,k) + the binomial pmf, in logs. d is small and n is a read length, so this is nowhere
// near a hot path -- it runs once per candidate parent, not per base.
double log_binom_pmf(int d, int n, double p) {
    if (n <= 0) return 0.0;
    p = std::clamp(p, 1e-9, 1.0 - 1e-9);
    return std::lgamma(n + 1.0) - std::lgamma(d + 1.0) - std::lgamma(n - d + 1.0) +
           d * std::log(p) + (n - d) * std::log1p(-p);
}

}  // namespace

namespace {

// Correction over one sorted, reduced barcode array. The array is the whole library on a resident
// counter and one range-partition bucket on a spilled one, which is the only difference between
// the two: everything the posterior needs that a bucket cannot see for itself -- the library's
// distinct count, its collision probability, its effective space and its error rate -- arrives
// through `params`, and `params.scan_from/scan_to` say which barcode positions this call owns.
// `library_sizes` is the MIG size distribution of the whole library, or null when this array is
// the whole library. Never: a bucket's own size distribution is NOT the library's -- it is the
// same shape drawn a few thousand times, and the posterior weighs an exact per-size probability
// against the error hypothesis, so estimating it per bucket moves borderline decisions. Measured:
// one merge in 547 differed from the resident answer until this was threaded through.
// How often two UNRELATED barcodes carry the same payload anyway -- the library's clonality, which
// is exactly what payload agreement is worth: log(1/clonality). Deterministic sampling, a fixed
// stride over the payloads it was handed, so the answer does not depend on an RNG seed and two
// runs of the pipeline agree.
double payload_agreement(const std::vector<char>& payload, int pw,
                         const CorrectionParams& params) {
    if (pw <= 0 || params.payload_null_samples <= 0) return 1.0;
    const size_t n_entries = payload.size() / static_cast<size_t>(pw);
    if (n_entries < 3) return 1.0;
    const size_t samples =
        std::min<size_t>(static_cast<size_t>(params.payload_null_samples), n_entries * 4);
    const size_t stride = std::max<size_t>(1, n_entries / 977 + 1);
    uint64_t same = 0, tried = 0;
    for (size_t s = 0; s < samples; ++s) {
        const size_t a = (s * 7919) % n_entries;
        const size_t b = (a + stride * (1 + s % 97)) % n_entries;
        if (a == b) continue;
        int mism = 0, cmp = 0;
        for (int j = 0; j < pw; ++j) {
            const char x = payload[a * static_cast<size_t>(pw) + static_cast<size_t>(j)];
            const char y = payload[b * static_cast<size_t>(pw) + static_cast<size_t>(j)];
            if (x == 0 || y == 0 || x == 'N' || y == 'N') continue;
            ++cmp;
            mism += x != y;
        }
        if (cmp < 8) continue;
        ++tried;
        same += static_cast<double>(mism) <= params.payload_same_fraction * cmp;
    }
    return tried ? std::max(static_cast<double>(same) / static_cast<double>(tried),
                            1.0 / static_cast<double>(tried + 1)) : 1.0;
}

// Payload drafts sampled evenly along the table, in key order, bounded.
//
// Never: sampling by INDEX INTO THE ARRAY measures the array, not the library. Two partitionings
// of one library hand the sampler different arrays, so an index rule draws different pairs, the
// clonality lands a few percent apart, and a borderline merge follows it -- which makes the output
// depend on the memory budget. Barcodes arrive in KEY order whatever the partition, so a rule
// written on that order sees the same barcodes either way, and a run that spilled and a run that
// did not agree byte for byte.
//
// Keeps every `stride_`-th barcode and doubles the stride whenever the buffer fills, so what is
// held is always "every k-th barcode in key order" for the k it ended on -- the same set however
// many times it doubled. Under kMax barcodes that is all of them.
class PayloadReservoir {
public:
    explicit PayloadReservoir(int width) : pw_(width) {}

    void add(const char* payload) {
        if (pw_ <= 0) return;
        const uint64_t idx = seen_++;
        if (idx % stride_) return;
        kept_.insert(kept_.end(), payload, payload + pw_);
        if (kept_.size() / static_cast<size_t>(pw_) < kMax) return;
        // Half of what is held is at an index divisible by twice the stride: keep those.
        const size_t w = static_cast<size_t>(pw_);
        size_t out = 0;
        for (size_t i = 0; i < kept_.size() / w; i += 2, ++out) {
            std::copy_n(kept_.begin() + static_cast<long>(i * w), w,
                        kept_.begin() + static_cast<long>(out * w));
        }
        kept_.resize(out * w);
        stride_ *= 2;
    }

    double clonality(const CorrectionParams& params) const {
        return payload_agreement(kept_, pw_, params);
    }

private:
    static constexpr size_t kMax = 8192;  // 256 kB at a 32 nt draft
    int pw_;
    uint64_t seen_ = 0, stride_ = 1;
    std::vector<char> kept_;
};

// One barcode's best parent, by KEY, from one pass over one array. This is what a bucketed run
// collects instead of merging: a barcode can have a candidate parent in each pass -- one
// substitution in the positions the partition leaves alone, one in the positions it hides -- and
// merging in the first pass would take the first candidate rather than the best.
struct Proposal {
    uint64_t child, parent;
    uint32_t child_count, parent_count;
    double posterior;
};

CorrectionResult correct_entries(const std::vector<UmiCounts::Entry>& m, int L,
                                 const CorrectionParams& params, const BarcodeEvidence& evidence,
                                 const std::map<uint32_t, uint64_t>* library_sizes = nullptr,
                                 std::vector<Proposal>* proposals = nullptr,
                                 const BarcodeSearch* supplied = nullptr,
                                 uint64_t* observed_d1 = nullptr) {
    CorrectionResult res;
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

    double eps = params.sequencing_error;
    if (eps <= 0.0) eps = 1e-4;  // a floor, so correction still runs on a clean small library
    res.estimated_error = eps;

    const double p_coll = params.library_collision;
    // The library's distinct count, which on a bucket is not this array's. Both are needed and
    // they are not interchangeable: `n` is how many molecules could have landed on a neighbouring
    // barcode, `n_local` is how many barcodes this call can see and is what the empirical MIG size
    // distribution below is estimated from.
    const double n =
        params.library_distinct ? static_cast<double>(params.library_distinct)
                                : static_cast<double>(n_entries);
    const double n_local = static_cast<double>(n_entries);
    const double space = params.library_space;
    res.saturated = space > 0.0 && n > 0.05 * space;

    // Barcode positions this call owns. All of them for a resident table; for a bucket, the ones
    // whose substitutions cannot have moved the barcode out of the bucket.
    const int j_from = std::max(0, params.scan_from);
    const int j_to = params.scan_to < 0 ? L : std::min(params.scan_to, L);

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
    const double n_sizes = library_sizes ? n : n_local;
    if (library_sizes) {
        for (const auto& kv : *library_sizes) {
            size_pmf[kv.first] = static_cast<double>(kv.second) / n_sizes;
        }
    } else {
        for (const UmiCounts::Entry& e : m) size_pmf[e.count] += 1.0;
        for (auto& kv : size_pmf) kv.second /= n_sizes;
    }
    const double size_floor = 1.0 / (n_sizes + 1.0);  // never claim a size is impossible

    // How often do two UNRELATED barcodes carry the same payload anyway? That is the library's
    // clonality, and it is exactly what payload agreement is worth: log(1/clonality). In a diverse
    // repertoire two random molecules never match and agreement is decisive; in a clonal library
    // they always match and it says nothing. Measured from the data rather than assumed, so the
    // evidence self-calibrates to the library it is given.
    const int pw = evidence.has_payload() ? evidence.payload_width : 0;
    double clonality = 1.0;
    if (params.library_clonality > 0.0) {
        // A bucket's own sample is a sample of the bucket: its barcodes share a key prefix, and
        // what payload agreement is worth is a property of the library. The bucketed driver
        // measures it over the whole partition and passes it down.
        clonality = params.library_clonality;
    } else if (pw > 0 && n_entries > 2) {
        clonality = payload_agreement(evidence.payload, pw, params);
    }
    res.payload_clonality = pw > 0 ? clonality : 0.0;

    // Order by count descending, then inspect native candidates from the smallest MIG
    // upwards. Indices, not copies of the entries: 4 bytes each rather than 16.
    std::vector<uint32_t> order(n_entries);
    for (size_t i = 0; i < n_entries; ++i) order[i] = static_cast<uint32_t>(i);
    std::sort(order.begin(), order.end(), [&m](uint32_t a, uint32_t b) {
        if (m[a].count != m[b].count) return m[a].count > m[b].count;
        return m[a].key < m[b].key;  // total order, so the result is reproducible
    });

    // ------------------------------------------------------------------ scan, in parallel
    //
    // Every child's best parent is a pure function of the table and the evidence: the loop below
    // reads `m`, `evidence`, `size_pmf` and the constants above, and NOTHING it writes is read by
    // another child. Candidate likelihoods and payload comparisons therefore parallelise
    // exactly, and the answer is identical at any thread
    // count. The decisions are applied afterwards, serially, in the original order.
    struct Decision {
        size_t parent = static_cast<size_t>(-1);
        double posterior = 0.0;
    };
    std::vector<Decision> decisions(n_entries);
    std::unique_ptr<BarcodeSearch> owned_search;
    if (!supplied) {
        owned_search = std::make_unique<BarcodeSearch>(m, L);
        supplied = owned_search.get();
    }
    auto canonical = [&](uint64_t key) {
        return params.coordinate_rotation ? rotate_barcode(key, L, L - params.coordinate_rotation) : key;
    };
    auto candidate_rank = [&](uint64_t child, uint64_t candidate) {
        child = canonical(child);
        candidate = canonical(candidate);
        const uint64_t difference = child ^ candidate;
        int position = 0;
        while (position < L && ((difference >> (62 - 2 * position)) & 3u) == 0) ++position;
        return std::make_pair(position, candidate);
    };
    auto consider = [&](size_t child_idx, const seqtree::Hit& hit) {
        const uint64_t child = m[child_idx].key;
        const uint32_t c_child = m[child_idx].count;
        auto& decision = decisions[child_idx];
        size_t& best_parent = decision.parent;
        double& best_post = decision.posterior;
        if (hit.n_subs == 0) return;
        const size_t cand_idx = hit.ref_id;
        const uint64_t cand = m[cand_idx].key;
        int j = 0;
        const uint64_t differences = child ^ cand;
        while (j < L && ((differences >> (62 - 2 * j)) & 3u) == 0) ++j;
        if (j < j_from || j >= j_to) return;
        const uint32_t c_par = m[cand_idx].count;

        // The payload likelihood ratio, before the count gates -- because it is what
        // decides whether those gates apply at all.
        double lr_payload = 1.0;
        bool payload_decisive = false, payload_refutes = false;
        if (pw > 0) {
            int mism = 0, cmp = 0;
            for (int q = 0; q < pw; ++q) {
                const char x = evidence.payload[child_idx * static_cast<size_t>(pw) +
                                                static_cast<size_t>(q)];
                const char y = evidence.payload[cand_idx * static_cast<size_t>(pw) +
                                                static_cast<size_t>(q)];
                if (x == 0 || y == 0 || x == 'N' || y == 'N') continue;
                ++cmp;
                mism += x != y;
            }
            if (cmp >= 8) {
                // Same molecule: two drafts of it disagree at ~2e. Independent molecule:
                // the same thing with probability `clonality`, and otherwise a different
                // sequence altogether.
                const double ll_same = log_binom_pmf(mism, cmp, params.payload_error);
                const double ll_diff = log_binom_pmf(mism, cmp, 0.75);
                const double ll_ind =
                    std::log(clonality * std::exp(ll_same) +
                             (1.0 - clonality) * std::exp(ll_diff) + 1e-300);
                lr_payload = std::exp(std::clamp(ll_same - ll_ind, -60.0, 60.0));
                payload_decisive = lr_payload > 10.0;
                payload_refutes = lr_payload < 0.1;
            }
        }
        // A payload that disagrees is not this molecule, whatever the counts say.
        if (payload_refutes) return;

        // The count gates exist because a child is smaller than its parent -- true, and
        // vacuous at 1-3 reads per UMI. They are lifted exactly when the reads themselves
        // say the two barcodes carry the same molecule.
        if (!payload_decisive) {
            if (c_par <= c_child) return;
            if (static_cast<double>(c_child) >
                params.max_child_fraction * static_cast<double>(c_par)) {
                return;
            }
        } else if (c_par < c_child) {
            return;  // still orient the merge into the larger of the two
        } else if (c_par == c_child && canonical(cand) > canonical(child)) {
            return;  // a tie folds into the lexicographically smaller key, once
        }

        // Two ways to be an error child. Sequencing miscalls land on one specific
        // alternative base, so the rate per neighbour is eps/3, not eps.
        // The barcode's OWN reported quality at the base that differs, when it is known:
        // a miscall carries a low Phred there and an early-PCR child carries a high one,
        // which is the distinction a single global rate has to average away.
        // A specific two-base miscall requires both errors, not either error.
        // This preserves the distance-one model exactly and prevents a single poor
        // base from licensing another confidently observed mismatch.
        double error_probability = 1.0, single_rate = eps;
        for (int position = 0; position < L; ++position) {
            if (((differences >> (62 - 2 * position)) & 3u) == 0) continue;
            double rate = eps;
            if (evidence.has_quality()) {
                const size_t at = child_idx * static_cast<size_t>(L) + static_cast<size_t>(position);
                if (at < evidence.position_error.size() && evidence.position_error[at] > 0.0f)
                    rate = std::clamp(static_cast<double>(evidence.position_error[at]), 1e-6, 0.75);
            }
            single_rate = rate;
            error_probability *= rate / 3.0;
        }
        const double lam = hit.n_subs == 1 ? static_cast<double>(c_par) * single_rate / 3.0
                                              : static_cast<double>(c_par) * error_probability;
        const double l_seq = poisson_pmf(c_child, lam);
        const double l_pol = ld_pmf(c_child, c_par, params.max_child_fraction);
        // Independent polymerase events are a model assumption for multi-base children.
        const double rho_events = std::pow(rho_pol, hit.n_subs);
        const double l_err = (1.0 - rho_events) * l_seq + rho_events * l_pol;

        // ...against being a real molecule that happens to sit one substitution away and
        // to have this many reads.
        auto sp = size_pmf.find(c_child);
        const double p_size = sp == size_pmf.end() ? size_floor
                                                   : std::max(sp->second, size_floor);
        const double l_ind = std::max(a_ind * p_size, 1e-300);

        // The payload evidence multiplies the error hypothesis, because it is a
        // likelihood ratio between exactly the two hypotheses already being weighed.
        const double post = (l_err * lr_payload) / (l_err * lr_payload + l_ind);
        if (post > best_post ||
            (post == best_post && best_parent != static_cast<size_t>(-1) &&
             candidate_rank(child, cand) < candidate_rank(child, m[best_parent].key))) {
            best_post = post;
            best_parent = cand_idx;
        }
    };
    // The verified development extension applies to bulk UMI cohorts. Preserve the
    // released cell-tagged model until its conditional background is implemented.
    const int substitutions = params.cell_length ? 1 : 2;
    supplied->consume(substitutions, params.threads, [&](size_t, size_t, seqtree::CountBuffer&& block) {
        // seqtree serializes consumer callbacks within its one native worker pool.
        // Every query's complete candidate set is consumed; only O(queries) decisions survive.
        for (size_t row = 0; row < block.data.size(); row += 6) {
            const size_t query = block.data[row], reference = block.data[row + 1];
            const auto substitutions = static_cast<uint16_t>(block.data[row + 3]);
            if (observed_d1 && substitutions == 1 && reference > query) ++*observed_d1;
            consider(query, seqtree::Hit{static_cast<uint32_t>(reference),
                static_cast<int32_t>(block.data[row + 2]), substitutions,
                static_cast<uint16_t>(block.data[row + 4]), static_cast<uint16_t>(block.data[row + 5])});
        }
    });

    // ------------------------------------------------------------------ apply, serially
    //
    // In the original order, because merges CHAIN: a child folds into a parent that has itself
    // already folded into someone else, and which root it lands on depends on what happened
    // before it. This part is cheap -- a union-find update per merged barcode -- and making it
    // parallel would change the answer, not the wall clock.
    for (auto it = order.rbegin(); it != order.rend(); ++it) {
        const size_t child_idx = *it;
        const uint32_t c_child = m[child_idx].count;
        const size_t best_parent = decisions[child_idx].parent;
        const double best_post = decisions[child_idx].posterior;

        if (best_post >= params.min_posterior && best_parent != static_cast<size_t>(-1)) {
            if (proposals) {
                // Scan only. The apply is the caller's, over every pass at once.
                proposals->push_back(Proposal{m[child_idx].key, m[best_parent].key, c_child,
                                              m[best_parent].count, best_post});
                continue;
            }
            if (m[best_parent].count <= c_child) ++res.merged_by_payload;
            // Follow the parent to its current root -- it may itself already have been merged.
            uint32_t root = static_cast<uint32_t>(best_parent);
            for (int guard = 0; guard < 64 && res.root[root] != root; ++guard) root = res.root[root];
            if (root == child_idx) continue;  // never make a cycle
            res.root[child_idx] = root;
            res.corrected[root] += res.corrected[child_idx];
            // This barcode's OWN reads, not its running total. Merges chain -- x folds into y and
            // y later folds into z -- and by then `corrected[y]` already carries x's reads, which
            // were counted when x moved. `merged_reads` is reads whose barcode changed, and each
            // read changes barcode once.
            res.merged_reads += m[child_idx].count;
            res.corrected[child_idx] = 0;
            ++res.merged;
        }
    }

    if (proposals) return res;  // the scalars are the caller's to compute, over the whole library

    // Flatten. A child merged early can point at a parent that was itself merged later in the
    // walk, so the invariant "root[root[i]] == root[i]" only holds after this pass. The read
    // counts were already correct -- they follow the chain as it forms -- but a consumer reading
    // `root` directly would otherwise land on an intermediate.
    for (size_t i = 0; i < n_entries; ++i) {
        uint32_t r = res.root[i];
        for (int guard = 0; guard < 64 && res.root[r] != r; ++guard) r = res.root[r];
        res.root[i] = r;
        // ...and the same thing by KEY, which is what survives a range partition. `root` indexes
        // `entries()`, and a bucketed run has no such array.
        if (r != i) {
            res.merges.push_back(
                CorrectionResult::Merge{m[i].key, m[r].key, m[i].count});
        }
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

// Collision probability and effective space of a whole library, from its composition.
double collision_of(const UmiComposition& comp, int L) {
    double p = 1.0;
    for (int j = 0; j < L; ++j) p *= comp.collision(j);
    return p;
}

// Correction over three disjoint prefix partitions. At most two substitutions cannot
// change every prefix, so every candidate pair co-occurs in at least one bucket.
// Proposals are deduplicated before merging; the distance-one census retains its
// separate two-pass ownership rule. Working memory includes the native index/hits
// for the current bucket and the temporary counters, not only the entry table.
CorrectionResult correct_spilled(const UmiCounts& counts, const CorrectionParams& params) {
    const int L = counts.length();
    const int pb = counts.spill_prefix_bases();
    // Two substitutions leave at least one of three DISJOINT prefixes untouched.
    // Older caller-selected partitions may be wider: re-partition without losing evidence.
    if (3 * pb > L) {
        const auto directory = (std::filesystem::path(counts.spill_dir()) / "three_prefixes").string();
        UmiCounts narrower(L);
        if (counts.carries_evidence()) narrower.carry_evidence(counts.payload_width(), true);
        const int bits = std::min(counts.spill_bits(), 2 * (L / 3));
        if (bits > 0) narrower.enable_spill(directory, counts.spill_budget(), bits);
        // With L < 3 the complete alphabet has at most 16 keys and is bounded without spilling.
        const size_t P = static_cast<size_t>(counts.payload_width());
        std::vector<float> sums(static_cast<size_t>(L));
        counts.for_each_bucket([&](const auto& bucket, const BarcodeEvidence& ev) {
            for (size_t i = 0; i < bucket.size(); ++i) {
                for (int j = 0; j < L && ev.has_quality(); ++j)
                    sums[static_cast<size_t>(j)] = ev.position_error[i * L + j] * bucket[i].count;
                narrower.add(bucket[i].key, bucket[i].count, ev.has_quality() ? sums.data() : nullptr,
                    ev.has_payload() ? std::string_view(ev.payload.data() + i * P, P) : std::string_view());
            }
        });
        auto result = correct_umis(narrower, params);
        result.root.clear();
        result.corrected.clear();
        std::error_code error;
        std::filesystem::remove_all(directory, error);
        return result;
    }
    CorrectionResult res;

    // One pass over the spilled buckets does four things at once, because each is a full read of
    // the partition off disk: the library composition, its distinct count, the census of the pairs
    // this partition can see, and the rotated copy the second pass needs.
    UmiComposition comp;
    comp.length = L;
    comp.freq.assign(static_cast<size_t>(L), {0.0, 0.0, 0.0, 0.0});
    double n = 0.0;
    uint64_t d1 = 0;
    std::map<uint32_t, uint64_t> sizes;  // ordered: the expectation sums over it
    const std::string rot_dir =
        (std::filesystem::path(counts.spill_dir()) / "rotated").string();
    UmiCounts rotated(L);
    if (counts.carries_evidence()) {
        rotated.carry_evidence(counts.payload_width(), true);
    }
    rotated.enable_spill(rot_dir, counts.spill_budget(), counts.spill_bits());
    const size_t E = static_cast<size_t>(L);
    const size_t P = static_cast<size_t>(counts.payload_width());
    PayloadReservoir reservoir(counts.payload_width());
    std::vector<float> sums(E);

    counts.for_each_bucket([&](const std::vector<UmiCounts::Entry>& bucket,
                               const BarcodeEvidence& ev) {
        for (size_t i = 0; i < bucket.size(); ++i) {
            const UmiCounts::Entry& e = bucket[i];
            for (int j = 0; j < L; ++j) {
                const size_t code = static_cast<size_t>((e.key >> (62 - 2 * j)) & 3u);
                comp.freq[static_cast<size_t>(j)][code] += 1.0;
            }
            n += 1.0;
            ++sizes[e.count];
            // The evidence rotates with the key. It arrives averaged over the barcode's reads and
            // is stored as a sum, so it is multiplied back on the way in -- otherwise the second
            // pass weighs a barcode quality a factor of its own depth too small.
            if (ev.has_quality()) {
                for (size_t j = 0; j < E; ++j) {
                    sums[j] = ev.position_error[i * E + (j + static_cast<size_t>(pb)) % E] * static_cast<float>(e.count);
                }
            }
            rotated.add(rotate_barcode(e.key, L, pb), e.count,
                        ev.has_quality() ? sums.data() : nullptr,
                        ev.has_payload()
                            ? std::string_view(ev.payload.data() + i * P, P)
                            : std::string_view());
            // Clonality is a property of the LIBRARY, and the buckets arrive in key order, so the
            // reservoir sees the same barcodes a resident table would.
            if (ev.has_payload()) reservoir.add(ev.payload.data() + i * P);
        }
        d1 += d1_census(bucket, L, pb, L, params.threads);
    });
    if (n > 0.0) {
        for (auto& row : comp.freq) {
            for (double& v : row) v /= n;
        }
    }
    rotated.for_each_bucket([&](const std::vector<UmiCounts::Entry>& bucket) {
        d1 += d1_census(bucket, L, L - pb, L, params.threads);
    });

    CorrectionParams sub = params;
    if (counts.payload_width() > 0) {
        sub.library_clonality = reservoir.clonality(params);
        res.payload_clonality = sub.library_clonality;
    }
    sub.library_distinct = static_cast<size_t>(n);
    sub.library_collision = collision_of(comp, L);
    sub.library_space = comp.effective_space();
    if (sub.sequencing_error < 0.0) {
        sub.sequencing_error =
            solve_umi_error(d1, n, comp, L, SizeHist(sizes.begin(), sizes.end()));
    }
    // Never: the same floor the buckets are about to apply, applied HERE too, or the rate reported
    // is not the rate used. A clean library has no distance-1 excess, the census returns 0, each
    // bucket quietly corrects at 1e-4 -- and the summary said 0.0, which then propagated into the
    // error budget as a ratio of zero. Caught by diffing a whole run against its resident twin.
    if (sub.sequencing_error <= 0.0) sub.sequencing_error = 1e-4;

    // Each pair within two substitutions shares at least one of the three prefixes.
    // Scan all positions in each bucket; combine duplicate proposals before any merge.
    // The established two-pass distance-one census above still counts each pair once.
    const auto third_dir = (std::filesystem::path(counts.spill_dir()) / "rotated_twice").string();
    UmiCounts third = counts.rotated_copy(2 * pb, third_dir);
    std::vector<Proposal> p1, p2, p3;
    sub.scan_from = 0;
    sub.scan_to = L;
    auto scan = [&](const UmiCounts& table, int rotation, std::vector<Proposal>& proposals) {
        sub.coordinate_rotation = rotation;
        table.for_each_bucket([&](const auto& bucket, const BarcodeEvidence& ev) {
            const size_t begin = proposals.size();
            correct_entries(bucket, L, sub, ev, &sizes, &proposals);
            if (rotation) {
                for (size_t i = begin; i < proposals.size(); ++i) {
                    proposals[i].child = rotate_barcode(proposals[i].child, L, L - rotation);
                    proposals[i].parent = rotate_barcode(proposals[i].parent, L, L - rotation);
                }
            }
        });
    };
    scan(counts, 0, p1);
    scan(rotated, pb, p2);
    scan(third, 2 * pb, p3);
    auto rank = [&](const Proposal& p) {
        const uint64_t difference = p.child ^ p.parent;
        int position = 0;
        while (position < L && ((difference >> (62 - 2 * position)) & 3u) == 0) ++position;
        return std::make_pair(position, p.parent);
    };
    std::unordered_map<uint64_t, Proposal> best;
    best.reserve(p1.size() + p2.size() + p3.size());
    for (const auto* proposals : {&p1, &p2, &p3}) {
        for (const Proposal& p : *proposals) {
            auto it = best.find(p.child);
            if (it == best.end() || p.posterior > it->second.posterior ||
                (p.posterior == it->second.posterior && rank(p) < rank(it->second))) best[p.child] = p;
        }
    }

    // Apply, serially, smallest MIG first and ties by descending key -- the order the resident
    // walk applies in, because merges CHAIN and which root a child lands on depends on what
    // happened before it.
    std::vector<Proposal> apply;
    apply.reserve(best.size());
    for (const auto& kv : best) apply.push_back(kv.second);
    std::sort(apply.begin(), apply.end(), [](const Proposal& a, const Proposal& b) {
        if (a.child_count != b.child_count) return a.child_count < b.child_count;
        return a.child > b.child;
    });

    // Running read counts, and the union-find, keyed rather than indexed. Both are O(barcodes
    // actually corrected): every key either of them touches appears in some proposal, which is
    // what makes seeding them from the proposals complete.
    std::unordered_map<uint64_t, uint32_t> live, own;
    std::unordered_map<uint64_t, uint64_t> parent_of;
    live.reserve(apply.size() * 4);
    for (const Proposal& p : apply) {
        live.emplace(p.child, p.child_count);
        live.emplace(p.parent, p.parent_count);
        own.emplace(p.child, p.child_count);
    }
    for (const Proposal& p : apply) {
        if (p.parent_count <= p.child_count) ++res.merged_by_payload;
        uint64_t root = p.parent;
        for (int guard = 0; guard < 64; ++guard) {
            auto it = parent_of.find(root);
            if (it == parent_of.end()) break;
            root = it->second;
        }
        if (root == p.child) continue;  // never make a cycle
        parent_of[p.child] = root;
        live[root] += live[p.child];
        // This barcode's OWN reads, not its running total: merges chain, and a descendant's reads
        // were counted when the descendant moved.
        res.merged_reads += p.child_count;
        live[p.child] = 0;
        ++res.merged;
    }
    for (const auto& kv : parent_of) {
        uint64_t r = kv.second;
        for (int guard = 0; guard < 64; ++guard) {
            auto it = parent_of.find(r);
            if (it == parent_of.end()) break;
            r = it->second;
        }
        res.merges.push_back(CorrectionResult::Merge{kv.first, r, own[kv.first]});
    }
    // Key order, so the answer does not depend on a hash table's iteration order.
    std::sort(res.merges.begin(), res.merges.end(),
              [](const CorrectionResult::Merge& a, const CorrectionResult::Merge& b) {
                  return a.child < b.child;
              });

    // Best effort: the rotated copy is a temporary of this call, and failing to remove it is not
    // a reason to lose a correction result that is already computed.
    std::error_code rm_ec;
    std::filesystem::remove_all(rot_dir, rm_ec);
    std::filesystem::remove_all(third_dir, rm_ec);

    // The per-bucket molecule counts mean nothing on their own -- a bucket is a slice of the
    // barcode space, not of the molecules -- so they are recomputed here from the totals. Every
    // merge zeroes exactly one barcode, which is what makes this a subtraction rather than a scan.
    res.estimated_error = sub.sequencing_error;
    res.molecules_observed = static_cast<size_t>(n) - res.merged;
    const double space = sub.library_space;
    const double m_obs = static_cast<double>(res.molecules_observed);
    res.saturated = space > 0.0 && n > 0.05 * space;
    if (space > 0.0 && m_obs < 0.9 * space) {
        res.molecules_corrected = space * -std::log1p(-m_obs / space);
    } else {
        res.molecules_corrected = m_obs;
        res.saturated = true;
    }
    return res;
}

}  // namespace

CorrectionResult correct_umis(const UmiCounts& counts, const CorrectionParams& params,
                              const BarcodeEvidence& evidence) {
    if (params.cell_length < 0 || params.cell_length > counts.length())
        throw MigecError("correct_umis: cell length is outside the composite barcode");
    if (counts.spilled()) {
        if (evidence.has_quality() || evidence.has_payload()) {
            throw MigecError(
                "correct_umis: a side BarcodeEvidence is indexed against entries(), which a "
                "spilled counter does not have. Carry the evidence with the counter "
                "(UmiCounts::carry_evidence) so it is partitioned with it, or drop it "
                "deliberately -- ignoring it here would report a merge count from a weaker model "
                "as if it came from the full one.");
        }
        return correct_spilled(counts, params);
    }
    CorrectionParams p = params;
    const auto& entries = counts.entries();
    const BarcodeSearch search(entries, counts.length());
    const bool fit_error = p.library_distinct == 0 && p.sequencing_error < 0.0;
    UmiComposition comp;
    if (p.library_distinct == 0) {
        comp = counts.composition(false);
        p.library_distinct = counts.distinct();
        p.library_collision = collision_of(comp, counts.length());
        p.library_space = comp.effective_space();
    }
    auto run = [&](const BarcodeEvidence& ev) {
        const bool full_quality = entries.size() > 1 && ev.has_quality() &&
            ev.position_error.size() == entries.size() * static_cast<size_t>(counts.length()) &&
            std::all_of(ev.position_error.begin(), ev.position_error.end(),
                        [](float rate) { return rate > 0.0f && std::isfinite(rate); });
        auto fit = [&](uint64_t d1) {
            const double estimate = solve_umi_error(d1, static_cast<double>(entries.size()),
                comp, counts.length(), size_hist_of(entries));
            return estimate <= 0.0 ? 1e-4 : estimate;
        };
        // With every mismatch quality observed, candidate likelihoods do not depend
        // on the fitted fallback rate: census and decisions share ONE native batch.
        // Missing-quality inputs first fit the census using the SAME owned index.
        if (fit_error && !full_quality)
            p.sequencing_error = fit(d1_from_search(entries, search, 0, counts.length(), params.threads));
        uint64_t d1 = 0;
        auto result = correct_entries(entries, counts.length(), p, ev, nullptr, nullptr, &search,
                                      fit_error && full_quality ? &d1 : nullptr);
        if (fit_error && full_quality) result.estimated_error = fit(d1);
        return result;
    };
    if (counts.carries_evidence()) {
        // A resident table that carries its own evidence is one bucket of a partition that
        // happened to fit. Same call, same evidence, so the answer does not depend on whether the
        // budget was reached.
        CorrectionResult out;
        counts.for_each_bucket([&](const std::vector<UmiCounts::Entry>& bucket,
                                   const BarcodeEvidence& ev) {
            if (ev.has_payload()) {
                PayloadReservoir reservoir(ev.payload_width);
                for (size_t i = 0; i < bucket.size(); ++i) {
                    reservoir.add(ev.payload.data() +
                                  i * static_cast<size_t>(ev.payload_width));
                }
                p.library_clonality = reservoir.clonality(p);
            }
            out = run(ev);
        });
        return out;
    }
    return run(evidence);
}

}  // namespace migec
