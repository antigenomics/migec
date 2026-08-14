#include "migec/consensus.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <unordered_map>

#include "migec/types.hpp"

namespace migec {

namespace {

// log(1-e) and log(e/3) depend only on the reported Phred, so they tabulate. The per-column work
// is then adds; only the posterior at the end needs a transcendental, three times per output base
// rather than once per read base.
struct LogTables {
    std::array<double, kMaxPhred + 1> match{}, mismatch{};
};

LogTables build_tables(const std::vector<float>& calibration);

// Built once. `assemble_group` runs per molecule, and rebuilding 61 log/log1p pairs for every
// group put ~122 transcendentals against the ~3 per emitted base that the posterior actually
// needs. Same reason pattern.cpp keeps a `nominal_table()`; the calibrated path still builds its
// own, because there the table depends on the data.
const LogTables& nominal_tables() {
    static const LogTables t = build_tables({});
    return t;
}

LogTables build_tables(const std::vector<float>& calibration) {
    LogTables t;
    for (int q = 0; q <= kMaxPhred; ++q) {
        double e = (static_cast<size_t>(q) < calibration.size())
                       ? static_cast<double>(calibration[static_cast<size_t>(q)])
                       : phred_error(static_cast<uint8_t>(q));
        // A calibration table can legitimately contain a zero (no mismatch was seen at that Q);
        // log(0) is not a likelihood, so bound it by the best the instrument ever claims.
        e = std::clamp(e, 1e-7, 0.75);
        t.match[static_cast<size_t>(q)] = std::log1p(-e);
        t.mismatch[static_cast<size_t>(q)] = std::log(e / 3.0);
    }
    return t;
}

int read_length(const ConsensusRead& r) {
    return static_cast<int>(std::min(r.seq.size(), r.qual.size()));
}

// The group's frame runs from the leftmost read start to the rightmost read end. With every
// offset at 0 -- the amplicon case -- that is the shortest read, which is the old behaviour: no
// column is emitted that some read does not reach.
int group_width(const std::vector<ConsensusRead>& reads) {
    if (reads.empty()) return 0;
    int begin = reads[0].offset, end = 0;
    for (const ConsensusRead& r : reads) {
        begin = std::min(begin, r.offset);
        end = std::max(end, r.offset + read_length(r));
    }
    bool all_anchored = true;
    for (const ConsensusRead& r : reads) all_anchored = all_anchored && r.offset == begin;
    if (all_anchored) {
        int w = read_length(reads[0]);
        for (const ConsensusRead& r : reads) w = std::min(w, read_length(r));
        return w;
    }
    return end - begin;
}

// The base a read contributes at frame column j, or 0 if it does not reach that far.
inline char base_at(const ConsensusRead& r, int j) {
    const int i = j - r.offset;
    return (i < 0 || i >= read_length(r)) ? '\0' : r.seq[static_cast<size_t>(i)];
}

inline char qual_at(const ConsensusRead& r, int j) {
    const int i = j - r.offset;
    return (i < 0 || i >= read_length(r)) ? '\0' : r.qual[static_cast<size_t>(i)];
}

// log C(n, k), for the hypergeometric tail. Not a hot path -- at most 28 pairs per group.
double log_choose(int n, int k) {
    if (k < 0 || k > n) return -std::numeric_limits<double>::infinity();
    return std::lgamma(n + 1.0) - std::lgamma(k + 1.0) - std::lgamma(n - k + 1.0);
}

// -log10 P(X >= c11) for X ~ Hypergeometric(n; a, b).
double hypergeom_sf(int c11, int n, int a, int b) {
    const double denom = log_choose(n, a);
    double acc = -std::numeric_limits<double>::infinity();
    for (int x = c11; x <= std::min(a, b); ++x) {
        const double term = log_choose(b, x) + log_choose(n - b, a - x) - denom;
        acc = (acc == -std::numeric_limits<double>::infinity())
                  ? term
                  : std::max(acc, term) + std::log1p(std::exp(-std::abs(acc - term)));
    }
    if (acc == -std::numeric_limits<double>::infinity()) return 0.0;
    return -acc / std::log(10.0);
}

struct MinorColumn {
    size_t position;
    char major;
    uint32_t minor_reads;
};

// Positions where a minority of the reads carry a different base, strongest first. A position
// where the "minor" allele is half the reads is still callable -- that is what a 50/50 doublet
// looks like -- but one carried by a single read is not, because a single read carries no
// linkage information at all.
std::vector<MinorColumn> minor_columns(const std::vector<ConsensusRead>& reads,
                                       const ConsensusParams& params) {
    const int width = group_width(reads);
    std::vector<MinorColumn> out;
    for (int j = 0; j < width; ++j) {
        std::array<uint32_t, 4> counts{};
        uint32_t called = 0;
        for (const ConsensusRead& r : reads) {
            const uint8_t c = base_code(base_at(r, j));
            if (c != kInvalidBase) { ++counts[c]; ++called; }
        }
        if (called == 0) continue;
        const size_t best = static_cast<size_t>(
            std::max_element(counts.begin(), counts.end()) - counts.begin());
        const uint32_t minor = called - counts[best];
        // Against the reads that actually COVER this column, not against the group: in a contig a
        // column reached by four reads must be judged on four.
        if (minor >= params.min_minor_reads && minor * 2 <= called) {
            out.push_back({static_cast<size_t>(j), base_char(static_cast<uint8_t>(best)), minor});
        }
    }
    std::sort(out.begin(), out.end(), [](const MinorColumn& a, const MinorColumn& b) {
        if (a.minor_reads != b.minor_reads) return a.minor_reads > b.minor_reads;
        return a.position < b.position;  // ties by position, so the result is reproducible
    });
    if (out.size() > static_cast<size_t>(params.max_split_columns)) {
        out.resize(static_cast<size_t>(params.max_split_columns));
    }
    return out;
}

// The reads x columns "carries the minor base here" matrix, as one bitmask per read.
std::vector<uint32_t> minor_matrix(const std::vector<ConsensusRead>& reads,
                                   const std::vector<MinorColumn>& cols) {
    std::vector<uint32_t> rows(reads.size(), 0);
    for (size_t i = 0; i < reads.size(); ++i) {
        for (size_t c = 0; c < cols.size(); ++c) {
            const char b = base_at(reads[i], static_cast<int>(cols[c].position));
            if (base_code(b) != kInvalidBase && b != cols[c].major) rows[i] |= 1u << c;
        }
    }
    return rows;
}

struct BestPair {
    double score = 0.0;
    int x = -1, y = -1;
    bool flip = false;  // the subclone carries the minor at x and the MAJOR at y
};

// Two-sided, because at a 50/50 split which allele is "major" is a coin toss taken separately per
// column: two columns of a genuine doublet then come out anti-correlated, and a one-sided test
// scores the strongest possible evidence as nothing at all. The excess of minor-with-MAJOR is the
// same evidence in the opposite phase, so both are tested and the count is halved for it.
BestPair best_linked_pair(const std::vector<uint32_t>& rows, size_t ncols) {
    BestPair best;
    if (ncols < 2) return best;
    const int n = static_cast<int>(rows.size());
    const double correction =
        std::log10(static_cast<double>(ncols * (ncols - 1) / 2)) + std::log10(2.0);
    std::vector<int> colsum(ncols, 0);
    for (uint32_t r : rows) {
        for (size_t c = 0; c < ncols; ++c) colsum[c] += (r >> c) & 1u;
    }
    for (size_t x = 0; x < ncols; ++x) {
        for (size_t y = x + 1; y < ncols; ++y) {
            int c11 = 0;
            for (uint32_t r : rows) c11 += ((r >> x) & 1u) && ((r >> y) & 1u);
            const int a = colsum[x], b = colsum[y];
            const double same = c11 >= 2 ? hypergeom_sf(c11, n, a, b) : 0.0;
            const double opposite =
                (a - c11) >= 2 ? hypergeom_sf(a - c11, n, a, n - b) : 0.0;
            const bool flip = opposite > same;
            const double s = std::max(same, opposite) - correction;
            if (s > best.score) best = {s, static_cast<int>(x), static_cast<int>(y), flip};
        }
    }
    return best;
}

Consensus call_consensus(const std::vector<ConsensusRead>& reads, const ConsensusParams& params,
                         const LogTables& tables) {
    const int width = group_width(reads);
    Consensus out;
    out.reads = static_cast<uint32_t>(reads.size());
    out.seq.resize(static_cast<size_t>(width));
    out.qual.resize(static_cast<size_t>(width));
    const uint8_t q_cap = static_cast<uint8_t>(
        std::min<double>(kMaxPhred, std::floor(-10.0 * std::log10(params.rt_floor))));
    double err_sum = 0.0;
    for (int j = 0; j < width; ++j) {
        std::array<double, 4> ll{0.0, 0.0, 0.0, 0.0};
        bool any = false;
        for (const ConsensusRead& r : reads) {
            const uint8_t c = base_code(base_at(r, j));
            if (c == kInvalidBase) continue;
            any = true;
            const size_t q = phred_from_char(qual_at(r, j));
            for (int b = 0; b < 4; ++b) {
                ll[static_cast<size_t>(b)] +=
                    (b == c) ? tables.match[q] : tables.mismatch[q];
            }
        }
        if (!any) {
            out.seq[static_cast<size_t>(j)] = 'N';
            out.qual[static_cast<size_t>(j)] = char_from_phred(0);
            continue;
        }
        // A tie is resolved by base order rather than by an N: the posterior is then 0.5 and the
        // emitted quality says so (~Q3). An N discards the information that it is one of two.
        const size_t best = static_cast<size_t>(
            std::max_element(ll.begin(), ll.end()) - ll.begin());
        double rest = 0.0;
        for (size_t b = 0; b < 4; ++b) {
            if (b != best) rest += std::exp(ll[b] - ll[best]);
        }
        const double p_err = rest / (1.0 + rest);
        err_sum += p_err;
        // The floor is added, not compared: an RT error is in every read, so the two failure
        // modes are independent and the emitted quality must carry both.
        const double total = p_err + params.rt_floor;
        int q = static_cast<int>(std::lround(-10.0 * std::log10(total)));
        q = std::clamp(q, static_cast<int>(params.min_quality), static_cast<int>(q_cap));
        out.seq[static_cast<size_t>(j)] = base_char(static_cast<uint8_t>(best));
        out.qual[static_cast<size_t>(j)] = char_from_phred(static_cast<uint8_t>(q));
    }
    out.mean_error = width ? err_sum / static_cast<double>(width) : 0.0;
    return out;
}

// Counting mode: the group's most frequent exact sequence, with each base carrying the best
// quality any read of that sequence reported for it.
//
// Never: the max is taken over the reads that carry the WINNING sequence, not over every read in
// the group. Maximising across variants would take its highest quality from exactly the reads that
// disagree at that position, which asserts most confidence where the evidence conflicts.
//
// The RT floor is added as it is in the full path, so this path claims no more than that one.
// What is missing by construction is error correction: a base is right because one read read it
// well, not because n reads agreed on it.
Consensus modal_consensus(const std::vector<ConsensusRead>& reads, const ConsensusParams& params) {
    Consensus out;
    out.reads = static_cast<uint32_t>(reads.size());

    // Ties by the sequence itself, so the output does not depend on read order within the group.
    std::unordered_map<std::string_view, uint32_t> votes;
    votes.reserve(reads.size() * 2);
    for (const ConsensusRead& r : reads) ++votes[r.seq];
    std::string_view best;
    uint32_t best_votes = 0;
    for (const auto& [seq, n] : votes) {
        if (n > best_votes || (n == best_votes && seq < best)) { best_votes = n; best = seq; }
    }
    out.support = best_votes;

    const size_t width = best.size();
    std::vector<uint8_t> quality(width, 0);
    for (const ConsensusRead& r : reads) {
        if (r.seq != best) continue;
        const size_t n = std::min(width, r.qual.size());
        for (size_t j = 0; j < n; ++j) {
            quality[j] = std::max(quality[j], phred_from_char(r.qual[j]));
        }
    }

    const uint8_t q_cap = static_cast<uint8_t>(
        std::min<double>(kMaxPhred, std::floor(-10.0 * std::log10(params.rt_floor))));
    out.seq.assign(best);
    out.qual.resize(width);
    double err_sum = 0.0;
    for (size_t j = 0; j < width; ++j) {
        const double p = phred_error(quality[j]);
        err_sum += p;
        int q = static_cast<int>(std::lround(-10.0 * std::log10(p + params.rt_floor)));
        q = std::clamp(q, static_cast<int>(params.min_quality), static_cast<int>(q_cap));
        out.qual[j] = char_from_phred(static_cast<uint8_t>(q));
    }
    out.mean_error = width ? err_sum / static_cast<double>(width) : 0.0;
    return out;
}

// Union-find carrying an offset relative to the component root, so a component's reads come out
// already placed against each other.
struct OffsetUnion {
    std::vector<int> parent;
    std::vector<int> delta;  // offset of this node relative to its parent

    explicit OffsetUnion(size_t n) : parent(n), delta(n, 0) {
        for (size_t i = 0; i < n; ++i) parent[i] = static_cast<int>(i);
    }

    // Returns the root and sets `offset` to this node's offset within the component.
    int find(int i, int& offset) {
        if (parent[static_cast<size_t>(i)] == i) { offset = 0; return i; }
        int up = 0;
        const int root = find(parent[static_cast<size_t>(i)], up);
        delta[static_cast<size_t>(i)] += up;      // path compression keeps the offset consistent
        parent[static_cast<size_t>(i)] = root;
        offset = delta[static_cast<size_t>(i)];
        return root;
    }

    // Join so that offset(a) - offset(b) == shift. Silently declines a join that contradicts one
    // already made: a seed hit that disagrees with the layout is evidence against itself, and
    // forcing it in would shift every read that was already placed.
    void join(int a, int b, int shift) {
        int oa = 0, ob = 0;
        const int ra = find(a, oa), rb = find(b, ob);
        if (ra == rb) return;
        parent[static_cast<size_t>(rb)] = ra;
        delta[static_cast<size_t>(rb)] = oa - ob - shift;
    }
};

// The modal offset between two reads, by exact seed votes. Returns false when the evidence is too
// thin -- no offset at all is the right answer for two reads of the same molecule that simply do
// not overlap, and inventing one is what glues a contig across a gap.
bool seed_offset(const ConsensusRead& a, const ConsensusRead& b, const ConsensusParams& params,
                 int& shift) {
    const int k = params.seed_length;
    const int la = read_length(a), lb = read_length(b);
    if (la < k || lb < k) return false;
    // ponytail: a plain map of seed -> first position, and O(n^2) pairs per group. Contig mode
    // runs on random-primed data where a barcode carries a handful of reads (X1: 1.5% of 10x
    // groups hold more than one), so the quadratic is over single digits. Index the group once if
    // a benchmark ever puts hundreds of reads on one barcode.
    std::unordered_map<std::string_view, int> seeds;
    seeds.reserve(static_cast<size_t>(la - k + 1));
    for (int i = 0; i <= la - k; ++i) {
        seeds.emplace(a.seq.substr(static_cast<size_t>(i), static_cast<size_t>(k)), i);
    }
    std::unordered_map<int, int> votes;
    for (int j = 0; j <= lb - k; ++j) {
        auto it = seeds.find(b.seq.substr(static_cast<size_t>(j), static_cast<size_t>(k)));
        if (it != seeds.end()) ++votes[it->second - j];
    }
    int best_votes = 0;
    for (const auto& [offset, n] : votes) {
        if (n > best_votes || (n == best_votes && offset < shift)) { best_votes = n; shift = offset; }
    }
    if (best_votes < params.min_seed_votes) return false;
    // The implied overlap has to be long enough to be evidence rather than a repeat.
    const int overlap = std::min(shift + la, lb) - std::max(shift, 0);
    return overlap >= params.min_overlap;
}

}  // namespace

std::vector<std::vector<ConsensusRead>> place_reads(const std::vector<ConsensusRead>& reads,
                                                    const ConsensusParams& params) {
    OffsetUnion uf(reads.size());
    for (size_t i = 0; i < reads.size(); ++i) {
        for (size_t j = i + 1; j < reads.size(); ++j) {
            int shift = 0;
            // seed_offset returns offset(j) - offset(i), so join(j, i, shift).
            if (seed_offset(reads[i], reads[j], params, shift)) {
                uf.join(static_cast<int>(j), static_cast<int>(i), shift);
            }
        }
    }
    std::unordered_map<int, std::vector<ConsensusRead>> components;
    for (size_t i = 0; i < reads.size(); ++i) {
        int offset = 0;
        const int root = uf.find(static_cast<int>(i), offset);
        ConsensusRead placed = reads[i];
        placed.offset = offset;
        components[root].push_back(placed);
    }
    std::vector<std::vector<ConsensusRead>> out;
    out.reserve(components.size());
    for (auto& [root, group] : components) {
        int begin = group[0].offset;
        for (const ConsensusRead& r : group) begin = std::min(begin, r.offset);
        for (ConsensusRead& r : group) r.offset -= begin;
        out.push_back(std::move(group));
    }
    // Largest first, and within a size by the component's smallest read. The tiebreak has to be a
    // property of the SET, not of its first element: components come out of a hash map, and the
    // order of reads inside one follows the input, so `a[0]` makes the result depend on both.
    auto smallest = [](const std::vector<ConsensusRead>& g) {
        std::string_view m = g[0].seq;
        for (const ConsensusRead& r : g) m = std::min(m, r.seq);
        return m;
    };
    std::sort(out.begin(), out.end(),
              [&smallest](const std::vector<ConsensusRead>& a, const std::vector<ConsensusRead>& b) {
                  if (a.size() != b.size()) return a.size() > b.size();
                  return smallest(a) < smallest(b);
              });
    return out;
}

double linkage_score(const std::vector<ConsensusRead>& reads, const ConsensusParams& params) {
    if (reads.size() < params.min_split_reads) return 0.0;
    const std::vector<MinorColumn> cols = minor_columns(reads, params);
    if (cols.size() < 2) return 0.0;
    return best_linked_pair(minor_matrix(reads, cols), cols.size()).score;
}

namespace {

// One overlap component -> one or two molecules. Split before consensus, never after.
std::vector<Consensus> assemble_component(const std::vector<ConsensusRead>& reads,
                                          const ConsensusParams& params,
                                          const LogTables& tables) {
    std::vector<Consensus> out;

    double score = 0.0;
    if (reads.size() >= params.min_split_reads) {
        const std::vector<MinorColumn> cols = minor_columns(reads, params);
        if (cols.size() >= 2) {
            const std::vector<uint32_t> rows = minor_matrix(reads, cols);
            const BestPair pair = best_linked_pair(rows, cols.size());
            score = pair.score;
            if (pair.score > params.linkage_threshold) {
                // ponytail: one split level. A MIG holding three molecules is rare enough that
                // recursing is not worth the extra false-positive surface; recurse here if a
                // benchmark ever shows otherwise.
                std::vector<ConsensusRead> minor, major;
                const uint32_t mask = (1u << pair.x) | (1u << pair.y);
                const uint32_t want = pair.flip ? (1u << pair.x) : mask;
                for (size_t i = 0; i < reads.size(); ++i) {
                    ((rows[i] & mask) == want ? minor : major).push_back(reads[i]);
                }
                if (minor.size() >= params.min_subclone_reads &&
                    major.size() >= params.min_subclone_reads) {
                    for (const std::vector<ConsensusRead>* part : {&major, &minor}) {
                        Consensus c = call_consensus(*part, params, tables);
                        c.linkage = pair.score;
                        out.push_back(std::move(c));
                    }
                    return out;
                }
            }
        }
    }
    Consensus c = call_consensus(reads, params, tables);
    c.linkage = score;
    out.push_back(std::move(c));
    return out;
}

}  // namespace

std::vector<Consensus> assemble_group(const std::vector<ConsensusRead>& reads,
                                      const ConsensusParams& params) {
    std::vector<Consensus> out;
    if (reads.empty()) return out;
    // Before the tables are touched: the fast path computes no likelihood at all.
    if (params.fast) {
        out.push_back(modal_consensus(reads, params));
        return out;
    }
    LogTables calibrated;
    if (!params.calibration.empty()) calibrated = build_tables(params.calibration);
    const LogTables& tables = params.calibration.empty() ? nominal_tables() : calibrated;

    if (!params.contig) {
        return assemble_component(reads, params, tables);
    }
    // Random-primed reads tile the molecule, so they are placed against each other and cut into
    // overlap components first. A component is a contig; the group's components are never merged,
    // because a barcode's reads failing to reach each other is exactly the case where no read
    // covers the sequence a single consensus would be asserting.
    const std::vector<std::vector<ConsensusRead>> components = place_reads(reads, params);
    for (size_t c = 0; c < components.size(); ++c) {
        for (Consensus& m : assemble_component(components[c], params, tables)) {
            m.component = static_cast<uint32_t>(c);
            m.components = static_cast<uint32_t>(components.size());
            out.push_back(std::move(m));
        }
    }
    return out;
}

std::vector<Consensus> assemble_pairs(const std::vector<ConsensusRead>& mates1,
                                      const std::vector<ConsensusRead>& mates2,
                                      const ConsensusParams& params) {
    if (mates2.empty()) return assemble_group(mates1, params);

    // Pairs enough to outvote a pair whose overlap happened to fall in an error. Eight, because
    // the offset has one degree of freedom and a wrong vote needs `min_seed_votes` exact seeds to
    // even be cast -- this is a tie-break, not an estimate.
    constexpr size_t kOffsetVoteSample = 8;
    std::unordered_map<int, int> votes;
    const size_t sample = std::min(mates2.size(), kOffsetVoteSample);
    for (size_t i = 0; i < sample; ++i) {
        if (i >= mates1.size() || mates2[i].seq.empty()) continue;
        int shift = 0;
        // Two pairs agreeing is the answer: a vote is only cast when `min_seed_votes` exact seeds
        // already agree within the pair, so two independent pairs saying the same thing settles a
        // one-degree-of-freedom layout. Stopping here is most of what makes merging cheap.
        if (seed_offset(mates1[i], mates2[i], params, shift) && ++votes[shift] >= 2) break;
    }
    int best = 0, best_votes = 0;
    for (const auto& [offset, n] : votes) {
        // Ties go to the smaller offset so the answer is the group's, never the iteration order's.
        if (n > best_votes || (n == best_votes && offset < best)) { best_votes = n; best = offset; }
    }

    LogTables calibrated;
    if (!params.calibration.empty()) calibrated = build_tables(params.calibration);
    const LogTables& tables = params.calibration.empty() ? nominal_tables() : calibrated;

    if (best_votes == 0) {
        // The mates do not reach each other. Two contigs, and the bases between them are asserted
        // by nobody -- which is the same rule `place_reads` follows and the reason this is
        // placement rather than a merge.
        std::vector<ConsensusRead> second;
        second.reserve(mates2.size());
        for (const ConsensusRead& r : mates2) {
            if (!r.seq.empty()) second.push_back({r.seq, r.qual, 0});
        }
        std::vector<Consensus> out;
        const uint32_t components = second.empty() ? 1u : 2u;
        for (Consensus& m : assemble_component(mates1, params, tables)) {
            m.component = 0;
            m.components = components;
            out.push_back(std::move(m));
        }
        for (Consensus& m : assemble_component(second, params, tables)) {
            m.component = 1;
            m.components = components;
            out.push_back(std::move(m));
        }
        return out;
    }

    // One frame for the whole group: mate 1 where checkout anchored it, mate 2 at the voted
    // offset. A negative offset means mate 2 starts first, so everything shifts to keep the frame
    // non-negative -- `base_at` indexes from zero.
    const int base1 = best < 0 ? -best : 0;
    const int base2 = best < 0 ? 0 : best;
    std::vector<ConsensusRead> placed;
    placed.reserve(mates1.size() + mates2.size());
    for (const ConsensusRead& r : mates1) placed.push_back({r.seq, r.qual, base1});
    for (const ConsensusRead& r : mates2) {
        if (!r.seq.empty()) placed.push_back({r.seq, r.qual, base2});
    }
    return assemble_component(placed, params, tables);
}

}  // namespace migec
