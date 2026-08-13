#include "migec/suggest.hpp"

#include <algorithm>
#include <cmath>

#include "migec/fastq.hpp"
#include "migec/pattern.hpp"
#include "migec/types.hpp"

namespace migec {

std::array<double, 4> CycleStats::frequencies() const {
    std::array<double, 4> f{};
    const double tot = static_cast<double>(counts[0] + counts[1] + counts[2] + counts[3]);
    if (tot <= 0.0) return f;
    for (int i = 0; i < 4; ++i) f[static_cast<size_t>(i)] = counts[static_cast<size_t>(i)] / tot;
    return f;
}

double CycleStats::entropy() const {
    double h = 0.0;
    for (double p : frequencies()) {
        if (p > 0.0) h -= p * std::log2(p);
    }
    return h;
}

double CycleStats::collision() const {
    double m = 0.0;
    for (double p : frequencies()) m += p * p;
    return m;
}

char CycleStats::consensus() const {
    const std::array<double, 4> f = frequencies();
    int best = 0;
    for (int i = 1; i < 4; ++i) {
        if (f[static_cast<size_t>(i)] > f[static_cast<size_t>(best)]) best = i;
    }
    return base_char(static_cast<uint8_t>(best));
}

double CycleStats::consensus_fraction() const {
    const std::array<double, 4> f = frequencies();
    return *std::max_element(f.begin(), f.end());
}

double CycleStats::mean_phred() const {
    return n ? static_cast<double>(phred_sum) / static_cast<double>(n) : 0.0;
}

double CycleStats::deviation_from_uniform() const {
    double d = 0.0;
    for (double p : frequencies()) d += std::fabs(p - 0.25);
    return d;  // 0 when flat, 1.5 when fixed -- reported halved below for a 0..0.75 scale
}

namespace {

// The k-mer counters: exact counts in a flat array, plus the position sum so a hit can say where
// in the read it sits. An adapter is at one end; a repeat is everywhere.
struct KmerCounts {
    static constexpr size_t kSize = size_t{1} << (2 * kKmerLength);
    std::vector<uint64_t> count = std::vector<uint64_t>(kSize, 0);
    std::vector<uint64_t> position_sum = std::vector<uint64_t>(kSize, 0);
    std::vector<uint32_t> last_read = std::vector<uint32_t>(kSize, 0);  // reads carrying it once+
    std::vector<uint64_t> reads_with = std::vector<uint64_t>(kSize, 0);
    uint64_t scanned = 0;

    // One pass, 2 bits per base, restarting the window at anything that is not ACGT -- an N is not
    // a base and a k-mer spanning one is not evidence of anything.
    void add(std::string_view seq, uint32_t read_index) {
        constexpr uint64_t mask = (uint64_t{1} << (2 * kKmerLength)) - 1;
        uint64_t code = 0;
        int run = 0;
        for (size_t i = 0; i < seq.size(); ++i) {
            const uint8_t b = base_code(seq[i]);
            if (b == kInvalidBase) { run = 0; code = 0; continue; }
            code = ((code << 2) | b) & mask;
            if (++run < kKmerLength) continue;
            const size_t at = static_cast<size_t>(code);
            ++count[at];
            position_sum[at] += i + 1 - kKmerLength;
            ++scanned;
            // read_index + 1, so that 0 means "not seen yet" without a second pass to clear.
            if (last_read[at] != read_index + 1) {
                last_read[at] = read_index + 1;
                ++reads_with[at];
            }
        }
    }
};

std::string kmer_string(uint64_t code) {
    std::string s(static_cast<size_t>(kKmerLength), 'A');
    for (int i = kKmerLength - 1; i >= 0; --i) {
        s[static_cast<size_t>(i)] = base_char(static_cast<uint8_t>(code & 3u));
        code >>= 2;
    }
    return s;
}

}  // namespace

CycleProfile profile_cycles(const std::string& fastq_path, int n_cycles, uint64_t max_reads) {
    if (n_cycles <= 0) throw MigecError("suggest: n_cycles must be positive");
    CycleProfile p;
    p.cycles.resize(static_cast<size_t>(n_cycles));
    for (int j = 0; j < n_cycles; ++j) p.cycles[static_cast<size_t>(j)].cycle = j;
    KmerCounts kmers;

    FastqReader r(fastq_path);
    FastqRecord rec;
    while ((max_reads == 0 || p.reads < max_reads) && r.next(rec)) {
        p.read_length = std::max(p.read_length, rec.seq.size());
        // Over the WHOLE read, not the profiled prefix: an adapter that survived trimming is at
        // the 3' end, which is exactly the part the cycle profile does not reach.
        kmers.add(rec.seq, static_cast<uint32_t>(p.reads));
        for (char ch : rec.seq) {
            const uint8_t b = base_code(ch);
            if (b != kInvalidBase) ++p.base_counts[b];
        }
        const size_t lim = std::min(rec.seq.size(), static_cast<size_t>(n_cycles));
        for (size_t j = 0; j < lim; ++j) {
            CycleStats& c = p.cycles[j];
            const uint8_t code = base_code(rec.seq[j]);
            ++c.counts[code == kInvalidBase ? 4 : code];
            if (j < rec.qual.size()) {
                c.phred_sum += phred_from_char(rec.qual[j]);
                ++c.n;
            }
        }
        ++p.reads;
    }
    if (p.reads == 0) throw MigecError("suggest: no reads in " + fastq_path);

    // Overrepresentation is against the reads' OWN base composition, never against a flat 1/4: a
    // 70% AT library makes every AT-rich k-mer look enriched fivefold against uniform, and the
    // report would then be a list of the library's GC content.
    p.kmers_scanned = kmers.scanned;
    const double total = static_cast<double>(
        p.base_counts[0] + p.base_counts[1] + p.base_counts[2] + p.base_counts[3]);
    std::array<double, 4> freq{0.25, 0.25, 0.25, 0.25};
    if (total > 0.0) {
        for (int b = 0; b < 4; ++b) {
            freq[static_cast<size_t>(b)] = static_cast<double>(p.base_counts[static_cast<size_t>(b)]) / total;
        }
    }
    // Seen fewer times than this and the ratio is a small-count artefact rather than a finding.
    const uint64_t min_count = std::max<uint64_t>(20, kmers.scanned / 100000);
    for (size_t code = 0; code < KmerCounts::kSize; ++code) {
        const uint64_t n = kmers.count[code];
        if (n < min_count) continue;
        double q = 1.0;
        for (int i = 0; i < kKmerLength; ++i) {
            q *= freq[static_cast<size_t>((code >> (2 * i)) & 3u)];
        }
        const double expected = q * static_cast<double>(kmers.scanned);
        if (expected <= 0.0) continue;
        KmerHit hit;
        hit.kmer = kmer_string(code);
        hit.count = n;
        hit.expected = expected;
        hit.ratio = static_cast<double>(n) / expected;
        hit.mean_position = static_cast<double>(kmers.position_sum[code]) / static_cast<double>(n);
        hit.read_fraction =
            static_cast<double>(kmers.reads_with[code]) / static_cast<double>(p.reads);
        p.kmers.push_back(std::move(hit));
    }
    std::sort(p.kmers.begin(), p.kmers.end(), [](const KmerHit& a, const KmerHit& b) {
        if (a.ratio != b.ratio) return a.ratio > b.ratio;
        return a.kmer < b.kmer;  // ties by sequence, so the table is reproducible
    });
    if (p.kmers.size() > 64) p.kmers.resize(64);
    return p;
}

Suggestion suggest_pattern(const CycleProfile& profile, double umi_deviation,
                           double constant_fraction, int min_run) {
    Suggestion out;
    out.profile = profile;
    const size_t n = profile.cycles.size();
    if (n == 0) return out;

    // Classify each cycle on its own, then smooth by requiring a run: one stray cycle in the
    // middle of a UMI block is sampling noise, not a fixed base, and a pattern that alternates
    // N and a literal every other position is a description of the noise.
    std::vector<CycleKind> kind(n, CycleKind::kVariable);
    for (size_t j = 0; j < n; ++j) {
        const CycleStats& c = profile.cycles[j];
        const double dev = c.deviation_from_uniform() / 2.0;  // 0 flat .. 0.75 fixed
        if (c.consensus_fraction() >= constant_fraction) {
            kind[j] = CycleKind::kConstant;
        } else if (dev <= umi_deviation) {
            kind[j] = CycleKind::kUmi;
        }
    }
    // A UMI run shorter than `min_run` is almost certainly a stretch of biological sequence that
    // happens to be even. Demote it.
    for (size_t j = 0; j < n;) {
        size_t k = j;
        while (k < n && kind[k] == kind[j]) ++k;
        if (kind[j] == CycleKind::kUmi && static_cast<int>(k - j) < min_run) {
            for (size_t t = j; t < k; ++t) kind[t] = CycleKind::kVariable;
        }
        j = k;
    }

    for (size_t j = 0; j < n;) {
        size_t k = j;
        while (k < n && kind[k] == kind[j]) ++k;
        Segment s;
        s.kind = kind[j];
        s.begin = static_cast<int>(j);
        s.end = static_cast<int>(k);
        double dev = 0.0;
        for (size_t t = j; t < k; ++t) {
            dev += profile.cycles[t].deviation_from_uniform() / 2.0;
            if (s.kind == CycleKind::kConstant) s.consensus.push_back(profile.cycles[t].consensus());
        }
        s.mean_deviation = dev / static_cast<double>(k - j);
        out.segments.push_back(s);
        j = k;
    }

    // The pattern ends at the last CONSTANT segment, not at the last non-payload one.
    //
    // Composition cannot tell a UMI from diverse payload: both show four lines near 1/4. What
    // separates them is that a barcode is anchored and payload is not -- a uniform run with no
    // constant sequence after it cannot be placed in a read, so claiming it as a barcode would
    // produce a pattern that matches everywhere. Anything past the final anchor is reported in
    // the note and left out of the pattern.
    size_t last = 0;
    for (size_t i = 0; i < out.segments.size(); ++i) {
        if (out.segments[i].kind == CycleKind::kConstant) last = i + 1;
    }
    int trailing_uniform = 0;
    for (size_t i = last; i < out.segments.size(); ++i) {
        if (out.segments[i].kind == CycleKind::kUmi) trailing_uniform += out.segments[i].length();
    }
    for (size_t i = 0; i < last; ++i) {
        const Segment& s = out.segments[i];
        if (s.kind == CycleKind::kUmi) {
            out.pattern.append(static_cast<size_t>(s.length()), 'N');
            out.umi_length += s.length();
        } else if (s.kind == CycleKind::kConstant) {
            // Lowercase: scored at half weight. A constant run recovered from the reads is an
            // adapter or primer, which is where mismatches are expected and tolerated.
            for (char c : s.consensus) out.pattern.push_back(static_cast<char>(std::tolower(c)));
            out.anchor_length += s.length();
        } else {
            // A variable stretch inside the pattern is real but unmodelled -- a sample tag we
            // cannot read off one file, most often. '.' matches anything and scores nothing.
            out.pattern.append(static_cast<size_t>(s.length()), '.');
        }
    }

    if (out.umi_length == 0 && trailing_uniform > 0) {
        out.note = "the only near-uniform run sits after the last constant sequence, with nothing "
                   "to anchor it. That is what diverse payload looks like. If the barcode really "
                   "is at the 3' end, raise --cycles so the primer past it is profiled too.";
    } else if (out.umi_length == 0) {
        out.note = "no UMI block found: no run of cycles is close enough to a flat 1/4 trace. "
                   "If the library is paired, the UMI is probably on the other mate.";
    } else if (trailing_uniform >= min_run) {
        out.note = std::to_string(trailing_uniform) +
                   " cycles past the anchor are also near-uniform. Payload across many molecules "
                   "looks like this, so they are excluded; if they are a second barcode, raise "
                   "--cycles to profile the sequence that anchors them.";
    } else if (out.anchor_length == 0) {
        out.note = "a UMI block with nothing constant after it, so there is nothing to anchor on. "
                   "Extend --cycles, or the barcode may sit at the 3' end.";
    } else if (out.umi_length > kMaxBarcodeLen) {
        out.note = "the UMI block is longer than the 32 bases the packed representation holds; "
                   "this is more likely two barcodes than one.";
    } else if (out.anchor_length > 30) {
        // A conserved locus is indistinguishable from an adapter by composition alone -- both are
        // one base at ~100% per cycle. It matters because --trim pattern discards everything the
        // pattern covers, and discarding conserved biology is not the same as discarding a primer.
        out.note = "the constant run is " + std::to_string(out.anchor_length) +
                   " nt. Composition cannot tell a primer from a conserved locus, and `--trim "
                   "pattern` discards everything the pattern covers -- check where the primer "
                   "actually ends before trimming that much.";
    }
    return out;
}

}  // namespace migec
