#include "migec/checkout.hpp"

#include <zlib.h>

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>

#include "migec/fastq.hpp"
#include "migec/mig_record.hpp"
#include "migec/parallel.hpp"
#include "migec/resource.hpp"
#include "migec/types.hpp"

namespace migec {

void QualityCalibration::merge(const QualityCalibration& o) {
    for (size_t p = 0; p < by_position.size() && p < o.by_position.size(); ++p) {
        for (size_t q = 0; q < 61; ++q) {
            by_position[p][q][0] += o.by_position[p][q][0];
            by_position[p][q][1] += o.by_position[p][q][1];
        }
    }
}

void QualityCalibration::fit(uint64_t min_bases, double max_excess) {
    // Per-position mismatch rate first. A position far above the median is not measuring the
    // instrument -- it is a base the pattern is wrong about, or one that is genuinely variable,
    // and leaving it in would put its variation into the intercept as a fake error floor.
    std::vector<double> rates;
    std::vector<uint64_t> totals(by_position.size(), 0), bad(by_position.size(), 0);
    for (size_t p = 0; p < by_position.size(); ++p) {
        for (size_t q = 0; q < 61; ++q) {
            totals[p] += by_position[p][q][0] + by_position[p][q][1];
            bad[p] += by_position[p][q][1];
        }
        if (totals[p] >= min_bases) {
            rates.push_back(static_cast<double>(bad[p]) / static_cast<double>(totals[p]));
        }
    }
    double median = 0.0;
    if (!rates.empty()) {
        std::sort(rates.begin(), rates.end());
        median = rates[rates.size() / 2];
    }
    counts = {};
    position_used.assign(by_position.size(), 0);
    positions_dropped = 0;
    for (size_t p = 0; p < by_position.size(); ++p) {
        if (!totals[p]) continue;
        const double rate = static_cast<double>(bad[p]) / static_cast<double>(totals[p]);
        if (totals[p] >= min_bases && median > 0.0 && rate > max_excess * median) {
            ++positions_dropped;
            continue;
        }
        position_used[p] = 1;
        for (size_t q = 0; q < 61; ++q) {
            counts[q][0] += by_position[p][q][0];
            counts[q][1] += by_position[p][q][1];
        }
    }

    // ...then e_hat(q) = eps_qi + a * 10^(-q/10), weighted by how many bases carried that Q. Weighting
    // matters: on a 2-colour instrument almost every base sits at one of four Q values, and an
    // unweighted fit would let a Q value seen a hundred times outvote one seen a billion.
    double sw = 0.0, sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0;
    int used = 0;
    bases = 0;
    for (size_t q = 0; q < counts.size(); ++q) {
        const uint64_t n = counts[q][0] + counts[q][1];
        bases += n;
        if (n < min_bases) continue;
        const double w = static_cast<double>(n);
        const double x = std::pow(10.0, -static_cast<double>(q) / 10.0);
        const double y = static_cast<double>(counts[q][1]) / w;
        sw += w; sx += w * x; sy += w * y; sxx += w * x * x; sxy += w * x * y;
        ++used;
    }
    const double det = sw * sxx - sx * sx;
    if (used < 2 || det <= 0.0) {
        fitted = false;
        return;
    }
    slope = (sw * sxy - sx * sy) / det;
    quality_independent = (sy - slope * sx) / sw;
    // A negative intercept is not an error rate. It means the two best Q values disagree with the
    // straight line by more than the line's own span, so report no floor rather than a fiction.
    if (quality_independent < 0.0) quality_independent = 0.0;
    if (slope < 0.0) slope = 0.0;
    fitted = true;
}

double QualityCalibration::error(int q) const {
    const double nominal = phred_error(static_cast<uint8_t>(std::clamp(q, 0, 60)));
    if (!fitted) return nominal;
    // Slope only. The intercept is the synthesised anchor's own defect rate; folding it in would
    // add ~4e-3 to every base likelihood in the pipeline on the strength of the primer's quality.
    return std::clamp(slope * nominal, 1e-7, 0.75);
}

void CheckoutCounters::merge(const CheckoutCounters& o) {
    total += o.total;
    assigned += o.assigned;
    unmatched += o.unmatched;
    ambiguous += o.ambiguous;
    short_payload += o.short_payload;
    bad_umi += o.bad_umi;
    normalised += o.normalised;
    if (per_sample.size() < o.per_sample.size()) per_sample.resize(o.per_sample.size(), 0);
    for (size_t i = 0; i < o.per_sample.size(); ++i) per_sample[i] += o.per_sample[i];
    calibration.merge(o.calibration);
    if (umi_phred.size() < o.umi_phred.size()) umi_phred.resize(o.umi_phred.size());
    for (size_t i = 0; i < o.umi_phred.size(); ++i) {
        for (size_t q = 0; q < 61; ++q) umi_phred[i][q] += o.umi_phred[i][q];
    }
    if (payload_len.size() < o.payload_len.size()) payload_len.resize(o.payload_len.size());
    for (size_t i = 0; i < o.payload_len.size(); ++i) {
        for (size_t L = 0; L < kPayloadHistLen; ++L) payload_len[i][L] += o.payload_len[i][L];
    }
    trimmed_bases += o.trimmed_bases;
    for (const auto& kv : o.index_pairs) index_pairs[kv.first] += kv.second;
    for (const auto& kv : o.tiles) tiles[kv.first] += kv.second;
}

Checkout::Checkout(const PatternSet& patterns, CheckoutParams params)
    : patterns_(patterns), params_(std::move(params)) {
    counters_.per_sample.assign(patterns.size(), 0);
    counters_.umi_phred.assign(patterns.size(), {});
    counters_.payload_len.assign(patterns.size(), {});
}

std::string Checkout::header_tags(const std::string& umi, const std::string& umi_qual,
                                  const std::string& sample, const std::string& cell,
                                  const std::string& cell_qual) {
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
    // CB/CY are the SAM standard tags for a cell barcode and its qualities, which is what Cell
    // Ranger, STARsolo and alevin all write and every downstream tool reads.
    if (!cell.empty()) {
        if (!out.empty()) out += "\t";
        out += "CB:Z:";
        out += cell;
        if (!cell_qual.empty()) {
            out += "\tCY:Z:";
            out += cell_qual;
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
    CheckoutPair p = process_pair(seq, qual, {}, {}, scratch_);
    CheckoutRead out;
    out.ok = p.ok;
    out.sample = p.sample;
    out.umi = std::move(p.umi);
    out.umi_qual = std::move(p.umi_qual);
    out.cell = std::move(p.cell);
    out.cell_qual = std::move(p.cell_qual);
    out.seq = p.seq1;
    out.qual = p.qual1;
    out.score = p.score;
    return out;
}

CheckoutPair Checkout::process_pair(std::string_view seq1, std::string_view qual1,
                                    std::string_view seq2, std::string_view qual2,
                                    CheckoutScratch& scratch) {
    ++counters_.total;
    CheckoutPair out;
    const bool paired = !seq2.empty();

    PatternSet::Assignment a = patterns_.assign(seq1, qual1, params_.match);
    bool normalised = false;

    // Only if R1 came up empty. An ambiguous placement is a different diagnosis -- the tag is
    // there and the sheet cannot resolve it -- and looking elsewhere would not help.
    if (!a.ambiguous && a.sample < 0) {
        if (paired) {
            PatternSet::Assignment b = patterns_.assign(seq2, qual2, params_.match);
            if (!b.ambiguous && b.sample >= 0) {
                a = b;
                normalised = true;
            }
        } else {
            scratch.seq1.assign(seq1);
            scratch.qual1.assign(qual1);
            if (scratch.qual1.size() == scratch.seq1.size()) {
                reverse_complement(scratch.seq1, scratch.qual1);
                PatternSet::Assignment b =
                    patterns_.assign(scratch.seq1, scratch.qual1, params_.match);
                if (!b.ambiguous && b.sample >= 0) {
                    a = b;
                    normalised = true;
                }
            }
        }
    }

    if (a.ambiguous) {
        ++counters_.ambiguous;
        return out;
    }
    if (a.sample < 0) {
        ++counters_.unmatched;
        return out;
    }

    // The mate that turned out to carry the tag becomes R1, so everything downstream sees one
    // orientation.
    std::string_view seq = seq1, qual = qual1, mate_seq = seq2, mate_qual = qual2;
    if (normalised) {
        if (paired) {
            seq = seq2;
            qual = qual2;
            mate_seq = seq1;
            mate_qual = qual1;
        } else {
            seq = scratch.seq1;
            qual = scratch.qual1;
        }
    }

    // The slave barcode, on the other mate. Both halves or nothing: a dual-end design that
    // accepts a read on the master alone emits half-length UMIs next to full-length ones, and
    // every collision estimate downstream is then computed over two barcode spaces at once.
    PatternMatch m = a.match;
    if (patterns_.has_slave(static_cast<size_t>(a.sample))) {
        if (!paired) {
            throw MigecError("checkout: sample '" + patterns_.samples()[static_cast<size_t>(a.sample)] +
                             "' declares a slave barcode, which needs a second mate to match "
                             "against -- run checkout with both R1 and R2");
        }
        const PatternMatch sm = patterns_.slave(static_cast<size_t>(a.sample))
                                    .match(mate_seq, mate_qual, params_.match);
        if (!sm.found) {
            ++counters_.unmatched;
            return out;
        }
        m.umi += sm.umi;
        m.umi_qual += sm.umi_qual;
    }
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

    size_t begin = 0;
    const size_t end = seq.size();
    // Everything 5' of the payload goes: adapter, sample tag, UMI. This is synthetic sequence and
    // leaving it in costs soft-clips at best and mismapping at worst.
    if (params_.trim == TrimMode::kPattern) begin = static_cast<size_t>(m.payload_begin);
    if (begin > end) begin = end;

    // The PAIR must carry payload, not this mate alone. On 10x the barcode read is 26 nt of cell
    // barcode and UMI and nothing else, so trimming leaves R1 empty while R2 holds the whole
    // cDNA -- checking R1 by itself drops 100% of a perfectly good library as "too short".
    const size_t payload = (end - begin) + mate_seq.size();
    if (static_cast<int>(payload) < params_.min_payload) {
        ++counters_.short_payload;
        return out;
    }

    out.ok = true;
    out.sample = a.sample;
    out.umi = m.umi;
    out.umi_qual = m.umi_qual;
    out.cell = m.cell;
    out.cell_qual = m.cell_qual;
    out.seq1 = seq.substr(begin, end - begin);
    out.qual1 = qual.empty() ? std::string_view() : qual.substr(begin, end - begin);
    // The mate is passed through whole, even when a slave pattern matched in it: trimming it would
    // need its own payload_begin carried alongside, and the mate's barcode bases are already in
    // the UMI by then. What the reader loses is a few synthetic bases at the mate's 5' end.
    out.seq2 = mate_seq;
    out.qual2 = mate_qual;
    out.normalised = normalised;
    out.score = m.score;
    ++counters_.assigned;
    if (normalised) ++counters_.normalised;
    ++counters_.per_sample[static_cast<size_t>(a.sample)];
    // The trim's own QC. A pattern placed one base off still matches and still trims -- it just
    // leaves every payload one base long or short, which no counter of matched reads can show.
    counters_.trimmed_bases += begin;
    ++counters_.payload_len[static_cast<size_t>(a.sample)]
                           [std::min(end - begin, kPayloadHistLen - 1)];
    for (char ch : m.umi_qual) ++counters_.umi_phred[static_cast<size_t>(a.sample)][phred_from_char(ch)];
    // The pattern's own constant bases are known sequence, so a disagreement there is an
    // instrument error and nothing else -- the only free calibration standard in the read.
    patterns_.pattern(static_cast<size_t>(a.sample))
        .calibrate(seq, qual, m.offset, counters_.calibration.by_position);
    return out;
}

// ---------------------------------------------------------------------------------------------
// The whole-file driver.
//
// Reads are independent, so the parallel part is trivial; what is not trivial is keeping the
// output byte-identical whatever the thread count, because a demultiplexer whose output depends
// on -t is a demultiplexer whose results cannot be compared between runs. The shape here is a
// round of chunks read serially, matched and compressed in parallel, and appended to the output
// files in chunk order. Matching AND compression happen on the workers; the serial stage does
// nothing but fwrite, which is what makes the scaling real rather than nominal.

namespace {

// Appends already-compressed bytes. Deliberately not FastqWriter: that one owns the compression,
// and here the compression has already happened on a worker thread.
struct BlockFile {
    std::FILE* f = nullptr;
    bool wrote = false;
    std::string path;   // for the error message: a path is what makes it actionable

    explicit BlockFile(const std::string& p) : f(std::fopen(p.c_str(), "wb")), path(p) {
        if (!f) throw MigecError("checkout: cannot open " + p);
    }
    // The destructor cannot throw -- it runs during unwinding -- so it drops the handle silently
    // and leaves the error to the explicit close() every success path makes.
    ~BlockFile() {
        if (f) std::fclose(f);
        f = nullptr;
    }
    BlockFile(const BlockFile&) = delete;
    BlockFile& operator=(const BlockFile&) = delete;

    void append(std::string_view b) {
        if (b.empty()) return;
        if (std::fwrite(b.data(), 1, b.size(), f) != b.size()) {
            throw MigecError("checkout: short write");
        }
        wrote = true;
    }
    void close() {
        if (!f) return;
        // Never: fclose is where a full disk shows up. Every append above went into the C
        // library's buffer, so ignoring the flush reports a successful run over a truncated
        // FASTQ -- and a truncated gzip member is a file the next stage refuses to read.
        const bool ok = std::fclose(f) == 0;
        f = nullptr;
        if (!ok) throw MigecError("checkout: could not flush " + path + " -- is the disk full?");
    }
};

// One chunk of input, owned. Fields are offsets into a single arena rather than four std::strings
// per read: at 100 M reads the difference is 800 M allocations.
struct Chunk {
    struct Rec {
        uint32_t off[4];  // name, comment, seq, qual
        uint32_t len[4];
    };
    std::string arena;
    std::vector<Rec> a, b;  // b is empty for single-end input

    void clear() {
        arena.clear();
        a.clear();
        b.clear();
    }
    std::string_view field(const Rec& r, int i) const {
        return std::string_view(arena.data() + r.off[i], r.len[i]);
    }
};

Chunk::Rec append_record(Chunk& c, const FastqRecord& r) {
    Chunk::Rec o{};
    const std::string_view src[4] = {r.name, r.comment, r.seq, r.qual};
    for (int i = 0; i < 4; ++i) {
        o.off[i] = static_cast<uint32_t>(c.arena.size());
        o.len[i] = static_cast<uint32_t>(src[i].size());
        c.arena.append(src[i]);
    }
    return o;
}

size_t read_chunk(FastqReader& r1, FastqReader* r2, size_t n, Chunk& c, uint64_t limit,
                  uint64_t& seen) {
    c.clear();
    FastqRecord a, b;
    // The read limit is applied HERE, at the intake, rather than by stopping the workers: a chunk
    // is the unit of both parallelism and output order, so cutting inside one would make the last
    // chunk's length depend on the thread count.
    while (c.a.size() < n && (!limit || seen < limit) && r1.next(a)) {
        ++seen;
        if (r2) {
            if (!r2->next(b)) throw MigecError("checkout: R2 ended before R1");
            Chunk::Rec ra = append_record(c, a);
            Chunk::Rec rb = append_record(c, b);
            c.a.push_back(ra);
            c.b.push_back(rb);
        } else {
            c.a.push_back(append_record(c, a));
        }
    }
    return c.a.size();
}

// One assigned read, staged for `.mig` output. The payload is copied into the worker's own arena
// rather than pointed at: `CheckoutPair`'s sequence views live in the scratch buffer that the very
// next read overwrites, and trimming means they are not the input bytes either.
struct MigStaged {
    uint64_t cell = 0, umi = 0, src_index = 0;
    uint32_t writer = 0;  // sample * n_buckets + bucket; also who owns writing it
    // seq1, qual1, seq2, qual2, then the BARCODE's own quality: the UMI's and the cell's. The
    // last two are what `refine`'s posterior weighs at the position that differs, and a `.mig`
    // that dropped them would hand refine a weaker model than the FASTQ route does.
    uint32_t off[6] = {0, 0, 0, 0, 0, 0};
    uint32_t len[6] = {0, 0, 0, 0, 0, 0};
    uint16_t flags = 0;
    uint8_t umi_minq = 0, cell_minq = 0;
};

struct Worker {
    std::unique_ptr<Checkout> co;
    CheckoutScratch scratch;
    std::vector<std::string> out1, out2;  // formatted FASTQ, one buffer per sample
    std::vector<std::string> z1, z2;      // ...and the same, compressed on this thread
    std::string un1, un2, zun1, zun2;
    // `.mig` mode only, and empty otherwise: the staged records of the chunk this worker matched,
    // in input order, and the arena their sequence and quality live in.
    std::vector<MigStaged> staged;
    std::string mig_arena;
    // (sample, packed UMI) per assigned read. Folded into the shared counters by the serial
    // stage, so there is one UMI counter per sample rather than one per thread per sample --
    // which at eight threads would be eight times the largest allocation in the process.
    std::vector<std::pair<uint32_t, uint64_t>> umis;
};

// An empty buffer still has a gzip member's worth of header and trailer, and there is one member
// per chunk per sample: on a 96-plex sheet a sample that is absent from most chunks accumulates
// megabytes of nothing. Emitting no bytes at all is equally valid -- a gzip stream is the
// concatenation of its members, and zero members decompresses to an empty file.
void gzip_member_or_nothing(std::string_view in, std::string& out, int level) {
    if (in.empty()) {
        out.clear();
        return;
    }
    gzip_member(in, out, level);
}

// How `.mig` output partitions this run. `on` false is the ordinary FASTQ path and everything else
// here is ignored.
struct MigLayout {
    bool on = false;
    int bits = 0;
    size_t n_buckets = 1;
};

// The (i7, i5) index pair out of an Illumina read header, or two empty strings.
//
// The comment is `<read>:<is filtered>:<control>:<index>`, and the index is `i7+i5` on a
// dual-indexed run, `i7` alone on a single-indexed one, and absent on anything that has been
// through a tool that rewrote the header. Never guess past that: a header migec does not recognise
// yields no pair and the table simply has nothing to say, which is the truth.
std::pair<std::string_view, std::string_view> index_pair(std::string_view comment) {
    const size_t last = comment.rfind(':');
    if (last == std::string_view::npos || last + 1 >= comment.size()) return {};
    std::string_view field = comment.substr(last + 1);
    // A tag block, not an index: the comment was already rewritten by something.
    if (field.find('\t') != std::string_view::npos) return {};
    const size_t plus = field.find('+');
    if (plus == std::string_view::npos) return {field, {}};
    return {field.substr(0, plus), field.substr(plus + 1)};
}

// Min Phred over a barcode's quality string, capped at 60 -- the `.mig` record's own field, and
// the evidence `refine` uses when the count ratio has nothing to say.
uint8_t min_phred(std::string_view qual) {
    int lo = 60;
    for (char ch : qual) lo = std::min(lo, static_cast<int>(ch) - 33);
    return static_cast<uint8_t>(std::clamp(lo, 0, 60));
}

void process_chunk(const Chunk& c, bool paired, bool write_unmatched, int gzip_level,
                   const std::vector<std::string>& ids, const std::vector<uint32_t>& file_of,
                   Worker& w, const MigLayout& mig = {}, uint64_t base_index = 0) {
    for (auto& s : w.out1) s.clear();
    for (auto& s : w.out2) s.clear();
    w.un1.clear();
    w.un2.clear();
    w.umis.clear();
    w.staged.clear();
    // Assigned into, never freed: the arena is reused chunk after chunk, so the payload copy costs
    // a memcpy and not an allocation.
    w.mig_arena.clear();

    std::string tags;
    for (size_t i = 0; i < c.a.size(); ++i) {
        const std::string_view n1 = c.field(c.a[i], 0), m1 = c.field(c.a[i], 1);
        const std::string_view s1 = c.field(c.a[i], 2), q1 = c.field(c.a[i], 3);
        std::string_view n2, m2, s2, q2;
        if (paired) {
            n2 = c.field(c.b[i], 0);
            m2 = c.field(c.b[i], 1);
            s2 = c.field(c.b[i], 2);
            q2 = c.field(c.b[i], 3);
        }

        // Before anything replaces the comment: the instrument's own index pair. Counted for
        // EVERY read, matched or not -- a hopped read that matches no pattern is still evidence of
        // hopping, and restricting the table to assigned reads would hide exactly the population
        // it exists to measure.
        {
            uint32_t lane = 0, tile = 0;
            if (parse_lane_tile(n1, &lane, &tile)) {
                ++w.co->counters_mutable().tiles[{lane, tile}];
            }
            const auto ix = index_pair(m1);
            if (!ix.first.empty()) {
                ++w.co->counters_mutable().index_pairs[{std::string(ix.first),
                                                        std::string(ix.second)}];
            }
        }

        CheckoutPair r = w.co->process_pair(s1, q1, s2, q2, w.scratch);
        if (!r.ok) {
            if (write_unmatched) {
                append_fastq(w.un1, n1, m1, s1, q1);
                if (paired) append_fastq(w.un2, n2, m2, s2, q2);
            }
            continue;
        }

        // Rows sharing a sample id share one output file and one UMI counter.
        const size_t s = file_of[static_cast<size_t>(r.sample)];
        bool umi_has_n = false, cell_has_n = false;
        const uint64_t umi_key = pack_barcode(r.umi, &umi_has_n);
        const uint64_t cell_key = r.cell.empty() ? 0 : pack_barcode(r.cell, &cell_has_n);
        // The molecule's whole key, cell then UMI, which is what `refine` and `assemble` group on.
        // `pack_barcode` puts base 0 in the top bits, so shifting the UMI down by the cell's width
        // lands it exactly where `pack_barcode(cell + umi)` would -- without the concatenation,
        // which would be an allocation per read.
        const uint64_t mol_key =
            r.cell.empty() ? umi_key : (cell_key | (umi_key >> (2 * r.cell.size())));
        if (mig.on) {
            MigStaged st;
            st.cell = cell_key;
            st.umi = umi_key;
            st.src_index = base_index + i;
            // Partition on the cell when there is one, exactly as `assemble` does: every read of a
            // cell then lands in one bucket, which is what makes a per-cell scope local.
            st.writer = static_cast<uint32_t>(
                s * mig.n_buckets +
                bucket_of(r.cell.empty() ? umi_key : cell_key, mig.bits));
            // What has ALREADY been applied, never what remains: a swapped pair is stored swapped
            // and a reverse-complemented single read is stored reverse-complemented, so `assemble`
            // must not re-orient anything.
            st.flags = static_cast<uint16_t>((paired ? 0 : kSingleEnd) |
                                             (umi_has_n ? kUmiHasN : 0) |
                                             (cell_has_n ? kCellHasN : 0) |
                                             (r.normalised ? (paired ? kMatesSwapped : kRevComp1)
                                                           : 0));
            st.umi_minq = min_phred(r.umi_qual);
            st.cell_minq = r.cell_qual.empty() ? 60 : min_phred(r.cell_qual);
            const std::string_view payload[6] = {r.seq1,     r.qual1,  r.seq2,
                                                 r.qual2,    r.umi_qual, r.cell_qual};
            for (int f = 0; f < 6; ++f) {
                st.off[f] = static_cast<uint32_t>(w.mig_arena.size());
                st.len[f] = static_cast<uint32_t>(payload[f].size());
                w.mig_arena.append(payload[f]);
            }
            w.staged.push_back(st);
        } else {
            tags = Checkout::header_tags(r.umi, r.umi_qual, ids[static_cast<size_t>(r.sample)],
                                         r.cell, r.cell_qual);
            // When the mates were swapped the names travel with them.
            std::string_view name1 = n1, name2 = n2;
            if (r.normalised && paired) std::swap(name1, name2);
            append_fastq(w.out1[s], name1, tags, r.seq1, r.qual1);
            if (paired) append_fastq(w.out2[s], name2, tags, r.seq2, r.qual2);
        }
        w.umis.emplace_back(static_cast<uint32_t>(s), mol_key);
    }
    // The per-sample buffers are empty in `.mig` mode -- nothing was formatted -- but unmatched
    // reads are FASTQ either way: they have no barcode, so there is no bucket to put them in.
    for (size_t s = 0; !mig.on && s < w.out1.size(); ++s) {
        gzip_member_or_nothing(w.out1[s], w.z1[s], gzip_level);
        if (paired) gzip_member_or_nothing(w.out2[s], w.z2[s], gzip_level);
    }
    if (write_unmatched) {
        gzip_member_or_nothing(w.un1, w.zun1, gzip_level);
        if (paired) gzip_member_or_nothing(w.un2, w.zun2, gzip_level);
    }
}

}  // namespace

CheckoutStats run_checkout(const PatternSet& patterns, const CheckoutParams& params,
                           const CheckoutRequest& request) {
    Stopwatch clock;
    const bool paired = !request.r2.empty();
    const size_t n_samples = patterns.size();
    if (n_samples == 0) throw MigecError("checkout: no barcode patterns");

    int nthreads = request.threads > 0 ? request.threads
                                       : static_cast<int>(hardware_threads());
    if (nthreads < 1) nthreads = 1;
    const size_t chunk = request.chunk_reads ? request.chunk_reads : 16384;

    const std::vector<std::string>& ids = patterns.samples();

    // Several rows may declare the same sample id -- that is how a sample sequenced with more than
    // one tag is written in a MIGEC barcode table. They are one sample: one output file, one UMI
    // counter. One file per row would open the same path twice and interleave two FILE* into it,
    // which does not even produce a valid gzip stream, and the summary would report success.
    CheckoutStats stats;
    std::vector<uint32_t> file_of(n_samples);
    // The first row that declared each sample. A sample with two tags has two rows and one file,
    // and the barcode lengths a `.mig` header needs are a property of the sample.
    std::vector<size_t> row_of_file;
    for (size_t i = 0; i < n_samples; ++i) {
        auto it = std::find(stats.sample_ids.begin(), stats.sample_ids.end(), ids[i]);
        if (it == stats.sample_ids.end()) {
            validate_sample_id(ids[i], "checkout");
            file_of[i] = static_cast<uint32_t>(stats.sample_ids.size());
            stats.sample_ids.push_back(ids[i]);
            row_of_file.push_back(i);
            // Never: the counter is keyed on the WHOLE barcode -- cell then UMI -- not on the UMI
            // alone. A molecule is sample + cell + UMI; UMIs repeat across cells by design, so a
            // UMI-keyed counter merges every cell's copy of one UMI into a single entry. On a 10x
            // library that read 221,026 distinct barcodes where there are 311,962 molecules, which
            // put the depth 1.41x high and the space 21% occupied against a true 3e-6 -- so the
            // saturation warning and `err_unreliable` fired on an artefact of the pooling.
            stats.umi_counts.emplace_back(patterns.cell_length(i) + patterns.umi_length(i));
            stats.sample_cell_length.push_back(patterns.cell_length(i));
        } else {
            file_of[i] = static_cast<uint32_t>(it - stats.sample_ids.begin());
            // The rows would otherwise write barcodes of two lengths into one counter, and the
            // composition, the collision statistics and the correction would all be nonsense.
            if (patterns.cell_length(i) + patterns.umi_length(i) !=
                stats.umi_counts[file_of[i]].length()) {
                throw MigecError("checkout: sample '" + ids[i] +
                                 "' is declared with two different barcode lengths");
            }
        }
    }
    const size_t n_files = stats.sample_ids.size();

    // Bound the counters. Past the budget each one range-partitions to disk and everything that
    // reads it -- the histogram, the composition, the correction -- streams a bucket at a time.
    //
    // Note: the budget is for the RUN and is divided by the samples, because a 96-plex sheet holds
    // 96 counters. The 16 MB floor is where partitioning stops paying for itself: a counter that
    // small is cheaper to hold than to write out and read back.
    //
    // Note: the partition is on the top bits of the key, so it is capped at half the barcode --
    // correction runs a second pass on keys rotated by the width of the prefix, and the two
    // prefixes have to be disjoint. A barcode too short to hold two prefixes is also too short to
    // produce a counter worth partitioning.
    if (request.umi_budget_bytes) {
        const std::string dir = request.umi_spill_dir.empty()
                                    ? request.out_prefix + ".umi_spill"
                                    : request.umi_spill_dir;
        stats.umi_spill = std::make_shared<UmiSpillDir>(UmiSpillDir{dir});
        // Never: the floor is capped by the budget itself. A fixed 16 MB floor silently overrode a
        // caller who asked for less, which makes the budget untestable at any corpus small enough
        // to run -- and a budget that only applies above 16 MB is not the property being claimed.
        const size_t floor = std::min<size_t>(size_t{16} << 20, request.umi_budget_bytes);
        const size_t per_counter = std::max(floor, request.umi_budget_bytes / n_files);
        for (size_t f = 0; f < n_files; ++f) {
            const int L = stats.umi_counts[f].length();
            const int bits = std::min(8, (L / 2) * 2);
            if (bits < 1) continue;
            stats.umi_counts[f].enable_spill(dir + "/" + stats.sample_ids[f], per_counter, bits);
        }
    }

    // `.mig` output: one range partition of the reads, per sample, written straight from the
    // workers -- which is the partition `assemble` would otherwise build for itself.
    //
    // Note: the open-file budget is for the RUN, not per sample. 256 writers is already more than
    // polite, and a 96-plex sheet holding 256 buckets each would be 24,576 of them; each sample of
    // such a sheet also holds a 96th of the reads, so a couple of buckets is the proportionate
    // answer rather than a compromise.
    MigLayout mig;
    mig.on = request.mig_output;
    if (mig.on) {
        int bits = request.mig_bucket_bits;
        if (bits <= 0) {
            bits = kMaxMigBucketBits;
            while (bits > 0 && (n_files << bits) > kMaxMigWriters) --bits;
        }
        mig.bits = bits;
        mig.n_buckets = size_t{1} << bits;
        stats.mig_output = true;
        stats.mig_bucket_bits = bits;
    }

    const std::string suffix1 = paired ? "_R1.fq.gz" : ".fq.gz";
    const std::string suffix2 = "_R2.fq.gz";

    std::vector<std::unique_ptr<BlockFile>> w1, w2;
    for (const std::string& id : stats.sample_ids) {
        if (mig.on) break;
        w1.push_back(std::make_unique<BlockFile>(request.out_prefix + id + suffix1));
        if (paired) w2.push_back(std::make_unique<BlockFile>(request.out_prefix + id + suffix2));
    }
    // One writer per (sample, bucket), opened lazily: a bucket that never receives a read is a
    // file that is never created, which on a fine partition of a small sample is most of them.
    std::vector<std::unique_ptr<MigWriter>> migw(mig.on ? n_files * mig.n_buckets : 0);
    std::vector<std::string> mig_path(migw.size());
    const size_t mig_block_bytes =
        migw.empty() ? 0
                     : std::clamp<size_t>(kMigWriterBudgetBytes / migw.size(), 256u << 10, 4u << 20);
    std::unique_ptr<BlockFile> u1, u2;
    if (request.write_unmatched) {
        u1 = std::make_unique<BlockFile>(request.out_prefix + "unmatched" + suffix1);
        if (paired) u2 = std::make_unique<BlockFile>(request.out_prefix + "unmatched" + suffix2);
    }

    stats.threads = nthreads;

    std::vector<Worker> workers(static_cast<size_t>(nthreads));
    for (Worker& w : workers) {
        w.co = std::make_unique<Checkout>(patterns, params);
        w.out1.resize(n_files);
        w.z1.resize(n_files);
        if (paired) {
            w.out2.resize(n_files);
            w.z2.resize(n_files);
        }
    }

    std::vector<Chunk> chunks(static_cast<size_t>(nthreads));
    for (Chunk& c : chunks) c.arena.reserve(chunk * 256);

    uint64_t seen = 0;
    FastqReader r1(request.r1);
    std::unique_ptr<FastqReader> r2;
    if (paired) r2 = std::make_unique<FastqReader>(request.r2);

    for (;;) {
        size_t filled = 0;
        // Where each chunk starts in the input. `src_index` is a read's position in the file and
        // is the sort tiebreak `.mig` depends on, so it cannot be handed out by a worker.
        std::vector<uint64_t> chunk_base(chunks.size(), 0);
        for (; filled < chunks.size(); ++filled) {
            chunk_base[filled] = seen;
            if (read_chunk(r1, r2.get(), chunk, chunks[filled], request.limit_reads, seen) == 0) {
                break;
            }
        }
        if (filled == 0) break;

        // An exception thrown on a worker would otherwise propagate out of the thread function and
        // call std::terminate -- an abort with no message, no stack and no output flushed, for
        // something as ordinary as a malformed pattern or a failed allocation.
        std::exception_ptr err;
        std::mutex err_mutex;
        auto guarded = [&](size_t t) {
            try {
                process_chunk(chunks[t], paired, request.write_unmatched, request.gzip_level, ids,
                              file_of, workers[t], mig, chunk_base[t]);
            } catch (...) {
                std::lock_guard<std::mutex> lock(err_mutex);
                if (!err) err = std::current_exception();
            }
        };

        // Chunk 0 stays on this thread; spawning a thread for it would only add a context switch.
        std::vector<std::thread> pool;
        pool.reserve(filled - 1);
        // Never: a failed spawn must not abort. std::system_error out of emplace_back would unwind
        // over threads that are still joinable, and ~thread() on a joinable thread is
        // std::terminate. Whatever could not be spawned is run on this thread instead: the chunks
        // are independent and the output order is fixed by the append loop below, so the bytes do
        // not move -- only the wall clock does.
        size_t spawned = 1;  // chunk 0 always runs here
        for (size_t t = 1; t < filled; ++t) {
            try {
                pool.emplace_back([&, t] { guarded(t); });
                ++spawned;
            } catch (...) {
                break;
            }
        }
        guarded(0);
        for (size_t t = spawned; t < filled; ++t) guarded(t);
        for (std::thread& th : pool) th.join();
        if (err) std::rethrow_exception(err);

        // Serial, in chunk order: this is what makes the output independent of the thread count.
        for (size_t t = 0; t < filled; ++t) {
            Worker& w = workers[t];
            for (size_t s = 0; !mig.on && s < n_files; ++s) {
                w1[s]->append(w.z1[s]);
                if (paired) w2[s]->append(w.z2[s]);
            }
            if (u1) u1->append(w.zun1);
            if (u2) u2->append(w.zun2);
            for (const auto& kv : w.umis) stats.umi_counts[kv.first].add(kv.second);
        }

        if (mig.on) {
            // Open the writers a record has just asked for. On the driver, because a `.mig` header
            // carries the barcode lengths and the sample id, and because two workers must never
            // race to create the same file.
            for (size_t t = 0; t < filled; ++t) {
                for (const MigStaged& st : workers[t].staged) {
                    if (migw[st.writer]) continue;
                    const size_t s = st.writer / mig.n_buckets;
                    const size_t b = st.writer % mig.n_buckets;
                    MigHeader header;
                    header.umi_len = static_cast<uint8_t>(stats.umi_counts[s].length());
                    header.cell_len = static_cast<uint8_t>(patterns.cell_length(row_of_file[s]));
                    header.bucket_index = static_cast<uint8_t>(b);
                    header.bucket_bits = static_cast<uint8_t>(mig.bits);
                    header.paired = paired;
                    header.barcode_quality = true;
                    header.sample_id = stats.sample_ids[s];
                    // Note: no quality calibration in the header. It is fitted from the whole run,
                    // and this file is opened while the run is still going -- `checkout.json`
                    // carries the fit, and a wrong table here would be worse than an absent one.
                    mig_path[st.writer] = request.out_prefix + stats.sample_ids[s] + "." +
                                          bucket_suffix(b) + ".mig";
                    migw[st.writer] =
                        std::make_unique<MigWriter>(mig_path[st.writer], header, mig_block_bytes);
                }
            }
            // Ownership, not locking: writer w is written by exactly one thread for the whole run,
            // so no writer state is ever shared. Every worker walks the chunks in input order and
            // each chunk forwards, so a record's position in its bucket is decided by the input and
            // never by who got there first -- which is what keeps `-t` out of the bytes.
            std::exception_ptr werr;
            std::mutex wmutex;
            parallel_for(static_cast<size_t>(nthreads), nthreads, [&](size_t owner, int) {
                try {
                    for (size_t t = 0; t < filled; ++t) {
                        const Worker& w = workers[t];
                        for (const MigStaged& st : w.staged) {
                            if (st.writer % static_cast<uint32_t>(nthreads) != owner) continue;
                            MigRecord rec;
                            rec.cell = st.cell;
                            rec.umi = st.umi;
                            rec.src_index = st.src_index;
                            rec.flags = st.flags;
                            rec.umi_minq = st.umi_minq;
                            rec.cell_minq = st.cell_minq;
                            rec.seq1 = std::string_view(w.mig_arena).substr(st.off[0], st.len[0]);
                            rec.qual1 = std::string_view(w.mig_arena).substr(st.off[1], st.len[1]);
                            rec.seq2 = std::string_view(w.mig_arena).substr(st.off[2], st.len[2]);
                            rec.qual2 = std::string_view(w.mig_arena).substr(st.off[3], st.len[3]);
                            rec.qual_umi = std::string_view(w.mig_arena).substr(st.off[4], st.len[4]);
                            rec.qual_cell =
                                std::string_view(w.mig_arena).substr(st.off[5], st.len[5]);
                            migw[st.writer]->write(rec);
                        }
                    }
                } catch (...) {
                    std::lock_guard<std::mutex> lock(wmutex);
                    if (!werr) werr = std::current_exception();
                }
            });
            if (werr) std::rethrow_exception(werr);
        }
        if (filled < chunks.size()) break;  // hit EOF part-way through the round
    }

    if (r2 && !(request.limit_reads && seen >= request.limit_reads)) {
        FastqRecord leftover;
        if (r2->next(leftover)) throw MigecError("checkout: R1 ended before R2");
    }

    // A file that received no chunk is zero bytes, and a zero-byte file is not a gzip stream --
    // `gzip -t` and `zcat` both reject it. A sample that got no reads gets exactly one empty
    // member, which reads as an empty FASTQ everywhere.
    std::string empty_member;
    gzip_member(std::string_view(), empty_member, request.gzip_level);
    auto finish = [&empty_member](BlockFile& f) {
        if (!f.wrote) f.append(empty_member);
        f.close();
    };
    for (auto& w : w1) finish(*w);
    for (auto& w : w2) finish(*w);
    if (u1) finish(*u1);
    if (u2) finish(*u2);

    // Note: a `.mig` bucket that received no read is not created at all, rather than created empty.
    // A bucket is addressed by the header inside it and by nothing else -- `assemble` reads the
    // files it is given, in bucket order -- so an absent bucket is an absent range of the key
    // space, which is exactly what it means. This is the opposite of the FASTQ case above, where
    // an empty file has to exist because the next tool globs for it.
    for (size_t i = 0; i < migw.size(); ++i) {
        if (!migw[i]) continue;
        migw[i]->close();
        stats.mig_paths.push_back(mig_path[i]);
    }

    for (const Worker& w : workers) stats.counters.merge(w.co->counters());
    stats.counters.per_sample.resize(n_samples, 0);
    stats.sample_reads.assign(n_files, 0);
    stats.sample_phred.assign(n_files, {});
    stats.counters.umi_phred.resize(n_samples);
    stats.counters.payload_len.resize(n_samples);
    stats.sample_payload_len.assign(n_files, {});
    for (size_t i = 0; i < n_samples; ++i) {
        stats.sample_reads[file_of[i]] += stats.counters.per_sample[i];
        for (size_t q = 0; q < 61; ++q) stats.sample_phred[file_of[i]][q] += stats.counters.umi_phred[i][q];
        for (size_t L = 0; L < kPayloadHistLen; ++L) {
            stats.sample_payload_len[file_of[i]][L] += stats.counters.payload_len[i][L];
        }
    }
    // Fit once, after every worker's counts are in: the intercept is a property of the run, and
    // fitting per chunk would give as many answers as there were chunks.
    stats.counters.calibration.fit();
    stats.wall_seconds = clock.seconds();
    stats.reads_per_second =
        stats.wall_seconds > 0.0 ? static_cast<double>(stats.counters.total) / stats.wall_seconds
                                 : 0.0;
    stats.peak_rss_bytes = peak_rss_bytes();
    for (const UmiCounts& u : stats.umi_counts) {
        stats.umi_memory_bytes += u.memory_bytes();
        stats.umi_spilled = stats.umi_spilled || u.spilled();
    }
    return stats;
}

bool parse_lane_tile(std::string_view name, uint32_t* lane, uint32_t* tile) {
    // Split on ':' and take the fields by position. Seven fields is Casava 1.8+ (lane 4, tile 5);
    // five or more with a '#' in the last is the older form (lane 2, tile 3).
    std::vector<std::string_view> f;
    size_t pos = 0;
    while (pos <= name.size() && f.size() < 12) {
        const size_t end = std::min(name.find(':', pos), name.size());
        f.push_back(name.substr(pos, end - pos));
        pos = end + 1;
    }
    auto number = [](std::string_view v, uint32_t* out) {
        if (v.empty() || v.size() > 9) return false;
        uint32_t n = 0;
        for (char c : v) {
            if (c < '0' || c > '9') return false;
            n = n * 10 + static_cast<uint32_t>(c - '0');
        }
        *out = n;
        return true;
    };
    if (f.size() >= 7) return number(f[3], lane) && number(f[4], tile);
    if (f.size() == 5) return number(f[1], lane) && number(f[2], tile);
    return false;
}

IndexHopping estimate_index_hopping(
    const std::map<std::pair<std::string, std::string>, uint64_t>& pairs, double min_share) {
    IndexHopping h;
    h.min_share = min_share;
    if (pairs.empty()) return h;
    std::map<std::string, uint64_t> by_i7, by_i5;
    uint64_t total = 0;
    bool dual = false;
    for (const auto& kv : pairs) {
        by_i7[kv.first.first] += kv.second;
        by_i5[kv.first.second] += kv.second;
        total += kv.second;
        if (!kv.first.second.empty()) dual = true;
    }
    h.i7_indices = by_i7.size();
    h.i5_indices = by_i5.size();
    // One index, or one of each: there are no combinations, so nothing can be off-diagonal.
    if (!dual || h.i7_indices < 2 || h.i5_indices < 2) return h;
    h.estimable = true;
    for (const auto& kv : pairs) {
        const double share_i7 =
            static_cast<double>(kv.second) / static_cast<double>(by_i7[kv.first.first]);
        const double share_i5 =
            static_cast<double>(kv.second) / static_cast<double>(by_i5[kv.first.second]);
        if (share_i7 >= min_share && share_i5 >= min_share) {
            ++h.declared_pairs;
            h.reads_declared += kv.second;
        } else {
            ++h.hopped_pairs;
            h.reads_hopped += kv.second;
        }
    }
    h.rate = total ? static_cast<double>(h.reads_hopped) / static_cast<double>(total) : 0.0;
    return h;
}

UmiSpillDir::~UmiSpillDir() {
    // Best effort: a temp directory that outlived its readers is not worth an exception out of a
    // destructor, and the alternative -- leaving gigabytes of buckets behind after a failed run --
    // is what this exists to prevent.
    std::error_code ec;
    std::filesystem::remove_all(path, ec);
}

}  // namespace migec
