#include "migec/checkout.hpp"

#include <zlib.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <mutex>
#include <thread>
#include <utility>

#include "migec/fastq.hpp"
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
}

Checkout::Checkout(const PatternSet& patterns, CheckoutParams params)
    : patterns_(patterns), params_(std::move(params)) {
    counters_.per_sample.assign(patterns.size(), 0);
    counters_.umi_phred.assign(patterns.size(), {});
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

    size_t begin = 0;
    const size_t end = seq.size();
    // Everything 5' of the payload goes: adapter, sample tag, UMI. This is synthetic sequence and
    // leaving it in costs soft-clips at best and mismapping at worst.
    if (params_.trim == TrimMode::kPattern) begin = static_cast<size_t>(m.payload_begin);
    if (begin > end) begin = end;

    if (static_cast<int>(end - begin) < params_.min_payload) {
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
    // The mate is passed through whole. Trimming it would need its own tag, and dual-end barcodes
    // are not implemented yet -- see ROADMAP M2.
    out.seq2 = mate_seq;
    out.qual2 = mate_qual;
    out.normalised = normalised;
    out.score = m.score;
    ++counters_.assigned;
    if (normalised) ++counters_.normalised;
    ++counters_.per_sample[static_cast<size_t>(a.sample)];
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

// A complete gzip member for `in`. Concatenated members are themselves a valid gzip stream (RFC
// 1952 s2.2), which is what lets each worker compress its own chunk and the writer merely append
// the bytes. That matters more than it sounds: zlib runs at ~7 MB/s on random DNA at level 6 and
// ~137 MB/s at level 1, so compression left on the serial path caps checkout at a fraction of what
// the matcher can do, however many threads are matching.
void gzip_member(std::string_view in, std::string& out, int level) {
    z_stream zs{};
    // 15 + 16: a 32 kB window with a gzip wrapper rather than a zlib one.
    if (deflateInit2(&zs, level, Z_DEFLATED, 15 + 16, 8, Z_DEFAULT_STRATEGY) != Z_OK) {
        throw MigecError("checkout: deflateInit2 failed");
    }
    out.resize(deflateBound(&zs, static_cast<uLong>(in.size())) + 32);
    zs.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in.data()));
    zs.avail_in = static_cast<uInt>(in.size());
    zs.next_out = reinterpret_cast<Bytef*>(out.data());
    zs.avail_out = static_cast<uInt>(out.size());
    const int rc = deflate(&zs, Z_FINISH);
    const size_t produced = out.size() - zs.avail_out;
    deflateEnd(&zs);
    if (rc != Z_STREAM_END) throw MigecError("checkout: deflate failed");
    out.resize(produced);
}

// Appends already-compressed bytes. Deliberately not FastqWriter: that one owns the compression,
// and here the compression has already happened on a worker thread.
struct BlockFile {
    std::FILE* f = nullptr;
    bool wrote = false;

    explicit BlockFile(const std::string& path) : f(std::fopen(path.c_str(), "wb")) {
        if (!f) throw MigecError("checkout: cannot open " + path);
    }
    ~BlockFile() { close(); }
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
        if (f) {
            std::fclose(f);
            f = nullptr;
        }
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

size_t read_chunk(FastqReader& r1, FastqReader* r2, size_t n, Chunk& c) {
    c.clear();
    FastqRecord a, b;
    while (c.a.size() < n && r1.next(a)) {
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

void emit(std::string& dst, std::string_view name, std::string_view comment, std::string_view seq,
          std::string_view qual) {
    dst += '@';
    dst += name;
    if (!comment.empty()) {
        dst += ' ';
        dst += comment;
    }
    dst += '\n';
    dst += seq;
    dst += "\n+\n";
    dst += qual;
    dst += '\n';
}

struct Worker {
    std::unique_ptr<Checkout> co;
    CheckoutScratch scratch;
    std::vector<std::string> out1, out2;  // formatted FASTQ, one buffer per sample
    std::vector<std::string> z1, z2;      // ...and the same, compressed on this thread
    std::string un1, un2, zun1, zun2;
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

void process_chunk(const Chunk& c, bool paired, bool write_unmatched, int gzip_level,
                   const std::vector<std::string>& ids, const std::vector<uint32_t>& file_of,
                   Worker& w) {
    for (auto& s : w.out1) s.clear();
    for (auto& s : w.out2) s.clear();
    w.un1.clear();
    w.un2.clear();
    w.umis.clear();

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

        CheckoutPair r = w.co->process_pair(s1, q1, s2, q2, w.scratch);
        if (!r.ok) {
            if (write_unmatched) {
                emit(w.un1, n1, m1, s1, q1);
                if (paired) emit(w.un2, n2, m2, s2, q2);
            }
            continue;
        }

        // Rows sharing a sample id share one output file and one UMI counter.
        const size_t s = file_of[static_cast<size_t>(r.sample)];
        tags = Checkout::header_tags(r.umi, r.umi_qual, ids[static_cast<size_t>(r.sample)],
                                     r.cell, r.cell_qual);
        // When the mates were swapped the names travel with them.
        std::string_view name1 = n1, name2 = n2;
        if (r.normalised && paired) std::swap(name1, name2);
        emit(w.out1[s], name1, tags, r.seq1, r.qual1);
        if (paired) emit(w.out2[s], name2, tags, r.seq2, r.qual2);
        w.umis.emplace_back(static_cast<uint32_t>(s), pack_barcode(r.umi));
    }

    for (size_t s = 0; s < w.out1.size(); ++s) {
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
    for (size_t i = 0; i < n_samples; ++i) {
        auto it = std::find(stats.sample_ids.begin(), stats.sample_ids.end(), ids[i]);
        if (it == stats.sample_ids.end()) {
            file_of[i] = static_cast<uint32_t>(stats.sample_ids.size());
            stats.sample_ids.push_back(ids[i]);
            stats.umi_counts.emplace_back(patterns.pattern(i).umi_length());
        } else {
            file_of[i] = static_cast<uint32_t>(it - stats.sample_ids.begin());
            // The rows would otherwise write UMIs of two lengths into one counter, and the
            // composition, the collision statistics and the correction would all be nonsense.
            if (patterns.pattern(i).umi_length() != stats.umi_counts[file_of[i]].length()) {
                throw MigecError("checkout: sample '" + ids[i] +
                                 "' is declared with two different UMI lengths");
            }
        }
    }
    const size_t n_files = stats.sample_ids.size();

    const std::string suffix1 = paired ? "_R1.fq.gz" : ".fq.gz";
    const std::string suffix2 = "_R2.fq.gz";

    std::vector<std::unique_ptr<BlockFile>> w1, w2;
    for (const std::string& id : stats.sample_ids) {
        w1.push_back(std::make_unique<BlockFile>(request.out_prefix + id + suffix1));
        if (paired) w2.push_back(std::make_unique<BlockFile>(request.out_prefix + id + suffix2));
    }
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

    FastqReader r1(request.r1);
    std::unique_ptr<FastqReader> r2;
    if (paired) r2 = std::make_unique<FastqReader>(request.r2);

    for (;;) {
        size_t filled = 0;
        for (; filled < chunks.size(); ++filled) {
            if (read_chunk(r1, r2.get(), chunk, chunks[filled]) == 0) break;
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
                              file_of, workers[t]);
            } catch (...) {
                std::lock_guard<std::mutex> lock(err_mutex);
                if (!err) err = std::current_exception();
            }
        };

        // Chunk 0 stays on this thread; spawning a thread for it would only add a context switch.
        std::vector<std::thread> pool;
        pool.reserve(filled - 1);
        for (size_t t = 1; t < filled; ++t) pool.emplace_back([&, t] { guarded(t); });
        guarded(0);
        for (std::thread& th : pool) th.join();
        if (err) std::rethrow_exception(err);

        // Serial, in chunk order: this is what makes the output independent of the thread count.
        for (size_t t = 0; t < filled; ++t) {
            Worker& w = workers[t];
            for (size_t s = 0; s < n_files; ++s) {
                w1[s]->append(w.z1[s]);
                if (paired) w2[s]->append(w.z2[s]);
            }
            if (u1) u1->append(w.zun1);
            if (u2) u2->append(w.zun2);
            for (const auto& kv : w.umis) stats.umi_counts[kv.first].add(kv.second);
        }
        if (filled < chunks.size()) break;  // hit EOF part-way through the round
    }

    if (r2) {
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

    for (const Worker& w : workers) stats.counters.merge(w.co->counters());
    stats.counters.per_sample.resize(n_samples, 0);
    stats.sample_reads.assign(n_files, 0);
    stats.sample_phred.assign(n_files, {});
    stats.counters.umi_phred.resize(n_samples);
    for (size_t i = 0; i < n_samples; ++i) {
        stats.sample_reads[file_of[i]] += stats.counters.per_sample[i];
        for (size_t q = 0; q < 61; ++q) stats.sample_phred[file_of[i]][q] += stats.counters.umi_phred[i][q];
    }
    // Fit once, after every worker's counts are in: the intercept is a property of the run, and
    // fitting per chunk would give as many answers as there were chunks.
    stats.counters.calibration.fit();
    stats.wall_seconds = clock.seconds();
    stats.reads_per_second =
        stats.wall_seconds > 0.0 ? static_cast<double>(stats.counters.total) / stats.wall_seconds
                                 : 0.0;
    stats.peak_rss_bytes = peak_rss_bytes();
    for (const UmiCounts& u : stats.umi_counts) stats.umi_memory_bytes += u.memory_bytes();
    return stats;
}

}  // namespace migec
