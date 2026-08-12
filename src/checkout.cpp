#include "migec/checkout.hpp"

#include <zlib.h>

#include <algorithm>
#include <cstdio>
#include <memory>
#include <thread>
#include <utility>

#include "migec/fastq.hpp"
#include "migec/resource.hpp"
#include "migec/types.hpp"

namespace migec {

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
}

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
    CheckoutPair p = process_pair(seq, qual, {}, {}, scratch_);
    CheckoutRead out;
    out.ok = p.ok;
    out.sample = p.sample;
    out.umi = std::move(p.umi);
    out.umi_qual = std::move(p.umi_qual);
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

void process_chunk(const Chunk& c, bool paired, bool write_unmatched, int gzip_level,
                   const std::vector<std::string>& ids, Worker& w) {
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

        const size_t s = static_cast<size_t>(r.sample);
        tags = Checkout::header_tags(r.umi, r.umi_qual, ids[s]);
        // When the mates were swapped the names travel with them.
        std::string_view name1 = n1, name2 = n2;
        if (r.normalised && paired) std::swap(name1, name2);
        emit(w.out1[s], name1, tags, r.seq1, r.qual1);
        if (paired) emit(w.out2[s], name2, tags, r.seq2, r.qual2);
        w.umis.emplace_back(static_cast<uint32_t>(s), pack_barcode(r.umi));
    }

    for (size_t s = 0; s < w.out1.size(); ++s) {
        gzip_member(w.out1[s], w.z1[s], gzip_level);
        if (paired) gzip_member(w.out2[s], w.z2[s], gzip_level);
    }
    if (write_unmatched) {
        gzip_member(w.un1, w.zun1, gzip_level);
        if (paired) gzip_member(w.un2, w.zun2, gzip_level);
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
    const std::string suffix1 = paired ? "_R1.fq.gz" : ".fq.gz";
    const std::string suffix2 = "_R2.fq.gz";

    std::vector<std::unique_ptr<BlockFile>> w1, w2;
    for (const std::string& id : ids) {
        w1.push_back(std::make_unique<BlockFile>(request.out_prefix + id + suffix1));
        if (paired) w2.push_back(std::make_unique<BlockFile>(request.out_prefix + id + suffix2));
    }
    std::unique_ptr<BlockFile> u1, u2;
    if (request.write_unmatched) {
        u1 = std::make_unique<BlockFile>(request.out_prefix + "unmatched" + suffix1);
        if (paired) u2 = std::make_unique<BlockFile>(request.out_prefix + "unmatched" + suffix2);
    }

    CheckoutStats stats;
    stats.threads = nthreads;
    for (size_t i = 0; i < n_samples; ++i) {
        stats.umi_counts.emplace_back(patterns.pattern(i).umi_length());
    }

    std::vector<Worker> workers(static_cast<size_t>(nthreads));
    for (Worker& w : workers) {
        w.co = std::make_unique<Checkout>(patterns, params);
        w.out1.resize(n_samples);
        w.z1.resize(n_samples);
        if (paired) {
            w.out2.resize(n_samples);
            w.z2.resize(n_samples);
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

        // Chunk 0 stays on this thread; spawning a thread for it would only add a context switch.
        std::vector<std::thread> pool;
        pool.reserve(filled - 1);
        for (size_t t = 1; t < filled; ++t) {
            pool.emplace_back([&, t] {
                process_chunk(chunks[t], paired, request.write_unmatched, request.gzip_level, ids,
                              workers[t]);
            });
        }
        process_chunk(chunks[0], paired, request.write_unmatched, request.gzip_level, ids,
                      workers[0]);
        for (std::thread& th : pool) th.join();

        // Serial, in chunk order: this is what makes the output independent of the thread count.
        for (size_t t = 0; t < filled; ++t) {
            Worker& w = workers[t];
            for (size_t s = 0; s < n_samples; ++s) {
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

    for (auto& w : w1) w->close();
    for (auto& w : w2) w->close();
    if (u1) u1->close();
    if (u2) u2->close();

    for (const Worker& w : workers) stats.counters.merge(w.co->counters());
    stats.counters.per_sample.resize(n_samples, 0);
    stats.wall_seconds = clock.seconds();
    stats.reads_per_second =
        stats.wall_seconds > 0.0 ? static_cast<double>(stats.counters.total) / stats.wall_seconds
                                 : 0.0;
    stats.peak_rss_bytes = peak_rss_bytes();
    for (const UmiCounts& u : stats.umi_counts) stats.umi_memory_bytes += u.memory_bytes();
    return stats;
}

}  // namespace migec
