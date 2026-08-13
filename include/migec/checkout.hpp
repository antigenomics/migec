// checkout: find the barcode pattern in each read, extract sample/UMI, trim, and hand the barcode
// to the read header so downstream tools can see it.

#ifndef MIGEC_CHECKOUT_HPP
#define MIGEC_CHECKOUT_HPP

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "migec/pattern.hpp"
#include "migec/umi_stats.hpp"

namespace migec {

enum class TrimMode {
    // Keep the read as it is. The UMI is in the header; the payload still carries the adapter.
    kNone,
    // Drop everything up to and including the matched pattern -- adapter, sample tag and UMI.
    // This is what you want before alignment: the tag is synthetic sequence and will not map.
    kPattern,
};

struct CheckoutParams {
    MatchParams match;
    TrimMode trim = TrimMode::kPattern;
    // Reads shorter than this after trimming are dropped as uninformative.
    int min_payload = 1;
    // Also require the UMI to have no ambiguous base. Off by default: an N in a UMI is a reason
    // to be less certain, not a reason to throw the molecule away.
    bool reject_umi_n = false;
    // Minimum Phred across the UMI bases. 0 disables. MIGEC used 15 and MAGERI 20, both as hard
    // drops; the default here is 0 because discarding a read for a low-quality UMI base loses
    // sequence information that the correction step can often recover.
    int min_umi_quality = 0;
};

struct CheckoutCounters {
    uint64_t total = 0;
    uint64_t assigned = 0;
    uint64_t unmatched = 0;
    uint64_t ambiguous = 0;
    uint64_t short_payload = 0;
    uint64_t bad_umi = 0;
    // Reads normalised onto the pattern's strand: the mate was swapped (paired) or the read was
    // reverse-complemented (single-end) because the pattern was only found the other way round.
    uint64_t normalised = 0;
    std::vector<uint64_t> per_sample;
    // Reported Phred over the *barcode* bases, per sample. 61 counters is 488 bytes a sample, and
    // it is what turns "the estimated error rate is 2.7e-4" into "...against 1.8e-3 predicted by
    // the quality the instrument reported", which is the comparison that says whether to believe
    // either number.
    std::vector<std::array<uint64_t, 61>> umi_phred;

    void merge(const CheckoutCounters& o);
};

struct CheckoutRead {
    bool ok = false;
    int sample = -1;
    std::string umi;
    std::string umi_qual;
    // Views into the input read, after trimming.
    std::string_view seq;
    std::string_view qual;
    double score = 0.0;
};

// Reusable buffers for the paired path. Strand normalisation has to reverse-complement, which
// cannot be done in a view, and a fresh allocation per read is a measurable fraction of the
// per-read cost -- so the caller keeps one of these and it stops allocating after the first read.
struct CheckoutScratch {
    std::string seq1, qual1, seq2, qual2;
};

struct CheckoutPair {
    bool ok = false;
    int sample = -1;
    std::string umi;
    std::string umi_qual;
    // Views into `scratch`, valid until the next call with the same scratch.
    std::string_view seq1, qual1, seq2, qual2;
    bool normalised = false;  // the mates were swapped, or the single read was rc'd
    double score = 0.0;
};

// Stateless apart from the pattern set and the counters.
class Checkout {
public:
    Checkout(const PatternSet& patterns, CheckoutParams params);

    CheckoutRead process(std::string_view seq, std::string_view qual);

    // Single-end when `seq2` is empty.
    //
    // The pattern is looked for in R1 first. On failure -- and only on failure, so the cost is
    // paid for reads that would otherwise be discarded -- it is looked for in R2, or in the
    // reverse complement for single-end input. When it turns up there the pair is swapped (or the
    // read flipped) so that everything downstream sees one orientation. A MIG holding both
    // orientations of the same molecule loses half its reads at consensus and nothing upstream
    // reports it, which is why this is not optional.
    CheckoutPair process_pair(std::string_view seq1, std::string_view qual1,
                              std::string_view seq2, std::string_view qual2,
                              CheckoutScratch& scratch);

    const CheckoutCounters& counters() const { return counters_; }

    // The SAM-style comment to append to a FASTQ header. TAB-separated, because `bwa mem -C` and
    // `minimap2 -y` copy the comment verbatim into the SAM record and it has to be conformant
    // there. `sample` may be empty.
    static std::string header_tags(const std::string& umi, const std::string& umi_qual,
                                   const std::string& sample);

private:
    const PatternSet& patterns_;
    CheckoutParams params_;
    CheckoutCounters counters_;
    CheckoutScratch scratch_;  // backs the single-read process() overload
};

// ---------------------------------------------------------------------------------------------
// The whole-file driver.

struct CheckoutRequest {
    std::string r1;                  // input FASTQ, plain or gzipped
    std::string r2;                  // empty for single-end
    std::string out_prefix;          // "<prefix><sample>_R1.fq.gz", or "<prefix><sample>.fq.gz"
    bool write_unmatched = false;
    // Level 1, not zlib's default 6. Read payload is close to incompressible, and on random DNA
    // level 6 runs at 7 MB/s against level 1's 137 MB/s for 13% more bytes. Paying 20x the CPU
    // for a tenth of the file is not a trade anyone would make deliberately.
    int gzip_level = 1;
    // 0 means one per hardware thread. Output is byte-identical whatever this is set to: reads are
    // processed in chunks and the chunks are written back in input order, so threading changes the
    // wall clock and nothing else.
    int threads = 0;
    // Bounds the per-worker buffers, and with them the only part of peak RSS that scales with
    // -t rather than with the library. 8192 x threads reads in flight costs ~5 MB per thread and
    // leaves the thread-spawn overhead per round under a tenth of a percent.
    size_t chunk_reads = 8192;
};

struct CheckoutStats {
    CheckoutCounters counters;  // per_sample is per pattern *row*, so a sample with two tags has two
    // Distinct sample ids, in first-appearance order. Several rows may carry the same id -- that is
    // how a sample sequenced with more than one tag is declared -- and they are one sample here:
    // one output file, one UMI counter. Opening a file per row instead means two handles on one
    // path, whose interleaved writes are not even a valid gzip stream.
    std::vector<std::string> sample_ids;
    std::vector<uint64_t> sample_reads;                    // parallel to sample_ids
    std::vector<std::array<uint64_t, 61>> sample_phred;    // parallel to sample_ids
    std::vector<UmiCounts> umi_counts;   // parallel to sample_ids
    double wall_seconds = 0.0;
    double reads_per_second = 0.0;
    size_t peak_rss_bytes = 0;
    size_t umi_memory_bytes = 0;  // the part of the above that is the UMI counters
    int threads = 1;
};

// Demultiplex, extract, trim and write. Throws MigecError on malformed input or an unwritable
// output path.
CheckoutStats run_checkout(const PatternSet& patterns, const CheckoutParams& params,
                           const CheckoutRequest& request);

}  // namespace migec

#endif  // MIGEC_CHECKOUT_HPP
