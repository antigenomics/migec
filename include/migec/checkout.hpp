// checkout: find the barcode pattern in each read, extract sample/UMI, trim, and hand the barcode
// to the read header so downstream tools can see it.

#ifndef MIGEC_CHECKOUT_HPP
#define MIGEC_CHECKOUT_HPP

#include <array>
#include <cstdint>
#include <memory>
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

// What the reported Phred is actually worth, measured against the pattern's own constant bases.
//
// Never: The nominal `e = 10^(-q/10)` is not the error rate. RTA3 emits about four distinct Q values,
// so the mapping is coarse to begin with, and it is wrong by an order of magnitude by cycle and
// context -- which matters because every likelihood in this pipeline is computed from it.
//
// Fitting the measured table as
//
//     e_hat(q) = eps_qi + a * 10^(-q/10)
//
// separates two things. The SLOPE calibrates the instrument: 1.0 means the reported Phred is
// exactly right, and it is what `error(q)` applies.
//
// Never: The INTERCEPT is NOT a sequencing floor, and using it as one would add a constant to every
// base likelihood in the pipeline. The standard being measured against is a SYNTHESISED oligo,
// and oligo synthesis carries roughly one defect per 200-500 bases -- so on a real primer the
// intercept comes out at ~4e-3, which is the primer's own quality and not the instrument's.
// Measured on SRR1763769: intercept 3.9e-3 spread evenly over all 23 anchor positions (none
// polymorphic), against an independently measured 0.55% rate of one-base-short barcodes from
// failed couplings in the same oligo. Same order, same cause. It is reported as a diagnostic of
// the primer and deliberately left out of `error()`.
// Note: It only works if the "constant" bases really are constant. A pattern position that is 97%
// conserved rather than 100% contributes 3% mismatch at every quality, and the fit reads that as a
// quality-independent floor. So the counts are kept per POSITION as well, and a position whose
// mismatch rate is far above its neighbours' is dropped before fitting -- it is polymorphic, or
// the pattern is wrong about it, and either way it is not measuring the instrument.
struct QualityCalibration {
    static constexpr size_t kMaxPositions = 64;
    // [position][q] = {bases seen, mismatches} at unambiguous scored pattern positions.
    std::vector<std::array<std::array<uint64_t, 2>, 61>> by_position;
    // Summed over the positions that survived the check above. Empty until fit() runs.
    std::array<std::array<uint64_t, 2>, 61> counts{};
    std::vector<uint8_t> position_used;
    size_t positions_dropped = 0;
    // The intercept: the anchor's own defect rate, not a sequencing floor. See above.
    double quality_independent = 0.0;
    double slope = 0.0;                // 1.0 would mean the reported Phred is exactly right
    uint64_t bases = 0;
    bool fitted = false;               // false when no Q value had enough bases to fit

    QualityCalibration() : by_position(kMaxPositions), position_used(kMaxPositions, 0) {}

    void merge(const QualityCalibration& o);
    // Drops polymorphic positions, then fits by weighted least squares over the Q values with at
    // least `min_bases` observations. `max_excess` is how far above the median mismatch rate a
    // position may sit before it is treated as variable rather than miscalled.
    void fit(uint64_t min_bases = 1000, double max_excess = 5.0);
    // The calibrated sequencing error at reported Phred q: `slope * 10^(-q/10)`, or the nominal
    // rate when nothing could be fitted. The intercept is deliberately excluded -- it belongs to
    // the primer that was used as the standard, not to the instrument.
    double error(int q) const;
};

// Exact payload lengths tracked, per sample. 512 covers every current read length; anything
// longer lands in the last bin, which is labelled as the catch-all it is.
inline constexpr size_t kPayloadHistLen = 512;

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
    // Payload length after trimming, per sample, exact up to kPayloadHistLen-1 with the last bin a
    // catch-all. This is the trim's own QC: a pattern matched at the wrong offset does not fail,
    // it succeeds and leaves the payload a fixed number of bases short or long, and the length
    // distribution is where that is visible. `trimmed_bases` is what was removed.
    std::vector<std::array<uint64_t, kPayloadHistLen>> payload_len;
    uint64_t trimmed_bases = 0;
    QualityCalibration calibration;

    void merge(const CheckoutCounters& o);
};

struct CheckoutRead {
    bool ok = false;
    int sample = -1;
    std::string umi;
    std::string umi_qual;
    std::string cell;       // empty unless the pattern captured X positions
    std::string cell_qual;
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
    std::string cell;
    std::string cell_qual;
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
                                   const std::string& sample, const std::string& cell = {},
                                   const std::string& cell_qual = {});

private:
    const PatternSet& patterns_;
    CheckoutParams params_;
    CheckoutCounters counters_;
    CheckoutScratch scratch_;  // backs the single-read process() overload
};

// ---------------------------------------------------------------------------------------------
// The whole-file driver.

// `.mig` output sizing. The writer count is a RUN budget: 256 open files is already more than
// polite, and it is the number of buckets across every sample rather than within one.
inline constexpr int kMaxMigBucketBits = 8;
inline constexpr size_t kMaxMigWriters = 256;
// Bytes shared across every open bucket writer. Split, not per writer: each writer accumulates a
// block before compressing it, so a fixed per-writer block would make a finer partition cost more
// memory, which is backwards. The same rule, and the same number, as `assemble`'s partition pass.
inline constexpr size_t kMigWriterBudgetBytes = 32u << 20;

struct CheckoutRequest {
    // Stop after this many input reads; 0 reads all of them. A smoke test, never a sample: the
    // first N reads of a FASTQ are one corner of one flowcell. `subsample` is the sampler.
    uint64_t limit_reads = 0;
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
    // Resident bytes the UMI counters may hold between them before they range-partition to disk.
    // This was the last allocation in the pipeline that grew with the library rather than with the
    // chunk: ~22 B per distinct barcode in one piece, 8.8 GB at NovaSeq scale. Past the budget the
    // counters spill and everything downstream of them streams instead.
    //
    // Note: it is a whole-run budget and is divided by the number of samples, because a 96-plex
    // sheet holds 96 of these. The floor keeps a heavily multiplexed run from spilling a counter
    // small enough to be free to hold.
    //
    // Note: 0 disables spilling. The counters then grow without bound and `umi_memory_bytes` is
    // the only thing that says so, which is where this was before the partition existed.
    size_t umi_budget_bytes = size_t{1} << 30;
    // Where the partition goes. Empty puts it next to the output, in `<out_prefix>.umi_spill`, and
    // it is removed when the statistics that read it are done with it.
    std::string umi_spill_dir;

    // Write `.mig` buckets -- `<prefix><sample>.<bbb>.mig` -- instead of one FASTQ per sample.
    //
    // The reads are range-partitioned on the same key `assemble` groups by, so `assemble` reads
    // them as its own buckets and skips its partition pass entirely. Opt-in: FASTQ stays the
    // default, because it is what every aligner, every existing pipeline and `docs/downstream.rst`
    // speak, and a `.mig` file is an intermediate that only migec reads.
    //
    // Never: RANGE, on the whole barcode. A hash would split a barcode from its 1-substitution
    // neighbours across buckets and correction could never be applied locally; and the partition
    // key is the cell when there is one, exactly as in `assemble`, because a molecule is
    // sample + cell + UMI and grouping on the UMI alone merges two molecules.
    bool mig_output = false;
    // 0 chooses from the number of samples: the open-file budget is for the RUN, not per sample,
    // so a 96-plex sheet gets a couple of buckets each rather than 96 x 256 open writers.
    int mig_bucket_bits = 0;
};

// Owns a directory of spilled UMI buckets for as long as anything can still read it.
//
// Never: the spill files outlive `run_checkout`. The per-sample statistics -- histogram,
// composition, correction -- stream the partition *after* the run returns, so deleting it at the
// end of the run would leave every counter pointing at nothing. It dies with the last copy of the
// stats instead, which is the object those readers hold.
struct UmiSpillDir {
    std::string path;
    ~UmiSpillDir();
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
    std::vector<std::array<uint64_t, kPayloadHistLen>> sample_payload_len;  // and so is this
    std::vector<UmiCounts> umi_counts;   // parallel to sample_ids
    double wall_seconds = 0.0;
    double reads_per_second = 0.0;
    size_t peak_rss_bytes = 0;
    size_t umi_memory_bytes = 0;  // the part of the above that is the UMI counters
    bool umi_spilled = false;     // at least one counter went past the budget and partitioned
    std::shared_ptr<UmiSpillDir> umi_spill;  // keeps the partition readable; see UmiSpillDir
    int threads = 1;
    // `.mig` output only: how the reads were partitioned, and the files that came out. The paths
    // are what `assemble` is handed to skip its own partition pass.
    bool mig_output = false;
    int mig_bucket_bits = 0;
    std::vector<std::string> mig_paths;  // sample-major, then bucket; empty buckets are omitted
};

// Demultiplex, extract, trim and write. Throws MigecError on malformed input or an unwritable
// output path.
CheckoutStats run_checkout(const PatternSet& patterns, const CheckoutParams& params,
                           const CheckoutRequest& request);

}  // namespace migec

#endif  // MIGEC_CHECKOUT_HPP
