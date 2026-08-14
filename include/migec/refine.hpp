// refine: correct the barcodes, then hand the reads on with the corrected barcode in the header.
//
// This is the stage that decides how many molecules there were. It reads checkout's tagged FASTQ,
// builds the barcode table, folds error-child barcodes into their parents, and rewrites the reads.
//
// What it holds is the barcode TABLE, never the reads: (key, count) plus this barcode's own
// evidence -- the mean error at each barcode position, and a short draft of its payload. The reads
// are streamed three times instead. The table is the same allocation checkout already carries, and
// it is bounded by the number of distinct barcodes rather than by the number of reads.
//
// When the reads carry a cell barcode the key is the CONCATENATION of cell and UMI, because a
// molecule is the whole barcode: the same UMI in two cells is two molecules, and correcting one
// against the other would merge them. Concatenating also means a 1-substitution neighbour is found
// whether the substitution landed in the cell barcode or in the UMI, which is the right
// neighbourhood -- both are real errors. Note: It is NOT a substitute for whitelisting a 10x cell
// barcode against the known list; that is a separate mechanism and is not implemented yet.
//
// Never: correction is not bucketable by a plain range partition. A range partition on the top b
// bits puts a barcode and its 1-substitution neighbours in the same bucket for every position
// EXCEPT the top b/2, and a neighbour that crosses a bucket boundary would never be found -- the
// partition would silently stop correcting. So it runs in two passes, the second on keys rotated
// past the partitioned prefix, and every pair is weighed in exactly one of them.
//
// Never: the table carries the EVIDENCE, not only the counts. A side array indexed against
// `entries()` cannot survive a partition, and dropping it would leave the bucketed run correcting
// on the count ratio alone -- which reports nothing at 1-3 reads per UMI, the regime the evidence
// exists for. `UmiCounts::carry_evidence` partitions the barcode quality and the payload draft
// with the key.

#ifndef MIGEC_REFINE_HPP
#define MIGEC_REFINE_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "migec/fastq.hpp"
#include "migec/umi_stats.hpp"
#include "migec/whitelist.hpp"

namespace migec {

struct RefineRequest {
    std::string input;       // a per-sample FASTQ written by checkout
    // ...or the `.mig` buckets `checkout --mig` wrote for one sample, in which case `input` is
    // ignored. The reads are then already partitioned on the barcode and already carry it packed,
    // so nothing here has to parse a tag -- and with a v2 file the barcode's own quality comes
    // with them, which is the evidence the posterior weighs at the position that differs.
    //
    std::vector<std::string> mig_inputs;
    // Write `.mig` buckets instead of a FASTQ. Requires `.mig` input: the buckets carry no read
    // name, and inventing one would put a name into the pipeline that no instrument produced.
    bool mig_output = false;
    std::string output_dir;
    std::string sample_id;   // taken from the BC tag when empty
    CorrectionParams correction;
    // Bases of payload kept per barcode as the draft used for payload agreement. 32 is enough to
    // tell two molecules apart -- two random sequences differ at ~24 of 32 -- and short enough
    // that the table stays small.
    int payload_width = 32;
    // Cell calling, when the reads carry a cell barcode. OrdMag, which is Cell Ranger's original
    // rule: take the 99th percentile of the top `expect_cells` barcodes by molecule count and
    // keep everything within a tenth of it. Never: EmptyDrops-style rescue of low-count cells is
    // deliberately NOT reproduced -- it is Cell Ranger's job, and pretending to match it would
    // make every comparison against their calls unreachable by construction rather than by
    // measurement.
    int expect_cells = 3000;
    // Snap cell barcodes to a list of the barcodes that were actually synthesised, before any of
    // the rest of this runs. Empty disables it.
    // Residual-FDR target for the reported MIG size threshold. Never: The threshold is REPORTED, never
    // applied: a molecule seen three times with no plausible parent is information, and cutting it
    // discards real sequence. Downstream may filter on it; refine does not.
    double target_fdr = 0.05;
    std::string cell_whitelist;
    WhitelistParams whitelist;
    // Stop early: a smoke test, never a sample. See IntakeLimit.
    IntakeLimit limit;
    // Turn the evidence off, to measure what the count ratio alone would have done.
    bool use_quality = true;
    bool use_payload = true;
    // Bytes the barcode table may hold before it range-partitions itself to disk, exactly as
    // `checkout`'s UMI counters do. Never a CLI flag: a run that needs it is one whose library
    // does not fit, which the stage can see for itself.
    //
    // Past the budget the table spills into `<out>/.refine_spill`, correction runs the same two
    // passes with the key rotated, and every table refine writes is streamed one bucket at a time.
    // Nothing else changes -- the resident case is one bucket, so there is one code path and not
    // two, and `nbuckets == 1` makes pass 2 a no-op.
    size_t table_budget_bytes = size_t{1} << 30;
    int table_bucket_bits = 8;
    // Level 1, not zlib's default 6. Measured on refine's own output: level 6 spent 1.78 s of a
    // 2.14 s run compressing 500,000 reads -- 83% of the wall clock -- against 0.34 s at level 1
    // for 21% more bytes. Read payload is close to incompressible (checkout measured 7 MB/s at
    // level 6 against 137 at level 1 on random DNA), so the extra CPU buys almost nothing, and
    // this file is an intermediate that the next stage decompresses immediately.
    int gzip_level = 1;
};

struct RefineStats {
    // `.mig` output only: the buckets that were written, in bucket order. Hand these to
    // `assemble` and it skips its partition pass, exactly as it does with checkout's.
    std::vector<std::string> mig_paths;
    uint64_t reads = 0;
    uint64_t reads_without_umi = 0;
    // True when --limit-read or --limit-umi stopped the intake: every number below then describes
    // a prefix of the file, not the library.
    bool limited = false;
    uint64_t barcodes = 0;          // distinct, before correction
    uint64_t merged = 0;            // barcodes folded into a parent
    uint64_t merged_reads = 0;
    uint64_t merged_by_payload = 0; // merges the count ratio alone would have refused
    uint64_t molecules = 0;         // distinct barcodes after correction
    double molecules_corrected = 0.0;
    double estimated_error = 0.0;
    double payload_clonality = 0.0;
    bool saturated = false;
    int umi_length = 0;
    int cell_length = 0;
    // Cell calling. All zero on a bulk library, which has no cells to call.
    uint64_t cells_observed = 0;
    uint64_t cells_called = 0;
    uint64_t molecules_in_called = 0;
    uint32_t cell_threshold = 0;   // molecules a cell needs, from OrdMag
    uint64_t knee_rank = 0;        // where the curve breaks, reported next to the threshold
    uint32_t knee_molecules = 0;
    WhitelistStats whitelist;
    // Smallest MIG size whose estimated residual false-molecule rate is at or below `target_fdr`,
    // and the rate at size 1. 0 when nothing was estimable.
    uint32_t mig_size_threshold = 0;
    double residual_fdr_at_one = 0.0;
    uint64_t suspected_residual = 0;
    // The barcode error rate read off the CHILDREN, per base per read, against the reads their
    // parents carried: a parent seen c times had c*L barcode bases to miscall, and its children
    // hold what was miscalled. Independent of the distance-1 census in `estimated_error`, which
    // counts distinct NEIGHBOURS and therefore saturates at 3L per parent as the space fills --
    // this one counts READS and does not. `<sample>.umi_errors.tsv` is the same measurement per
    // depth, and `error_at_depth` is it read where correction is near-complete rather than
    // averaged over the shallow molecules whose parents were never sequenced.
    double error_from_children = 0.0;   // library-wide, all depths; a LOWER bound at 1-3 reads/UMI
    double error_at_depth = 0.0;        // the same, restricted to parents at or above `error_depth`
    uint32_t error_depth = 0;           // the depth it was read at, 0 when no depth qualified
    // -10 log10(error_at_depth), so it can be put beside the barcode's own reported Phred.
    double error_phred = 0.0;
    std::string sample_id;
    // Bytes held by the barcode table. Reported for the same reason checkout reports its
    // counters: it is what decides whether a run fits.
    uint64_t table_bytes = 0;
    // ...and what it actually held while it ran, which is the budget rather than the library once
    // the partition fires. Reporting only the first would say a bounded run needed 40 GB;
    // reporting only the second would hide that the library is that big.
    uint64_t table_resident_bytes = 0;
    // True when the table went past `table_budget_bytes` and partitioned itself to disk. The
    // answers are the same either way; the wall clock is not, so the flag travels with them.
    bool table_spilled = false;
    // MIG size histogram after correction, power-of-two bins.
    std::vector<uint64_t> size_histogram;
    double wall_seconds = 0.0;
    // The three passes, separately, because they scale with different things and only one of them
    // threads. `table` and `rewrite` stream the reads; `correct` walks the barcode neighbourhood
    // and is the part `--threads` speeds up. Reporting one number would hide which to fix next.
    double table_seconds = 0.0;
    double correct_seconds = 0.0;
    double rewrite_seconds = 0.0;
    int threads = 1;
};

RefineStats refine(const RefineRequest& request);

}  // namespace migec

#endif  // MIGEC_REFINE_HPP
