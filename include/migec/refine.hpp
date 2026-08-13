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
// neighbourhood -- both are real errors. ⚠ It is NOT a substitute for whitelisting a 10x cell
// barcode against the known list; that is a separate mechanism and is not implemented yet.
//
// ⚠ Correction is not bucketable by a plain range partition. A range partition on the top b bits
// puts a barcode and its 1-substitution neighbours in the same bucket for every position EXCEPT
// the top b/2, and a neighbour that crosses a bucket boundary can never be found. Doing it in
// buckets needs two passes with the key rotated, so that every pair shares a bucket in at least
// one of them. Until that lands the table is held whole and its size is reported, exactly as
// checkout does with its counters.

#ifndef MIGEC_REFINE_HPP
#define MIGEC_REFINE_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "migec/umi_stats.hpp"

namespace migec {

struct RefineRequest {
    std::string input;       // a per-sample FASTQ written by checkout
    std::string output_dir;
    std::string sample_id;   // taken from the BC tag when empty
    CorrectionParams correction;
    // Bases of payload kept per barcode as the draft used for payload agreement. 32 is enough to
    // tell two molecules apart -- two random sequences differ at ~24 of 32 -- and short enough
    // that the table stays small.
    int payload_width = 32;
    // Cell calling, when the reads carry a cell barcode. OrdMag, which is Cell Ranger's original
    // rule: take the 99th percentile of the top `expect_cells` barcodes by molecule count and
    // keep everything within a tenth of it. ⛔ EmptyDrops-style rescue of low-count cells is
    // deliberately NOT reproduced -- it is Cell Ranger's job, and pretending to match it would
    // make every comparison against their calls unreachable by construction rather than by
    // measurement.
    int expect_cells = 3000;
    // Turn the evidence off, to measure what the count ratio alone would have done.
    bool use_quality = true;
    bool use_payload = true;
    int gzip_level = 6;
};

struct RefineStats {
    uint64_t reads = 0;
    uint64_t reads_without_umi = 0;
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
    std::string sample_id;
    // Bytes held by the barcode table. Reported for the same reason checkout reports its
    // counters: it is what decides whether a run fits.
    uint64_t table_bytes = 0;
    // MIG size histogram after correction, power-of-two bins.
    std::vector<uint64_t> size_histogram;
    double wall_seconds = 0.0;
};

RefineStats refine(const RefineRequest& request);

}  // namespace migec

#endif  // MIGEC_REFINE_HPP
