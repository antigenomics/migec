// assemble: reads carrying the same UMI are reads of one molecule, so collapse them.
//
// The input is checkout's tagged FASTQ (the UMI in RX:Z:/QX:Z:), the output is ordinary FASTQ,
// one record per molecule, with the barcode still in the header.
//
// Grouping needs the reads of a molecule together, and there is no hash map keyed by barcode
// anywhere in this pipeline -- at NovaSeq scale that is 19 GB. Instead the reads are RANGE
// partitioned on the packed UMI into `.mig` buckets, and each bucket is sorted in RAM on its own.
// Range, never hash: a barcode and its 1-substitution neighbours must land in the same bucket or
// UMI correction cannot be applied locally, and the two halves of a split molecule each look like
// a well-formed MIG, so nothing downstream detects it. Bucket order is also key order, so the
// output is sorted by UMI for free.
//
// One bucket is resident at a time, which is what bounds the memory; the bucket count comes from
// the input size rather than from a flag.

#ifndef MIGEC_ASSEMBLE_HPP
#define MIGEC_ASSEMBLE_HPP

#include <cstdint>
#include <string>
#include <vector>

#include "migec/consensus.hpp"
#include "migec/umi_stats.hpp"

namespace migec {

// Resident bytes a bucket is allowed before the partition is split further. 1 GiB holds roughly
// 4 M 150 nt reads with their qualities and record overhead, which is a bucket that sorts in
// under a second and leaves room for the rest of a laptop.
inline constexpr uint64_t kBucketBudgetBytes = 1ull << 30;
inline constexpr int kMaxBucketBits = 8;  // 256 open temp files is already more than polite
// Bytes shared across ALL open bucket writers while partitioning. Each writer accumulates a block
// before compressing it, so a fixed per-writer block size would make pass 1 cost grow with the
// bucket count -- which is backwards, since more buckets exist precisely to use less memory. The
// budget is split instead, so the partition pass costs the same whatever it is cut into.
inline constexpr uint64_t kWriterBudgetBytes = 32ull << 20;
inline constexpr size_t kMinBlockBytes = 256u << 10;
inline constexpr size_t kMaxBlockBytes = 4u << 20;   // the writer's own default

// Reads of one barcode that enter the consensus. Past this the column posterior has long since
// saturated -- the 10,001st read moves no call and no emitted quality, because the quality is
// capped by the RT floor at Q40 and a hundred reads already clear it -- while the group still
// costs time and resident memory proportional to its size. 10x caps the same way and says why:
// "Very high coverage (greater than 10,000 reads) of transcripts can be problematic because it
// degrades computational performance and adds little information."
//
// Never: the cap applies to the reads that are CONSENSED, never to the reads that are COUNTED.
// `reads` in the table and `cD` in the FASTQ stay the true depth of the molecule, because the
// count is the other half of what this pipeline produces and capping it would silently flatten
// the abundance of exactly the most-amplified molecules.
inline constexpr size_t kMaxReadsPerGroup = 10000;

struct AssembleRequest {
    std::string input;       // a per-sample FASTQ written by checkout
    std::string output_dir;
    std::string sample_id;   // taken from the BC tag when empty
    ConsensusParams consensus;
    int gzip_level = 6;
    // Molecules below this are still assembled and still written -- a molecule seen three times
    // is information. Raise it only with a reason.
    uint32_t min_reads = 1;
    // 0 means "choose from the input size".
    int bucket_bits = 0;
};

struct AssembleStats {
    uint64_t reads = 0;
    uint64_t reads_without_umi = 0;
    uint64_t groups = 0;      // distinct UMIs
    uint64_t molecules = 0;   // consensuses emitted; exceeds `groups` when a group splits
    uint64_t groups_split = 0;
    // Groups whose reads did not all reach each other: contig mode only, and the number that says
    // whether the library was random-primed enough to need it.
    uint64_t groups_fragmented = 0;
    uint64_t contigs = 0;
    uint64_t reads_dropped = 0;  // in groups below min_reads
    // Groups over kMaxReadsPerGroup, and the reads in them that did not enter the consensus. Both
    // are still counted as reads of their molecule.
    uint64_t groups_capped = 0;
    uint64_t reads_over_cap = 0;
    int umi_length = 0;
    int cell_length = 0;
    int buckets = 0;
    std::string sample_id;
    // MIG size histogram, power-of-two bins: [0] is 1 read, [1] is 2-3, [2] is 4-7, ...
    std::vector<uint64_t> size_histogram;
    // Mean emitted Phred, and the mean posterior error the consensus itself achieved before the
    // RT floor was added. The gap between them is what the chemistry costs.
    double mean_quality = 0.0;
    double mean_consensus_error = 0.0;
    // Counting mode only: the share of a group's reads that carried the sequence that was emitted.
    // Well below 1 means the reads of a molecule disagree, which is what the full path resolves
    // per column and this one resolves by majority vote over whole strings.
    double mean_support = 0.0;
    // The birthday arithmetic over the barcodes this run actually saw. A short UMI cannot tag
    // every input molecule distinctly by design, so a group is EXPECTED to hold more than one
    // molecule -- `expected_molecules_per_group` says how many, and `molecules / groups` says how
    // many were recovered. Contig assembly is what suffers first: two fragments of two different
    // molecules on one barcode have no sequence in common and are indistinguishable from two
    // fragments of one.
    BarcodeSpace space;
    double expected_molecules_per_group = 1.0;
    double wall_seconds = 0.0;
    double partition_seconds = 0.0;
};

AssembleStats assemble(const AssembleRequest& request);

}  // namespace migec

#endif  // MIGEC_ASSEMBLE_HPP
