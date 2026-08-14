// The .mig intermediate format: every stage reads and writes this and nothing else.
//
// Why this exists, and why it looks like this:
//
//  * Raw ASCII sequence, not 2-bit packed. Packing the sequence saves 0.75 B/base but the quality
//    string is the same length and is near-incompressible, so packing buys ~13% of the record and
//    costs a pack/unpack on every read. Measured on a 2x150 amplicon block: 197 B/pair raw+zstd-1
//    vs 227 B/pair packed+zstd-1 — packing came out WORSE, because zstd finds the cross-read
//    redundancy in ASCII amplicon sequence that packing destroys.
//
//  * Column-major within a block (all seq1, then all seq2, then all qual1, then all qual2).
//    Interleaving sequence and quality per record mixes two very different symbol distributions
//    and costs zstd 10-20% on the same data.
//
//  * src_index is u64, not u32. u32 caps at 4.29e9 read pairs and a NovaSeq X run exceeds that.
//    It is the sort tiebreak, so on overflow the "identical output at 1 and 8 threads" guarantee
//    fails nondeterministically -- the worst failure mode available.
//
//  * Buckets are RANGE partitions of the key (top bits), not hash partitions. A hash sends a
//    barcode and its Hamming-1 neighbour to uncorrelated buckets, which makes UMI correction
//    impossible to apply locally; range partitioning costs the same (the packed key is near
//    uniform) and additionally makes bucket order == key order, i.e. the on-disk sort by
//    sample/cell/UMI is a property of the layout rather than a separate pass.
//
// Layout: [Header][Block]*[Footer]. A block is a zstd-1 frame over the column-major payload,
// preceded by a plaintext BlockHeader carrying the CRC32C of the *uncompressed* payload.

#ifndef MIGEC_MIG_RECORD_HPP
#define MIGEC_MIG_RECORD_HPP

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace migec {

// Bumped whenever the on-disk layout changes. A reader refuses a file it does not know -- but it
// still reads every older version it can, because the alternative is that an intermediate written
// by yesterday's build becomes unreadable rather than merely older.
//
// v2 added the BARCODE's own quality string, per record. v1 stored only `umi_minq`/`cell_minq`,
// the minimum over the barcode, which is not what a correction posterior needs: it weighs the
// reported quality AT THE POSITION THAT DIFFERS, and a minimum over the whole barcode says every
// position is as bad as the worst one. That overstates the error everywhere, which makes merges
// easier -- the wrong direction, because a wrong merge destroys a molecule and a missed one only
// inflates a count. So `refine` reading a v1 file falls back to the global rate, exactly as it
// does for a FASTQ with no QX tag, and says so rather than using a worse number silently.
inline constexpr uint16_t kMigFormatVersion = 2;
inline constexpr uint16_t kMigMinReadableVersion = 1;
inline constexpr char kMigMagic[4] = {'M', 'I', 'G', 'B'};

// Flags describe what has ALREADY been applied to the stored sequence, never what remains to be
// done. In particular kRevComp1/2 mean "this mate is stored reverse-complemented relative to the
// input FASTQ" -- assemble must never re-orient anything.
enum RecordFlags : uint16_t {
    kRevComp1   = 1u << 0,
    kRevComp2   = 1u << 1,
    kMatesSwapped = 1u << 2,  // the primary pattern matched R2, so the mates were exchanged
    kUmiHasN    = 1u << 3,
    kCellHasN   = 1u << 4,
    kCellCorrected = 1u << 5, // cell barcode was changed by a whitelist lookup
    kUmiCorrected  = 1u << 6, // umi was changed by refine
    kSingleEnd  = 1u << 7,
};

// One read (pair). Barcodes are 2-bit packed, 2 bits per base, base 0 in the low bits; an N is
// stored as A with the corresponding kUmiHasN/kCellHasN flag set, so the packed value stays a
// plain integer and the ambiguity is recorded out of band. Lengths live in the file header
// because they are constant for a whole file.
struct MigRecord {
    uint64_t cell = 0;
    uint64_t umi = 0;
    uint64_t src_index = 0;  // input order; sort tiebreak, and the determinism guarantee
    uint16_t flags = 0;
    uint8_t umi_minq = 0;   // min Phred over the umi bases, capped at 60
    uint8_t cell_minq = 0;
    std::string_view seq1, qual1, seq2, qual2;  // views into the reader's block arena
    // The barcode's own quality, one Phred+33 character per base, `umi_len`/`cell_len` long. Empty
    // on a v1 file and whenever the header says the file does not carry it. Fixed width per
    // record, because both lengths are file constants -- so it costs no length field and stays a
    // column like the others.
    std::string_view qual_umi, qual_cell;
};

// Per-file constants. Written once, in plaintext, so that a truncated file can still be
// identified. `provenance` is a JSON blob (command line, version, sample id, pattern) -- it is
// documentation, not something the reader interprets.
struct MigHeader {
    uint16_t format_version = kMigFormatVersion;
    uint8_t umi_len = 0;
    uint8_t cell_len = 0;
    uint8_t bucket_index = 0;
    uint8_t bucket_bits = 0;  // 0 == a single bucket, i.e. the whole sample in key order
    bool paired = true;
    // Whether every record carries `qual_umi`/`qual_cell` (v2 and later). A per-FILE decision
    // rather than a per-record one: the columns are fixed width, so a file either has them for
    // every record or for none, and the writer refuses a record that disagrees with the header.
    bool barcode_quality = false;
    std::string sample_id;
    std::string provenance;
    // Measured error rate per reported Phred value: quality_calibration[q] is the observed error
    // frequency for bases the instrument called q. Empty means "not measured, use 10^(-q/10)".
    // This is a per-file constant because it is estimated once, by checkout, from mismatches
    // against the constant segments of the pattern.
    std::vector<float> quality_calibration;
};

// Writes records in the order given. The caller is responsible for the ordering; the writer only
// records, in the footer, whether the file it produced is sorted (so a reader can tell a
// checkout-ordered file from a refine-sorted one without scanning it).
class MigWriter {
public:
    MigWriter(const std::string& path, const MigHeader& header, size_t block_bytes = 4u << 20);
    ~MigWriter();
    MigWriter(MigWriter&&) noexcept;
    MigWriter& operator=(MigWriter&&) noexcept;

    void write(const MigRecord& rec);
    // Flushes the pending block and writes the footer. Called by the destructor, but call it
    // explicitly if you want the error.
    void close();

    uint64_t records_written() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

// Streams a file back. Sequence views point into an internal arena that stays valid until the
// next next() call.
class MigReader {
public:
    explicit MigReader(const std::string& path);
    ~MigReader();
    MigReader(MigReader&&) noexcept;
    MigReader& operator=(MigReader&&) noexcept;

    const MigHeader& header() const;
    // Returns false at clean end of file. Throws MigecError on a truncated or corrupt block --
    // a truncated file is an error, never a short read.
    bool next(MigRecord& out);

    // Total record count from the footer, or 0 if the footer is missing.
    uint64_t records_declared() const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

}  // namespace migec

#endif  // MIGEC_MIG_RECORD_HPP
