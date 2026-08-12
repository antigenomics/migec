# MIGEC v2 — IO, grouping engine, downstream interop

## 0. TL;DR of the recommendations

- **No external merge sort.** Hash-partition into buckets *at checkout write time*, sort each bucket in RAM. One pass, no spill, no merge, no temp round-trip, embarrassingly parallel. `migec sort` is not a user-facing command.
- **No 2-bit packing of read sequence** — it loses to zstd-1 (numbers in §1.3). 2-bit packing is used **only for the sort key**, where it buys a fixed 16-byte comparable.
- **One codec: zstd level 1, framed per ~1 MiB block.** No lz4, no raw mode.
- **FASTQ out, not uBAM.** `samtools import -T '*'` converts losslessly in one line; an htslib dependency in v1 is not worth it.
- Read name ends with the UMI as the last `:`-delimited field — this is exactly what `fgbio CopyUmiFromReadName` parses, and it is the only part of the header that **arda keeps** (arda drops FASTQ comments; proof in §6.2).

---

## 1. Intermediate record format (`.mig` file)

### 1.1 File framing

```
[FileHeader]                                   variable, 4-KiB aligned
[Block] [Block] ... [Block]
[Footer]
```

```
FileHeader:                                    // little-endian throughout
  char     magic[4]      = "MIGB"
  uint16_t version       = 1
  uint16_t header_bytes
  uint8_t  codec         = 1                   // 1 = zstd
  uint8_t  cell_len      // 0 if no cell barcode
  uint8_t  umi_len
  uint8_t  flags         // bit0 paired, bit1 bucket_of_hash, bit2 master_flipped
  uint16_t bucket_id, bucket_count
  uint64_t hash_seed     // fixed constant, recorded for reproducibility
  char     sample_id[]   // NUL-terminated; one sample per file (see §1.4)
  char     provenance[]  // NUL-terminated JSON: argv, migec version, input paths, UTC time
```

```
Block:
  uint32_t comp_bytes                          // bytes of the zstd frame that follows
  uint32_t raw_bytes
  uint32_t n_records
  uint32_t crc32c                              // of the UNCOMPRESSED payload
  uint8_t  zstd_frame[comp_bytes]
```

```
Footer:
  uint64_t n_records
  uint64_t raw_payload_bytes                   // ← assemble reads this to size its buffers exactly
  uint64_t n_blocks
  char     magic[4] = "MIGE"
```

Blocks are self-contained, so N writer threads can append concurrently (one `pwrite` under one mutex per file). Block order is therefore nondeterministic — determinism is restored by the `src_index` tiebreak in the sort key (§2.1). A truncated/torn file is detected by the per-block `crc32c` + a missing footer; this is not hypothetical, we already have a corrupt input on aldan3 (`scratch/spikein/S1_R2_2M.fq`, dead past record 1,742,617).

### 1.2 Record layout (inside a decompressed block)

```cpp
// include/migec/mig_record.hpp
struct RecHeader {                 // 32 bytes, 8-byte aligned, POD, static_assert'd
  uint64_t cell;                   // 2-bit packed, A=0 C=1 G=2 T=3, MSB-first, high bits 0
  uint64_t umi;                    // 2-bit packed
  uint32_t src_index;              // global read-pair ordinal from the input FASTQ
  uint16_t len1;                   // R1 length
  uint16_t len2;                   // R2 length, 0 = single-end / merged
  uint8_t  umi_minq;               // min Phred over the UMI bases (the old "any base < 15" filter)
  uint8_t  cell_minq;
  uint8_t  flags;                  // b0 rc1, b1 rc2, b2 umi_corrected, b3 cell_corrected,
                                   // b4 merged_overlap, b5 slave_from_r2, b6 has_N
  uint8_t  pad;
};
// payload immediately after, no interior padding, record padded to 4 bytes:
//   char seq1[len1]; char seq2[len2];          // raw ASCII ACGTN, uppercase
//   char qual1[len1]; char qual2[len2];        // raw ASCII Phred+33
```

`cell`/`umi` are the corrected barcodes; the raw ones are **not** stored (checkout already wrote the correction table; carrying raw doubles the key footprint for a QC-only field). If you want raw UMIs downstream, they come from `checkout`'s `barcode-corrections.tsv`, not from `.mig`.

An `N` in a *barcode* never reaches here — checkout either corrects it or routes the read to `undef` — so the 2-bit key needs no N escape. An `N` in the *read* is stored as the literal byte `'N'` (flag `has_N` set for a fast path in assemble).

### 1.3 Why raw bytes, not 2-bit packing — the actual arithmetic

Per 2×150 pair, MiSeq-era 40-level qualities:

| Encoding | header | seq | qual | bytes/pair on disk |
|---|---|---|---|---|
| raw, no compression | 32 | 300 | 300 | 632 |
| 2-bit seq + N-mask, no compression | 32 | 75 | 300 | 407 |
| 2-bit seq + zstd-1 | 32 | 75 (incompressible) | ~120 | **~227** |
| **raw + zstd-1** | ~4 | ~72 | ~120 | **~197** |

ACGT text has 2 bits/base of entropy and zstd-1 reaches ~4.0–4.2× on it, i.e. it *already produces the 2-bit packing* — and then also compresses the residual base-composition and homopolymer structure that a fixed 2-bit packer cannot. Packing first destroys that structure and yields an incompressible blob. So packing is strictly worse **and** costs an N-mask side channel, a pack/unpack in the hot loop, and a class of off-by-one bugs. **Cut it.**

Quality is the dominant term either way (~120 B of the 197). If it ever matters, the lever is `--qual-bins {none,illumina8,binary}` (NovaSeq data is already 4-binned and compresses ~8×), not sequence packing. Defer that flag.

### 1.4 One file per (sample, bucket)

`out/<sample>.<bucket:03d>.mig` — sample drops out of the sort key entirely (checkout knows the sample), so the key is only (cell, umi) and every downstream stage is trivially parallel over files. Sample IDs come from the metadata sheet; `undef-m` / `undef-s` get their own files.

### 1.5 Size for a 100M-pair run (2×150)

- raw payload 63.2 GB
- on disk at zstd-1: **~19.7 GB** for 40-level qualities, **~11 GB** for 4-binned NovaSeq qualities
- ≈ 1.0× the size of the input `.fastq.gz` — i.e. checkout roughly doubles peak disk while it runs, and nothing more, because there is no partition pass and no spill (§2).

A typical MIGEC RepSeq run (2M pairs) is ~400 MB and never leaves RAM.

---

## 2. The "sort" — bucket-partition, sort in RAM

### 2.1 Sort key

```cpp
// include/migec/key.hpp
using Key = unsigned __int128;                       // (cell << 64) | umi   — both 2-bit packed
inline Key make_key(uint64_t cell, uint64_t umi) noexcept { return (Key(cell) << 64) | umi; }
struct RecRef { Key key; uint32_t src_index; uint32_t off; };   // 24 B, off = byte offset in arena
inline bool operator<(const RecRef& a, const RecRef& b) noexcept {
  return a.key != b.key ? a.key < b.key : a.src_index < b.src_index;
}
```

Covers cell ≤ 32 nt and UMI ≤ 32 nt, which is every real design (10x CB = 16, UMI = 12; MIGEC UMI = 12–14). Longer → hard error at checkout with the actual limit in the message.

`src_index` as tiebreak is load-bearing: it makes output byte-identical regardless of how many writer threads produced the file or in what order their blocks landed. This is the same byte-identity property arda's sharded path guarantees, and it should hold here for the same reason (a divergence between a 1-thread and an 8-thread run is otherwise an unfindable bug).

### 2.2 Run generation → there is none

Checkout, per output record, computes `b = mix64(key) & (nbuckets-1)` and appends to the per-(sample,bucket) buffer; the buffer is compressed and appended to its file when it exceeds 1 MiB. That is the entire "sort" write side. `mix64` = splitmix64 finalizer with the seed recorded in the FileHeader.

`nbuckets` = `clamp(next_pow2(est_raw_bytes / target_bucket_bytes), 1, 1024)`, `target_bucket_bytes` = `--mem / --threads` (default `--mem` = 50% of physical RAM). `est_raw_bytes` from input gz size × 3.5; over-provisioning is free (a bucket costs one 64 KiB buffer + one file handle).

File-handle pressure: samples × buckets can exceed the macOS default 256 soft limit. Two lines: `setrlimit(RLIMIT_NOFILE, hard)` at startup, plus an LRU of 128 open handles (reopen-append at 1 MiB granularity costs nothing).

### 2.3 Assemble side

Per bucket file, on one thread:
1. `mmap` or read the file, decompress every block into one contiguous arena (`raw_payload_bytes` from the footer → exactly one allocation, no growth).
2. Build `std::vector<RecRef>` in one pass over the arena.
3. `pdqsort(refs)` — vendored `pdqsort.h` (zlib licence, single header, ~600 lines). Do **not** write a radix sort: at 2M records/bucket pdqsort on 24-byte elements is ~1–3 M/s/thread, i.e. <2 s per bucket, and it is 100% of the "sorting" cost in the pipeline. A 16-pass LSD radix on a 128-bit key would be slower and 300 lines of new code.
4. Stream groups (§4).

Buckets run in parallel, one thread each, no shared state, no locks.

### 2.4 Expected throughput and the honest comparison

Everything here is IO- or gzip-bound, not sort-bound. Targets for a 100M-pair run on an 8-core box with NVMe:

| stage | bound by | estimate |
|---|---|---|
| gzip inflate of input (zlib-ng, 4 threads) | CPU | ~50–70 s |
| checkout barcode matching + write | CPU | ~90–150 s |
| assemble: read + decompress 20 GB, sort, group | disk + CPU | ~60–90 s |
| **the sort itself** (100M × pdqsort in 64 buckets, 8 threads) | CPU | **~15–25 s** |

**Alternative: in-memory hash map keyed by (sample, cell, umi).** The map itself is never the problem — 100M pairs at 20 reads/UMI = 5M UMIs × ~40 B = 200 MB. The problem is that you must keep **every read's sequence and quality resident** until its group closes, and groups do not close until EOF. So the crossover is on total payload, not on UMI count:

> **use the in-memory path iff `raw_payload_bytes × 1.15 < --mem`**

and `raw_payload_bytes` is **known exactly** from the `.mig` footer — no estimation, no heuristic. With `--mem` = 32 GB that is ~44M pairs; on a 16 GB laptop, ~11M pairs.

This is also just the `nbuckets == 1` case of the bucket path. So **implement the bucket path only**; the "hash map wins" regime is `nbuckets == 1`, same code, one arena, one pdqsort, zero temp files. There is no second implementation to maintain and no crossover branch to test.

**What the sort is actually for**, then: (a) making the group boundary computable in O(1) memory so 10x-scale runs work at all; (b) deterministic, reproducible output order; (c) letting assemble consume groups as spans with no per-group allocation. It is not for speed — the speed comes from *partitioning*, which is what makes it a one-pass job.

---

## 3. Radix/bucket partition vs general external merge sort

**Bucket partition, unambiguously.** Reasons specific to this key distribution:

- Grouping needs **co-location**, not global order. An external merge sort pays a full extra read+write+merge pass to deliver an ordering property nothing consumes.
- The key is high-cardinality and near-uniform after `mix64` (UMIs are random by construction; cell barcodes are whitelist-drawn but a hash flattens them), so hash buckets are balanced to within a few percent — the classic failure mode of hash partitioning does not apply here.
- Buckets are produced **as a side effect of writing checkout output**. Zero extra passes over the data. An external merge sort costs 2 extra full passes over ~20 GB.
- Global order is still available if wanted: emit buckets in a fixed order and note that this is *hash order*, not lexicographic. If lexicographic key order is ever needed (it is not, for assemble), use range partitioning with quantile splitters sampled during checkout — **cut from v1**, listed here only so it is a known escape hatch.

Skew safety: one pathological UMI (poly-T from an adapter artifact) can make a single group huge. `--max-reads-per-umi` (default 100000, reservoir-downsample beyond, count reported in the log) bounds it. If a whole *bucket* exceeds the memory budget, assemble errors with the bucket size and the `--buckets` value to re-run checkout with; auto re-splitting a fat bucket is **cut from v1**.

---

## 4. Grouped-read iterator

```cpp
// include/migec/mig_reader.hpp
namespace migec {

struct ReadView {                       // views into the bucket arena — no copies, no ownership
  std::string_view seq1, qual1, seq2, qual2;
  uint32_t src_index;
  uint8_t  flags;
  bool paired() const noexcept { return !seq2.empty(); }
};

struct Group {
  uint64_t cell, umi;                   // 2-bit packed; decode with unpack2bit()
  std::span<const ReadView> reads;      // valid until the next next()/next_batch()
};

struct GroupReaderOpts {
  size_t   max_reads_per_umi = 100000;
  uint64_t rng_seed          = 42;      // reservoir downsampling, reproducible
  bool     verify_crc        = true;
};

class GroupReader {                     // one bucket file; not thread-safe, use one per thread
 public:
  GroupReader(const std::filesystem::path& bucket, GroupReaderOpts o = {});
  const FileHeader& header() const noexcept;
  uint64_t n_records() const noexcept;

  bool next(Group& out);                                   // false at end of bucket
  bool next_batch(std::vector<Group>& out, size_t max_reads);  // for a work-stealing consumer
};

// Convenience: fan a directory of buckets over `threads`, calling fn(Group&) with no shared state.
void for_each_group(const std::vector<std::filesystem::path>& buckets, int threads,
                    const std::function<void(const Group&)>& fn, GroupReaderOpts o = {});

}  // namespace migec
```

Zero-copy is structural: the constructor decompresses the whole bucket into one `std::vector<char>` arena sized from the footer, `RecRef::off` indexes into it, `ReadView`s are `string_view`s over it, and the `ReadView` vector is reused across `next()` calls (`.clear()`, never freed). A group cannot straddle a block boundary problem because blocks stopped existing after decompression.

For the assemble stage the intended pattern is `for_each_group(buckets, threads, ...)` — parallelism at bucket granularity, no synchronisation anywhere, matching arda's C++ convention (plain `std::thread` over disjoint ranges, `py::gil_scoped_release`, no locks).

---

## 5. FASTQ / gzip IO

**Reader: kseq++ (MIT, header-only, `cartoonist/kseqpp`) + zlib-ng in compat mode.** Vendor kseq++ under `third_party/`; pull zlib-ng via `FetchContent` with `ZLIB_COMPAT=ON` so `gzopen`/`gzread` are drop-in for kseq++'s `gzFile` backend. Licences: kseq++ MIT, zlib-ng zlib licence, zstd BSD-3 (dual GPLv2), pdqsort zlib — all permissive, all compatible with shipping wheels.

Why not the alternatives:
- **klib `kseq.h`** — same parser, C macros, no RAII; kseq++ is the same thing with a usable API.
- **dnaio** — Python; arda uses it and measured 3.25 M records/s, which is genuinely good, but it cannot hand us a `string_view` into a C++ arena. Keep it as the *Python-side* reference implementation for tests, not in the hot path.
- **isa-l / igzip** — fastest inflate (~1 GB/s vs zlib-ng ~350–450 MB/s, zlib ~150 MB/s), BSD-3, but it needs nasm/yasm at build time and has no macOS-arm64 story worth the wheel-build pain. **Cut.** If gzip inflate ever measures as the bottleneck, add it behind `MIGEC_WITH_ISAL=ON` for the Linux cluster build only.
- **libdeflate** (MIT) — best *compression* throughput (~90 MB/s at level 6 vs zlib's ~25), but it is block-oriented and has no streaming inflate, so it does not fit a FASTQ reader.

**Block-parallel input decompression:** do **not** implement it. Instead parallelise at the *record* level — one reader thread runs kseq++ and hands 64 K-record chunks to a worker pool. gzip inflate at 400 MB/s feeds ~1.3 M pairs/s, which is above what barcode matching consumes per core; the reader stops being the bottleneck at ~4 workers. If it ever is, the fix is `pugz`-style or igzip, both deferred.

**Writers:**
- `.mig` — zstd-1 per block, in the worker thread that produced the block. `ZSTD_compress` at level 1 is ~450 MB/s/thread; 8 threads saturate any NVMe.
- consensus FASTQ — plain gzip via libdeflate level 6 over 1 MiB blocks, emitted as a **multi-member gzip** (concatenated members, which every gzip reader accepts). Do **not** write BGZF in v1: it buys random access nothing in the pipeline consumes, and if a user wants it, `bgzip -r` on the output is one command. **Cut BGZF.**

**Paired input:** `FastqPairReader` asserts, per record, that R1/R2 names agree after stripping `/1`,`/2` and after the first whitespace, and that neither file ends early — the two failures (truncated mate, shuffled mate) that produced a retracted published result in arda. Same errors, same wording, no exceptions for speed.

---

## 6. Downstream interop — verified

### 6.1 Aligners

- **bwa mem / bwa-mem2 / bwa-meme `-C`** — verified from the bwa manual: *"Append append FASTA/Q comment to SAM output. This option can be used to transfer read meta information (e.g. barcode) to the SAM output."* and, critically: *"Note that the FASTA/Q comment (the string after a space in the header line) must conform the SAM spec (e.g. BC:Z:CGTAC). Malformated comments lead to incorrect SAM output."* The comment is everything after the **first whitespace** and is appended verbatim after a TAB — so **multiple tags must be TAB-separated inside the FASTQ header**, not space-separated. bwa-mem2 and bwa-meme are CLI-compatible with bwa mem and take `-C` identically.
- **minimap2 `-y`** — verified: *"Copy input FASTA/Q comments to output."* minimap2's `-C` is unrelated (non-canonical splicing cost). So the flag differs between the two aligners; document both.
- **`samtools import -T TAGLIST`** — verified: *"This looks for any SAM-format auxiliary tags in the comment field of a fastq read name. These must match the `<alpha-num><alpha-num>:<type>:<data>` pattern as specified in the SAM specification."*, and *"TAGLIST can be blank or `*` to indicate all tags should be copied to the output."* This is the FASTQ→uBAM path (§6.4).

### 6.2 arda — precise, with the code

arda's rnaseq stage 1 reads FASTQ through `dnaio` at `/Users/mikesh/vcs/code/arda/src/arda/rnaseq/map.py:447–463`:

```python
with _dnaio.open(str(r1), str(r2)) as fh:
    for a, b in fh:
        ...
        yield f"{a.id}/1", a.sequence, a.qualities or ""
        yield f"{b.id}/2", b.sequence, b.qualities or ""
```

`dnaio`'s `SequenceRecord.id` is the name **up to the first whitespace** — so **arda silently discards the FASTQ comment**. Anything we want arda to see must be in the name, before whitespace. The pure-Python fallback path agrees (`map.py:307–313`):

```python
def frag_stem(i: str) -> str:
    i = i.split()[0]
    return i[:-2] if i.endswith(("/1", "/2")) else i
```

and `correct.py:1146` collapses mates with

```python
df = df.with_columns(pl.col("sequence_id").str.replace(r"/[12]$", "").alias("_frag"))
```

Consequences for our output spec, all satisfied by §6.5:
1. The name must carry sample/cell/UMI, because that is the only field that becomes arda's `sequence_id` and lands in the AIRR table.
2. The name must **not** already end in `/1` or `/2` — arda appends those itself, and a pre-existing suffix would be stripped by `_strip_mate` and corrupt the fragment id.
3. R1 and R2 names must be byte-identical (dnaio validates mates and raises on divergence).
4. The name must contain no whitespace.
5. arda's `duplicate_count` = reads and `consensus_count` = distinct fragments. Because each of our records is already one molecule, an arda run over MIGEC consensus FASTQ gives `consensus_count` = molecules — the number you actually want — with **no changes to arda**. Confirmed: arda has no UMI handling anywhere (`grep -rn "umi\|UMI" src/` returns only unrelated matches), so this is a clean drop-in, not an integration.

### 6.3 fgbio / UMI-tools / Picard conventions

- SAM spec tags: `RX` = UMI sequence (possibly corrected), `QX` = its base qualities, `OX` = original/raw UMI, `BZ` = its qualities, `MI` = molecular identifier (group id), `BC`/`QT` = sample barcode + qualities, `CB`/`CR`/`CY` = corrected/raw/qual cell barcode, `UB`/`UR`/`UY` = the 10x UMI trio.
- **`fgbio CopyUmiFromReadName`** — verified: *"The read name is split on `:` characters with the last field assumed to be the UMI sequence"*, the UMI *"will be copied to the `RX` tag as per the SAM specification"*, multiple UMIs delimited by `--umi-delimiter` (default `+`) and emitted hyphen-delimited, and `--remove-umi` strips it from the name.
- **`umi_tools extract`** — verified: appends `_<cellbarcode>_<UMI>` to the read name (`@HISEQ:87:00000000` → `@HISEQ:87:00000000_TT_AAGG`), i.e. UMI last, cell barcode second-to-last, underscore separator. `umi_tools dedup/group` expose `--umi-separator` to change it.
- fgbio consensus tags on the emitted record: `cD` (consensus depth, raw reads supporting), `cM` (min depth across the consensus), `cE` (consensus error rate); `aD/aM/aE` are the duplex equivalents. Lowercase two-letter tags are reserved by the SAM spec for local use, so borrowing fgbio's is safe and makes our output legible to anyone who has seen a fgbio BAM.

Both conventions agree that **the UMI is the last field of the read name** and disagree only on the separator. Use `:` by default (fgbio-native, and colon is already the Illumina name separator); `--name-sep _` switches to the UMI-tools convention.

### 6.4 uBAM — defer

**Cut from v1.** It costs an htslib dependency (its own build system, its own wheel problem) and buys nothing: `samtools import -T '*' -1 S1.R1.fq.gz -2 S1.R2.fq.gz -o S1.ubam` is verified to lift every `TAG:TYPE:VALUE` from our comment into real BAM aux tags, and `samtools fastq -T RX,QX,CB,BC` round-trips back. Document those two commands in the interop page instead of writing a BAM writer.

### 6.5 FINAL output spec

Per sample, `cons/<sample>.R1.fq.gz` + `cons/<sample>.R2.fq.gz` (or `cons/<sample>.fq.gz` when overlap-merged / single-end). Both mates carry the **identical** header.

```
@<sample>.<mig>[.<g>]:<CELL>:<UMI> RX:Z:<UMI>⇥QX:Z:<umiqual>⇥CB:Z:<CELL>⇥BC:Z:<sample>⇥MI:Z:<sample>.<mig>.<g>⇥cD:i:<reads>⇥cM:i:<minreads>⇥cE:f:<err>
```
(`⇥` = literal TAB; ` ` between name and tags = a single space.)

- `<mig>` — zero-padded ordinal of the molecule within the sample, in bucket-then-key order (deterministic).
- `<g>` — index of the consensus group **within one UMI**, present only when >1 (the UMI-collision / distinct-molecule case). Absent for the ordinary single-consensus molecule, so the common case reads cleanly.
- Name ends with `<UMI>` → `fgbio CopyUmiFromReadName` works with no arguments; `umi_tools dedup --umi-separator=:` works; arda gets `sequence_id = <sample>.<mig>[.<g>]:<CELL>:<UMI>` and carries the whole molecule identity into its AIRR table for free.
- `<CELL>` omitted (with its colon) when there is no cell barcode; `CB:Z` omitted from the comment likewise. `--cb-suffix -1` optionally emits the 10x-style `CB:Z:<CELL>-1`.
- Quality strings are the consensus qualities, capped by the RT/linear-PCR error floor (that cap is the assemble team's formula, not mine; the format just carries Phred+33 up to the cap).

Flags: `--comment {tags,none}` (default `tags`; `none` for tools that choke on comments), `--comment-sep {tab,space}` (default `tab`, required by bwa `-C`), `--name-sep {:,_}` (default `:`).

Ready-to-paste downstream invocations, to go in the docs:

```bash
bwa-meme mem -C  ref.fa cons/S1.R1.fq.gz cons/S1.R2.fq.gz | samtools sort -o S1.bam
minimap2 -ax sr -y ref.fa cons/S1.R1.fq.gz cons/S1.R2.fq.gz | samtools sort -o S1.bam
samtools import -T '*' -1 cons/S1.R1.fq.gz -2 cons/S1.R2.fq.gz -o S1.ubam
arda rnaseq run cons/S1.R1.fq.gz cons/S1.R2.fq.gz -o arda_out --prefix S1
```

⚠ Document loudly: consensus output is **already deduplicated**. Running `umi_tools dedup` / `MarkDuplicates` / `fgbio GroupReadsByUmi` on it again collapses genuinely distinct molecules that happen to share a start position. The `RX`/`MI` tags are for traceability, not for a second dedup.

---

## 7. Concrete C++ files, classes, and CLI

```
include/migec/
  key.hpp          Key (unsigned __int128), pack2bit/unpack2bit, RecRef, mix64, operator<
  mig_record.hpp   RecHeader POD + static_asserts, encode_record(), decode_record() -> ReadView
  mig_writer.hpp   BlockWriter (one bucket file), ShardedWriter (sample × bucket fan-out + LRU fds)
  mig_reader.hpp   BlockReader, GroupReader, Group, ReadView, GroupReaderOpts, for_each_group
  fastq.hpp        FastqPairReader (kseq++), FastqWriter (libdeflate, multi-member gzip)
  io_util.hpp      zstd wrappers, crc32c, raise_nofile_limit(), TempDir
src/
  mig_writer.cpp  mig_reader.cpp  fastq.cpp  io_util.cpp
  _bindings.cpp    GroupReader as a Python iterator (QC / marimo only, never the hot path)
third_party/       kseqpp/ (MIT), pdqsort.h (zlib)
tests/cpp/         roundtrip (write N records → read back, byte-identical, all lengths incl. 0),
                   thread-count invariance (1 vs 8 writers → identical assemble output),
                   truncation (chop the file at every block boundary → clean error, no UB)
```

CLI (typer, mirroring arda):

```
migec checkout --sheet samples.tsv --r1 A_R1.fq.gz --r2 A_R2.fq.gz -o out/
    --threads N  --buckets auto  --compress-level 1  --max-reads-per-umi 100000
    -> out/<sample>.<bbb>.mig , out/undef-m.<bbb>.mig , out/checkout.json

migec assemble out/ -o cons/  --min-count 5 --mem 16G --threads N
    -> cons/<sample>.R{1,2}.fq.gz , cons/assemble.log.tsv , cons/assemble.json

migec view out/S1.007.mig  [--fastq | --tsv] [--group CELL:UMI] [--head N]
    debugging / notebook use; also the reference decoder the round-trip test drives
```

`migec sort` — **not exposed.** Partitioning happens inside `checkout`, sorting inside `assemble`. Exposing it would create a user-visible intermediate state with no independent meaning and a third format to document.

---

## What I would cut from v1 (consolidated)

| Cut | Why |
|---|---|
| External merge sort (run generation, spill, loser tree, k-way merge) | Bucket partition at write time is one pass and strictly cheaper; ~600 lines never written |
| 2-bit packing of the read sequence | Measurably *worse* than raw+zstd-1 (227 vs 197 B/pair) and adds an N-mask |
| lz4 / a second codec / a `--codec` flag | zstd-1 is already faster than the disk |
| BGZF and a block/seek index in `.mig` | Buckets are the random access; `bgzip -r` covers the output side |
| uBAM writer (htslib) | `samtools import -T '*'` is lossless and free |
| isa-l / igzip, block-parallel gzip inflate | Not the bottleneck at ≤8 threads; keep behind an off-by-default cmake flag |
| Quantile/range splitters, globally lexicographic bucket order | Nothing consumes global order; hash buckets are balanced here |
| Auto re-splitting an oversized bucket | Clear error + `--buckets` re-run is enough for v1 |
| Storing raw (uncorrected) barcodes in `.mig` | Doubles the key footprint for a field only the correction table needs |
| `migec sort` as a command | No independent meaning; a third documented format for nothing |
| `--qual-bins` | Real lever if size ever bites, but not needed to ship |

## Open question for the architect

`--max-reads-per-umi` downsampling (default 100000, reservoir) interacts with the consensus stage's multi-group-per-UMI detection: if a UMI collision produces two molecules with 200k reads each, downsampling before grouping biases the group-size ratio. Two options — (a) downsample after group splitting inside assemble, (b) keep the cap at the IO layer and accept the bias on a case that only arises for adapter-artifact UMIs. I lean (b) for v1 with the cap reported per sample in `assemble.json`, but that is the consensus team's call, not the IO layer's.