# checkout / barcode-extract — design

## 0. One-line summary

`migec checkout` compiles a per-sample **pattern program** (explicit named captures, MIGEC grammar accepted as a legacy dialect), matches it with a **bit-parallel matcher scored by a quality-aware log-likelihood ratio**, corrects sample/cell barcodes against whitelists with a Cell-Ranger posterior, and writes a **compact sortable binary record stream** (`.mgc`) plus, on demand, SAM-tagged FASTQ.

---

## 1. Barcode pattern grammar

### Decision

Move to **explicit named captures**; keep MIGEC's uppercase/lowercase only inside an auto-translated legacy dialect. Reason: `N` is irreconcilable across the three conventions (IUPAC "any base" vs MIGEC "UMI" vs umi_tools "UMI"); making every extracted field an explicit parenthesised capture removes the ambiguity, and named+numbered tags are what makes dual barcodes composable (`UMI1`/`UMI2` from two mates concatenate by declaration order — no "master/slave" vocabulary needed).

### Grammar (EBNF)

```
pattern   ::= elem+
elem      ::= anchor | capture | literal | gap
anchor    ::= '^'                          # read start (offset must be 0)
            | '$'                          # read end
capture   ::= '(' TAG (':' LEN)? ')'
TAG       ::= ('UMI'|'CB'|'SB'|'R'|'_') DIGIT?      # R = payload, _ = discard
LEN       ::= INT                          # required except for R and _
literal   ::= IUPAC+ ('{e<=' INT '}')?     # ACGTRYSWKMBDHVN, case-insensitive
gap       ::= '.{' INT ',' INT '}'         # variable uncaptured spacer
```

* `N` in a literal = IUPAC any-base, **uncaptured** (a spacer). UMI is *only* `(UMI:k)`.
* Case in a literal sets the mismatch weight: **UPPERCASE = hard (w=1.0), lowercase = soft (w=`--soft-weight`, default 0.5)**. Only affects the mismatch penalty (§2); it is a weight, not a boolean gate. This preserves MIGEC's intent without its dangling-else bug.
* Exactly ≤1 `(R)` per stream; if absent, `(R)` is implicit after the last element.
* Same-tag captures **concatenate in declaration order across streams** (order `r1, r2, i1, i2`), overridable with `--umi-order`, `--cb-order`.
* Streams: `--p1/--p2/--pi1/--pi2` (or sheet columns), one pattern each.

### Legacy dialect (`--pattern-dialect migec`, auto-detected for headerless `barcodes.txt`)

`N|n` → `(UMI:1)` (adjacent runs merged); uppercase run → hard literal; lowercase run → soft literal; whole pattern implicitly `^`-unanchored (searched anywhere, as before); trailing implicit `(R)`.

### Worked examples

| # | design | pattern |
|---|---|---|
| 1 | MIGEC interleaved UMI, master in R1 (legacy `TAGCGTNNNNtNNNNtNNNNctcgtat`) | `--p1 "TAGCGT(UMI:4)t(UMI:4)t(UMI:4)ctcgtat(R)"` |
| 2 | 10x 3' v3 (CB16+UMI12, positional) | `--p1 "^(CB:16)(UMI:12)" --p2 "(R)"` |
| 3 | Illumina dual index i5/i7 | `--pi1 "^(SB1:10)" --pi2 "^(SB2:10)"` (sample = `SB1+SB2` → whitelist/sheet) |
| 4 | Duplex-seq (Kennedy lab): 12nt tag + 5nt spacer on both mates | `--p1 "^(UMI1:12)AGTGGT(R)" --p2 "^(UMI2:12)AGTGGT(R)" --umi-canonical duplex` |
| 5 | MAGERI / SMART (internal exp1): UMI + SMART adapter in R2 | `--p2 "^(UMI:12)GTGGTATCAACGCAGAG(R)" --p1 "(R)" --rc-mask 0:1` |
| 6 | 10x 5' with variable TSO drift | `--p1 "^(CB:16)(UMI:10).{0,3}TTTCTTATATGGG(R)"` |

Sample identity comes from **either** the row's literal pattern (MIGEC: demux = which row matched) **or** an `(SB)` capture + whitelist (Illumina: demux = whitelist lookup). Both supported simultaneously.

---

## 2. Matching algorithm

Compile to a `PatternProgram` and pick one of three engines at compile time:

1. **Positional** — pattern fully `^`-anchored, no gaps, no literals with tolerance (10x). Direct slice, O(1)/read. Distinct code path; this is the throughput case.
2. **Anchored-with-drift** — anchored, offsets scanned over `[-d,+d]`, `--anchor-drift` (default 0; 2 for chemistries with a leading spacer).
3. **Seed-and-extend** — take the longest hard-literal run `L*` (require `L* ≥ 5`; else the pattern is rejected at compile time as unsearchable). Candidate offsets from a **bit-parallel scan** of that run over the read, then full-pattern scoring at each candidate.

**Bit-parallel matcher, not seqtree.** Patterns are ≤64 bases → one 64-bit word. Precompute `Peq[c]`, bit *i* set iff read char *c* ∈ IUPAC set of pattern position *i*; run Baeza-Yates–Gonnet shift-or with mismatch counting (no indels) or Myers' bit-vector edit-distance (indels enabled). One word-op per read base. seqtree indexes a *reference set* — the wrong shape for "one short pattern in one short read"; it is used only for whitelists at ≥2 edits (§4) and by the later refine/assemble stages. **Do not add an index here.**

### Accept/reject: quality-aware LLR (replaces the good/bad-mismatch rule)

For a candidate offset *o*, over literal positions *i* only (captures unscored), with read base `b_i`, Phred `q_i` at read position `o+i` (quality indexed at the *match* offset — this is what the original got wrong), error `e_i = 10^(-q_i/10)`, pattern IUPAC set `S_i` of size `m_i`:

```
S(o) = Σ_i  [ b_i ∈ S_i ]  · log2( 4(1-e_i) / m_i )
     + Σ_i  [ b_i ∉ S_i ]  · w_i · log2( 4 e_i / (3 (4 - m_i)) )
```

(numerator = P(base | adapter present), denominator = P(base | random sequence); `w_i` = 1.0 hard / 0.5 soft.)
Behaviour for `m=1`: a match gives **+2.00 bits**; a mismatch at q30 **−11.0 bits**, at q10 **−4.5 bits**, at q2 (an `N`) **−1.8 bits**. Low-quality mismatches are near-free, high-quality ones are lethal — MIGEC's intent, continuous and correctly indexed.

**Threshold** (E-value / Bonferroni over the search space):

```
S_min = log2( L_read · N_patterns / α ),   α = --fdr (default 0.01)
```

e.g. `L=150, N=96, α=0.01 → S_min ≈ 20.5 bits ≈ 10 clean matched bases`. Accept iff `S(o*) ≥ max(S_min, --min-score)` **and** `S(o*) − S(o₂) ≥ Δ` (`--min-score-gap`, default 5 bits) — otherwise the placement is ambiguous and the read goes to the `ambig` bucket. Cross-row ambiguity (two sample rows both match) uses the same Δ.

**Capture-level filters** (post-placement): reject UMI if `#N > --umi-max-n` (default 0), or homopolymeric (`max_b count_b ≥ 0.9·L`, on by default for 2-colour chemistries). **Do not** apply MIGEC's hard `min q < 15 → drop`; record `umi_minq` and let refine decide — dropping sequence information at checkout is exactly what the user asked to avoid.

**Indels.** Default `--indels 0`: Illumina indel rate ≈1e-6/base, no pattern needs it. With `--indels k` (k≤2) the engine switches to Myers + a banded (band=k) traceback that recovers *capture boundaries* — an indel 5' of a capture shifts it, so captured lengths become `L±k`; those are padded to the declared length and flagged `ZI:i:<indel>`. **v1 cut: allow indels only in the trailing adapter run used for trimming, not in the capture-bearing region.**

---

## 3. Dual barcoding, orientation, undef buckets

* **Declaration is per stream, not per role**: `pattern_r1 / pattern_r2 / pattern_i1 / pattern_i2`. Any may be empty.
* `--orientation fixed` (default with `-o/--oriented`): each pattern searched in its own stream only.
* `--orientation auto` (= MIGEC master/slave): the **primary** pattern (`--primary r1`, else first non-empty) is searched in R1, then in R2; if it hits R2 the mates are **swapped** so the primary is always emitted as mate 1. The secondary pattern is then searched in the other mate only.
* `--rc-mask 0:1` (MIGEC-compatible) RCs mate 2 after orientation. `--rc-search` additionally searches each pattern's reverse complement in each mate (amplicons sequenced both ways); strand recorded as `ZS:A:+/-`. Subsumes the old `--rc-barcodes`.
* **Concatenation**: all `UMI*` captures in (stream order) × (declaration order) → one UMI; likewise `CB*`, `SB*`; quality strings concatenate identically.
* **Duplex** (`--umi-canonical duplex`, exactly two UMI captures α,β): emit `UMI = min(α∥β, β∥α)` lexicographically, strand bit in `ZS`. The two single-strand families remain recoverable.

**Undef buckets** (`--undef {none,merged,split}`, counts always reported):

| bucket | meaning |
|---|---|
| `undef-primary` | no placement of the primary pattern above `S_min` |
| `undef-secondary` | primary placed, secondary not |
| `undef-sample` | patterns matched, SB/CB whitelist lookup failed |
| `undef-umi` | patterns matched, UMI failed N/homopolymer filter |
| `ambig` | ≥2 sample rows matched with ΔS < Δ (index hopping) |

---

## 4. Whitelists

**Format**: one barcode per line, optional TSV columns 2 = translated barcode (10x multiome `3M-february_2018.txt.gz` is 2-column), 3 = prior weight. `#` comments. `.txt`/`.tsv`/`.gz`. 10x `737K-*.txt.gz` load unchanged. Declared per row (`whitelist_cb`, `whitelist_sb`) or globally (`--whitelist CB=path`).

**Matching**:
1. Exact hit → `ZC:i:0`.
2. Else enumerate the 3·L single-substitution neighbours and probe a `flat_hash_set` of 2-bit-packed barcodes (48 probes for L=16 — faster than any index). For `--whitelist-max-edits ≥ 2`, build a `seqtree::Index` once and use `Searcher::search` with `SearchParams{max_substitutions=2}`.
3. **Posterior tie-break (Cell Ranger rule, made explicit)**: for each candidate `c` differing at position *i*,

```
P(c | obs) ∝ π_c · e_i / 3 ,      e_i = 10^(-q_i/10)
π_c = (n_c + 1) / (Σ_c' n_c' + |W|)      # Laplace-smoothed exact-hit abundance, pass 1
```
   assign `argmax_c P(c|obs)`; accept iff normalised posterior ≥ `--whitelist-min-post` (default **0.975**). Needs two passes (pass 1 = exact counts) — cheap because pass 2 reads the `.mgc` intermediate, not gzip FASTQ. `--whitelist-prior uniform` gives a single-pass mode (π ≡ 1, pure quality tie-break).

**Missing-barcode reporting**:
* `checkout.whitelist_<tag>.tsv` — `barcode, in_whitelist, exact_reads, corrected_reads, total_reads, rank`
* `checkout.whitelist_<tag>.missing.tsv` — whitelist entries with `total_reads < --min-reads-present` (default 1), plus a summary row `n_expected, n_observed, n_missing, fraction_missing`.
* When the sheet enumerates expected sample barcodes, the same report is keyed by `sample_id` → a dropped-out library is a single line.

---

## 5. Sample sheet and output manifest

`--sheet FILE` — TSV/CSV, delimiter autodetected, `#` comments, header required (except legacy):

| column | meaning |
|---|---|
| `sample_id` | **required**, unique, filesystem-safe |
| `pattern_r1`, `pattern_r2`, `pattern_i1`, `pattern_i2` | pattern per stream |
| `index1`, `index2` | literal i7 / i5; shorthand for `pattern_i1 = ^<seq>` |
| `preset` | `10x-3p-v3`, `10x-3p-v4`, `10x-5p-v2`, `10x-multiome-rna`, `migec-legacy`, `duplexseq`, `smartseq-umi` — fills patterns + whitelist + rc_mask |
| `whitelist_cb`, `whitelist_sb` | whitelist paths |
| `r1`, `r2`, `i1`, `i2` | input FASTQ paths (comma-separated for multiple lanes, or repeat the row) |
| `lane` | label, carried as `RG` |
| `rc_mask` | e.g. `0:1` |
| `orientation` | `fixed`\|`auto` |
| `expected_cells` | int, consumed by the barcode-rank/refine stage only |
| `notes` | ignored |

**Legacy autodetect**: no header line, 2–5 columns, column 2 matches `^[ACGTNacgtn]+$` → MIGEC `barcodes.txt` mapped to `sample_id, pattern_r1, pattern_r2, r1, r2` with `--pattern-dialect migec`, `--orientation auto`.
**Cut for v1**: native Illumina `SampleSheet.csv` v1/v2 parsing. Ship `migec sheet from-illumina` as a ~40-line Python helper — no C++.

**Output manifest** `checkout.manifest.tsv` (one row per output stream; this is the contract for sort/refine/assemble):

```
sample_id  file  format(mgc|fastq)  mates  umi_len  cb_len  has_cb  reads  reads_pf
n_undef  patterns_used  whitelist  qual_bins  checksum_xxh3
```
Plus `checkout.report.json` / `.tsv` per sample: `reads_total, reads_matched, reads_umi_ok, reads_wl_exact, reads_wl_corrected, reads_wl_failed, reads_ambig, mean_umi_qual, umi_effective_bits, est_error_rate`, and per-bucket undef counts.

---

## 6. Placement suggestion tool

```
migec suggest --r1 a.fq.gz [--r2 b.fq.gz] [--i1 …] [-n 200000] [--max-len 60]
              [--whitelist PATH ...] [--json suggest.json] [--plot DIR]
```

Streams the first *N* reads (default 2·10⁵), no correction. Steps:

1. **Per-cycle composition** `f_{p,b}`. Report entropy `H_p = -Σ_b f_{p,b} log2 f_{p,b}` ∈ [0,2], information `I_p = 2 − H_p − 3/(2 ln2 · N)` (Schneider small-sample correction), GC, N-fraction, mean Phred, frac q<20.
2. **Segmentation** by class: `constant` if `I_p ≥ 1.8`, `random` if `I_p ≤ 0.1`, else `structured`. Merge maximal same-class runs, drop runs <3 cycles, emit intervals. *(A proper HMM here is over-engineering — cut.)*
   * `random` run of length 8–14 immediately 5′ of a `constant` run → **UMI candidate**.
   * `structured` run at read start, few distinct values (≤ few hundred) → **sample-index candidate**; 10⁴–10⁶ distinct values → **cell-barcode candidate**.
3. **Conserved k-mer / adapter discovery**: count all k=12-mers with their start position. Candidate iff `max_p count(kmer,p) ≥ 0.2·N` (fixed) or `Σ_p count ≥ 0.2·N` with positional spread (drifting). Greedily extend the top 10 into a maximal consensus while per-base frequency ≥ 0.8 → recovers e.g. `GTGGTATCAACGCAGAG` verbatim.
4. **Whitelist probe** (highest-value output): for each bundled whitelist and each offset *o*, exact-hit fraction of `read[o:o+16]`. Hit fraction > 0.5 identifies the 10x chemistry unambiguously.
5. **Emit a paste-ready pattern** per stream with per-element confidence.

**Reported statistics** (`suggest.json` + TSVs for gnuplot):
`cycles.tsv`: `mate,pos,A,C,G,T,N,entropy_bits,info_bits,mean_q,frac_q20,class`
`segments.tsv`: `mate,start,end,class,mean_info,suggested_tag`
`kmers.tsv`: `kmer,count,modal_pos,positional_fraction,extended_consensus`
`whitelist_probe.tsv`: `whitelist,mate,offset,len,hit_fraction,distinct_observed`
`suggested_pattern_r1/r2/i1/i2`.

**Same code path yields the two estimators the requirements ask for:**

* **UMI quality (deviation from uniform)** over the UMI cycles:
  `D_UMI = (1/L) Σ_p (2 − H_p)` bits (0 = perfectly uniform), plus
  `X² = Σ_p Σ_b (n_{p,b} − N/4)² / (N/4)`, df = 3L, with p-value. PWM written to `umi.pwm.tsv` (successor of MIGEC `pwm.txt`).
* **Effective UMI space and birthday collisions**: `umi_effective_bits = Σ_p H_p`, `S = 2^(Σ_p H_p)`; for M observed molecules, `P(a molecule's UMI is shared) ≈ 1 − exp(−(M−1)/S)`, expected #multi-occupied UMIs `≈ S(1 − (1 + M/S)e^{−M/S})`. Reported so the user sees immediately that a biased UMI has a much smaller effective space than 4^L.
* **Empirical error rate / quality recalibration**: over `constant` segments, mismatch rate vs consensus stratified by reported Phred: `ê(q) = mismatches(q)/bases(q)`, written to `checkout.qual_calibration.tsv` next to the nominal `10^(-q/10)`. Feeds `ScoreModel` (§2) and the consensus quality cap in assemble.

---

## 7. Output

**Recommend both; binary is the default.**

* `--out-bin DIR` (**default**) — one `.mgc` shard per sample. The next stage sorts by (sample, CB, UMI) on disk; re-parsing gzip FASTQ for that is the dominant cost. `.mgc` is ~2.5–3× smaller than gzip FASTQ (2-bit bases, binned quality) and needs zero parsing.
* `--out-fastq DIR` — demultiplexed FASTQ(.gz), SAM-tagged headers, for handing to arda / minimap2 / bwa-meme.

### FASTQ header convention (default `--header-style sam`)

```
@<orig-id> RG:Z:<sample_id> BC:Z:<sample_bc> QT:Z:<sample_bc_qual> CR:Z:<raw_cb> CY:Z:<cb_qual> UR:Z:<raw_umi> UY:Z:<umi_qual> [CB:Z:<corrected_cb>] [ZS:A:+]
```

These are the standard SAM tags (`BC/QT` sample barcode+qual, `CR/CY/CB` raw/qual/corrected cell barcode, `UR/UY/UB` raw/qual/corrected UMI). **Why this exact format:** `bwa mem -C` and `minimap2 -y` copy everything after the first whitespace of the FASTQ comment verbatim into the SAM record, so already-valid `TAG:TYPE:VALUE` fields become valid SAM tags with no conversion step; `samtools fastq -T` round-trips them. Checkout writes only the *raw* tags (`CR/CY/UR/UY/BC/QT/RG`); `CB`/`UB` are added by refine (checkout emits `CB` only when a whitelist correction was applied).
Compatibility: `--header-style migec` → old ` UMI:<seq>:<qual>` suffix; `--header-style umitools` → `_<CB>_<UMI>` appended to the read *name* (what `umi_tools --extract-umi-method=read_id` and fgbio/fastp users expect).

### `.mgc` record layout (little-endian, zstd-framed blocks, streamable + seekable)

```
file header (64 B):
  "MGCR" | u16 version | u16 flags | u16 umi_len | u16 cb_len | u8 n_mates
  | u8 qual_bins | u16 sample_id_len | sample_id bytes  (padded to 64)

block: u32 comp_size | u32 raw_size | u32 n_records | u32 xxh3(raw)   [~1 MiB raw, zstd-1]
  record:
    u32 rec_len                       # skip/seek
    u64 key                           # 2-bit CB then 2-bit UMI, MSB-first  <- sort key
    u8  cb_codes[ceil(cb_len/4)]      # 2-bit packed (full codes; key is a prefix if >32 nt)
    u8  umi_codes[ceil(umi_len/4)]
    u8  cb_minq, u8 umi_minq
    u8  flags                         # b0 rc, b1 wl_corrected, b2 umi_has_N, b3 duplex_strand
    u16 name_len | u16 seq_len[n_mates]
    bytes name                        # only with --keep-names (default off, ~40% of bytes)
    u8  seq[mate]                     # 2-bit packed + N-mask bitmap when needed
    u8  qual[mate]                    # binned to qual_bins (default 8)
```

*Key width*: CB(16)+UMI(12)=28 nt=56 bits < 64 → the full sort key fits one `u64` for every 10x chemistry and every MIGEC design. `migec sort` therefore does an external LSD radix sort on `key` alone (blocks → runs → k-way merge) and never touches sequence bytes until the final write.
*Quality binning*: default 8 levels — NovaSeq/NextSeq RTA3 already emits only 4 distinct Q values, so this is lossless there and <1 Phred elsewhere. `--qual-bins 0` = full 8-bit, for benchmark parity runs.

---

## 8. C++ files / classes and the CLI

```
include/migec/
  types.hpp      Nt2/Nt4 codecs, std::array<int8_t,256> LUTs, phred tables, MigecError
  fastq.hpp      FastqRecord (views into a block buffer), FastqReader, FastqWriter
  reader.hpp     BlockReader: record-aligned blocks from plain/gzip/zstd, lockstep mates
  pattern.hpp    Tag, PatternElement, Pattern, PatternProgram, PatternDialect,
                 compile_pattern(), compile_legacy_migec()
  matcher.hpp    ScoreModel{fdr, soft_weight, qual_calibration}, Matcher (shift-or / Myers),
                 MatchResult{offset, score, score_gap, captures, ok}
  barcode.hpp    PackedBarcode, Whitelist, WhitelistCorrector, CorrectResult
  sheet.hpp      SampleRow, SampleSheet, Preset, load_sheet(), builtin_presets()
  checkout.hpp   CheckoutConfig, CheckoutEngine, CheckoutStats
  record.hpp     MgcRecord, MgcWriter, MgcReader, sort_key()
  profile.hpp    CycleProfile, KmerProfile, SuggestReport, suggest_patterns()
  report.hpp     counters -> json / tsv
src/*.cpp (one per header) + src/_bindings.cpp   ->  migec_core (STATIC) + _core (pybind11)
python/migec/  cli.py (typer), sheet.py (from-illumina), plots.py (gnuplot emit)
```

Key signatures:
```cpp
Pattern compile_pattern(std::string_view spec, PatternDialect = PatternDialect::Modern);
bool Matcher::match(std::string_view seq, std::string_view qual, MatchResult& out) const;
CheckoutStats CheckoutEngine::run();   // owns the pool, py::gil_scoped_release
```
Errors as `throw std::invalid_argument("compile_pattern: unbalanced '(' at offset 12")`.

### CLI

```
migec checkout
  --sheet FILE                       sample sheet (TSV/CSV; legacy barcodes.txt autodetected)
  --r1/--r2/--i1/--i2 FILE           inputs (override sheet; repeatable for lanes)
  --p1/--p2/--pi1/--pi2 PATTERN      single-sample mode, no sheet
  --preset NAME                      10x-3p-v3 | 10x-5p-v2 | migec-legacy | duplexseq | ...
  --pattern-dialect {modern,migec}
  -o/--oriented                      == --orientation fixed
  --orientation {fixed,auto}         auto == MIGEC master/slave with mate swap
  --primary {r1,r2}
  --rc-mask 0:1                      --rc-search
  --fdr 0.01  --min-score BITS  --min-score-gap 5  --soft-weight 0.5  --anchor-drift 0
  --indels 0                         (>0: adapter-run only in v1)
  --umi-max-n 0  --umi-min-qual 0  --umi-homopolymer/--no-umi-homopolymer
  --umi-order UMI1,UMI2  --cb-order ...  --umi-canonical {none,duplex}
  --whitelist CB=PATH --whitelist SB=PATH
  --whitelist-max-edits 1  --whitelist-min-post 0.975  --whitelist-prior {abundance,uniform}
  --min-reads-present 1
  -t/--trim  --max-trim-nts 10       trim consumed prefix/suffix
  --out-bin DIR                      (default)   --out-fastq DIR   --gzip
  --header-style {sam,migec,umitools}
  --keep-names / --no-keep-names     --qual-bins 8
  --undef {none,merged,split}
  -p/--threads 0  --first -1  --deterministic/--no-deterministic
  --report-json checkout.report.json

migec suggest   (see §6)
migec sheet from-illumina SampleSheet.csv -o sheet.tsv     # python helper
```

---

## 9. Threading and IO

* **Pipeline**: one reader thread per input file produces **record-aligned ~4 MiB blocks**, mates advanced in lockstep by record count (block boundaries chosen on R1 record boundaries), into a bounded SPMC queue. *W* workers each take one block-tuple, match, and fill a per-worker output buffer. A writer thread drains, with one mutex per output sample. No per-read synchronisation — the arda style (plain `std::thread` over disjoint ranges, GIL released). Blocks carry a sequence number; `--deterministic` (default on) reorders at the writer, costing one small buffer.
* **Decompression: libdeflate** for gzip (2–3× zlib inflate, already in homebrew, MIT, no nasm), **zstd** for `.mgc`. `find_package` fallback to zlib. isa-l/igzip is ~1.5× faster still but is x86-first and pulls in nasm — not worth it on an M-series-first project (**cut**). gzip members are not splittable, so inflate is one thread per *file*; matching is what scales.
* **Targets** (assert these in the CI benchmark tier):
  * positional 10x pattern, uncompressed input, 8 threads: **≥ 20 M read-pairs/min** (matching is a few ns/read; the limit is memory bandwidth, ~3 GB/s of FASTQ text);
  * seed-and-extend pattern: ~1 µs/read/thread → **≥ 8 M reads/min/thread**;
  * **gzip input is the real ceiling**: libdeflate ≈ 350 MB/s/stream on M3 → a 30 GB gzip R1+R2 pair in **4–6 min** end-to-end.

---

## 10. Cut from v1 (explicit)

| cut | why / replacement |
|---|---|
| indels inside the capture region | Illumina indel rate ~1e-6; only allow them in the trailing adapter run for trimming |
| seqtree for 1-mismatch whitelist lookup | 3L hash probes are faster; keep seqtree for `--whitelist-max-edits ≥ 2` and for refine/assemble |
| native Illumina `SampleSheet.csv` parsing in C++ | `migec sheet from-illumina`, ~40 lines of Python |
| HMM segmentation in `suggest` | run-merging on the entropy classes is sufficient |
| paired-read **overlap/merge** at checkout (MIGEC `--overlap*`) | belongs in assemble (where contigs are built anyway) or is delegated to fastp/bbmerge; four flags and a whole matcher removed |
| `-e` template-switch trimming (`^T{0,3}G{3,7}`) as a special flag | expressible as a pattern element `.{0,3}G{3,7}`; no special case |
| `--rc-barcodes` as its own flag | subsumed by `--rc-search` |
| `--max-mismatch-hq` hard cap | redundant with the LLR |
| writing FASTQ **and** binary by default | binary only unless `--out-fastq` |
| read names in `.mgc` by default | ~40% of record bytes; `--keep-names` when provenance is needed |
| MIGEC's hard "drop UMI if any base Q<15" | record `umi_minq`, decide in refine — keeping low-coverage/low-Q UMIs is an explicit requirement |