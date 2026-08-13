# 1. CUT LIST (ranked by complexity removed)

| # | Component | Design | Why cut / deferred | Covered instead by |
|---|---|---|---|---|
| 1 | **Two of the three intermediate binary formats** (`.mgc` §7 checkout, `.migb` §0 assemble) | checkout, assemble | Three teams specced three incompatible on-disk records for the same data. Only one can exist. | `.mig` (sort design §1) — it is the only one with framing, CRC, footer, provenance, and a written truncation story. |
| 2 | **External merge sort / run generation / k-way merge / spill files** (`extsort.hpp`, `ExternalSorter`, `--tmp-dir`) | repo, checkout | ~600 lines, two extra full passes over 20 GB, for an ordering property nothing consumes. | Hash-bucket partition at checkout write time + in-RAM pdqsort per bucket (sort §2). `nbuckets==1` *is* the in-memory case; no second code path. |
| 3 | **`migec suggest`'s HMM-free-but-still-large machinery: k-mer discovery, greedy consensus extension, whitelist probe, 4 TSV outputs, Schneider correction** | checkout §6 | User asked for "tools to SUGGEST the actual UMI/primer placement". Per-cycle entropy + run segmentation + a printed pattern satisfies that. The k-mer/adapter discovery and whitelist probe are the useful 20%. | Keep: per-cycle composition/entropy/mean-Q, run segmentation, top-10 conserved 12-mer extension, suggested pattern string. Cut: Schneider small-sample correction, `whitelist_probe.tsv` as a separate product (fold into one `suggest.json`), positional-spread scoring. |
| 4 | **Contig assembly (`--contigs`, `contig.hpp`, `extend_frame`, `--contig-iters`, `--contig-split-gaps`, `--trim-depth`)** | assemble §5 | User marked it OPTIONAL. The union-frame at `--max-offset>0` already gives it for free; the greedy re-extension is the only new code and it is untested against any truth set. | Union frame only. Ship `--max-offset`; document that overlapping read stacks extend the consensus. Defer `--contigs` to v2.1. |
| 5 | **Doublet detection's ambient guard + `hop_rate` 1/64 sampled estimator + `tags.mtx` + `--doublet-key cluster`** | assemble §4 | User asked for "detect cell doublets in scRNA-seq". `chain_multiplicity` (n_major_variants > 2) delivers that. The binomial ambient test needs a cross-cell frequency table and a Bonferroni sweep; the hop-rate estimator is a diagnostic nobody asked for; `tags.mtx` is CellPlex scope creep. | `cells.tsv` with `n_major_variants` + `doublet_call`. Defer ambient guard and hop-rate to v2.1. Bench design already cuts real-data doublet validation. |
| 6 | **The `migec stats` command as a separate stage** | refine | It is `refine --dry-run`: identical estimators, identical TSVs, no rewrite. A second command to document, test, and keep in sync. | `migec refine --no-apply` (writes all TSVs, no record rewrite). |
| 7 | **`migec correct` as a separate stage** (repo §2 CLI lists checkout/stats/correct/refine/assemble) | repo | Correction *is* refine. Five stages where three exist. | Fold into `refine`. |
| 8 | **`--indels` / `--umi-indels` / `--indel-mode banded` everywhere** | checkout, refine, assemble | Illumina indel rate ~1e-6/base. Three separate designs each carry an indel escape hatch nobody can benchmark in v1 (the only indel dataset, D14 IonTorrent, is itself cut). | Nothing. Ungapped/substitution-only throughout. Reinstate when ONT/IonTorrent lands. |
| 9 | **seqtree in the UMI 1-substitution path** | refine §2.1 (already self-cut), checkout §4 (already self-cut) | Both designs independently concluded 3L XOR-probe beats a trie for fixed-length ≤1-sub. Ratify it globally. | 2-bit hash + `3L` neighbour enumeration. seqtree used **only** for whitelist ≥2 edits — and even that is `--whitelist-max-edits 2`, off by default. |
| 10 | **seqtree C++ linkage (FetchContent/submodule) and the upstream install/export patch** | repo §3.1 | Build fragility for zero measured gain; one batched `search_batch` per stage. | Plain `seqtree>=0.6` runtime dep. Actually — with cut #9, v1 may not import seqtree at all; make it an **extra**, not core, and only if `--whitelist-max-edits 2` is used. |
| 11 | **`migec simulate` as a shipped user-facing command with 20 SimConfig fields, 7 truth files, MANIFEST.tsv** | repo §7.1 | This is the single most valuable test asset (keep it) but it does not need to be a public CLI command with chimera/doublet/lognormal-CV/PCR-efficiency knobs on day one. | `tests/synthetic/_sim.py`, ~150 lines, fields: n_molecules, umi_len, coverage, seq_error, umi_error, collision_rate, pcr_error, seed. Emits `truth_reads.tsv` + `truth_consensus.fasta`. Promote to `migec simulate` in v2.1 once the assertions stabilise. |
| 12 | **`migec subsample` as a C++-backed command with `--order random-umi`, `--per-sample`, exit assertion, `subsample.json` sha256s** | repo §8.1 | The *rule* (never subsample reads) is non-negotiable and must be in CLAUDE.md. The *command* is 40 lines of Python over the checkout output, used only to build notebook fixtures. | Keep as `migec subsample` but Python-only, two flags (`--umis N`, `-o`), hash-partition selection (bench §2's `blake2b mod` — unbiased, beats "first N lexicographic"). No `_core.subsample_by_umi`. |
| 13 | **`migec report` HTML** | repo | Already self-cut. Marimo notebooks cover it. | Notebooks. |
| 14 | **`migec view`** | sort §7 | Debug tool. | `migec refine --no-apply` + a pytest fixture decoder. Add `view` when someone actually needs it. |
| 15 | **Paired-read overlap/merge at checkout** | checkout §10 (self-cut) | Ratify. | fastp/bbmerge upstream; or nothing. |
| 16 | **uBAM/htslib writer, BGZF, `--codec` flag, lz4, libdeflate, isa-l, `--qual-bins`** | sort, repo | All self-cut in their designs. Ratify all. | zlib for gzip in, zstd-1 for `.mig`, plain multi-member gzip out. `samtools import -T '*'` for BAM. |
| 17 | **Illumina `SampleSheet.csv` parsing, `migec sheet from-illumina`** | checkout §5 | Nobody asked. 40 lines of Python is still 40 lines + a doc page + a test. | Write the TSV by hand. Defer. |
| 18 | **`--pattern-dialect migec` legacy grammar translator** | checkout §1 | **Note: Judgement call, flagged not cut.** It is ~60 lines and it is what lets G5/G7/G10 compare against MIGEC v1.2.9 on the same `barcodes.txt`. Keep. But cut the *autodetection* of headerless barcodes.txt — require `--pattern-dialect migec` explicitly. | Explicit flag. |
| 19 | **`--legacy-directional` (umi_tools rule reimplemented in our C++)** | refine | Benchmarking-only feature living in shipped code. umi_tools 1.1.5 is a pip install; run the real thing. | Run actual umi_tools in the bench repo. |
| 20 | **Chimera flag, `--tie-eps`, `--fast-dedup`, `--min-cons-qual`, `--ambient-inflation`, `--offset-margin` exposed as CLI** | assemble §8 | Six knobs, each one a doc line, a test, and a support question. | Hard-code the defaults as constants in `consensus.hpp` with the measurement that justifies them (arda style). Promote to flags when a benchmark demands it. |
| 21 | **D9, D10, D14, 40/42 runs of D6; Duplex-Seq-Pipeline; starcode/cd-hit/rainbow/du-novo; running Cell Ranger** | bench §8 | Already self-cut. Ratify. | D1,D3,D4,D5,D7,D12,D13 only. |
| 22 | **`isalgo/umi_data` `synthetic/calib_sim` tier (multi-GB LFS)** | bench §2 | Calib simulate is a one-line regenerate. Storing GBs of regenerable reads on HF is exactly what SOURCES.md exists to avoid. | `SOURCES.md` entry with the exact `calib simulate` command + seed. Store only the *truth* file if the command is non-deterministic. |
| 23 | **`--emit-cb-umi-graph`, `--emit-barcode-table`, `ambient.tsv`** | refine | Three optional 10M-row outputs, all default-off, none consumed by anything in v1. | Defer all three. |
| 24 | **`--whitelist-prior uniform` / two-pass abundance prior** | checkout §4 | The two-pass exact-count pass exists anyway (refine sees the sorted store); the *flag* to disable it does not. | Always abundance-prior. One flag fewer. |
| 25 | **`--orientation auto` mate-swap AND `--rc-search` AND `--rc-mask` AND `--primary`** | checkout §3 | Four orthogonal orientation controls. MIGEC had two. | `--rc-mask` (MIGEC-compatible) + `--orientation {fixed,auto}`. Cut `--rc-search` and `--primary` (primary = first non-empty pattern). |

**Not cut, explicitly, despite looking cuttable** — these are user-explicit or on the statistical path:
dual-end/duplex UMI concatenation; the quality-aware LLR acceptance score; the birthday-collision prior `A_ind`; the neighbour method-of-moments error estimate; the ZTNB mixture and its FDR threshold `m_hi`; the keep-low-coverage-orphan rule and its `gamma_mig` derate; the RT-error quality cap; multi-consensus-per-UMI splitting with `d_split` + ΔBIC; UMI PWM / KL / `L_eff`; whitelists + missing-barcode report; sample sheets; coverage plots; on-disk sort in C++.

---

# 2. CONFLICTS AND DUPLICATION — with decisions

### C1. Three intermediate record formats
- checkout: `.mgc`, magic `MGCR`, **2-bit packed seq**, binned quality, `u64 key` = CB‖UMI, LSD radix sort.
- sort: `.mig`, magic `MIGB`, **raw ASCII seq+qual + zstd-1**, `unsigned __int128` key, pdqsort, hash buckets.
- assemble: `.migb`, packed2bit, group headers with `f32 p_umi`.

**DECISION: `.mig` (sort design) wins, with two amendments.** Sort's arithmetic is the only one that measures the alternative (227 B/pair packed vs 197 B/pair raw+zstd-1) and it is the only format with CRC, footer, and provenance JSON. Amendments: (a) add `u8 umi_minq, cell_minq` — already present, good; (b) **add `f32 p_umi` to the record**, which assemble needs and sort's layout omits — but it is written by *refine*, not checkout, so it lands as a per-group side file `<sample>.<bucket>.pumi` (u32 group ordinal → f32) rather than bloating every record. Consequence: checkout writes `.mig`, refine writes `.pumi` + the correction map, assemble reads both. No re-write pass, no second format.

### C2. Grouping unit: buckets vs global sort order
- sort: hash-partition, **no global order**, `migec sort` not exposed.
- refine: assumes a `(sample_id, cell, umi)`-sorted store; applies the correction map "during the merge phase of the external sort".
- repo: `extsort.hpp` / `ExternalSorter` / `RecordKey{sample,cell,umi}`.

**DECISION: hash buckets, no global order, no `migec sort` command.** Refine's write-back moves from "merge phase" to "bucket rewrite": each bucket is read, corrected keys applied, re-sorted in RAM, rewritten. Refine's correction is per-scope (`cell` or `sample`), and a hash on `(cell,umi)` does **not** co-locate the 1-mm neighbourhood — so refine must build its `BarcodeTable` from a **streaming pass over all buckets collecting only (key,count,meanq)**, which is `n × (8+4+L)` bytes = ~400 MB at 10M distinct, well within budget. Correction decisions are global; only the *record rewrite* is per-bucket. Delete `extsort.hpp` from the repo tree.

### C3. Two header tag conventions
- checkout: `RG:Z BC:Z QT:Z CR:Z CY:Z UR:Z UY:Z CB:Z`, space-separated, name = original read id.
- sort: name `@<sample>.<mig>[.<g>]:<CELL>:<UMI>`, **TAB**-separated tags, `RX QX CB BC MI cD cM cE`.
- assemble: name `@<sample>:<cb>:<umi>:<ci>`, TAB tags, `RX QX CB CR MI cD cI cK cM cU cW cS`.

**DECISION: sort's spec, verbatim, is the *output* convention** — it is the only one backed by verified evidence (bwa `-C` requires SAM-conformant TAB-separated comment; dnaio drops comments so arda only sees the *name*; fgbio `CopyUmiFromReadName` splits on `:` and takes the last field). Two amendments: (a) assemble's `cI`/`cK` replace sort's `<g>`-in-name-only, since a downstream `cK==1` filter is genuinely useful — so name carries `[.<g>]` **and** tags carry `cI`/`cK`; (b) drop `cM cU cW cS` from the header (they are in `mig.tsv`), keep `cD` and `cE`. Checkout's *intermediate* FASTQ (`--out-fastq`, the demultiplexed-but-not-assembled path) uses the raw-tag set `RG BC QT CR CY UR UY` — different stage, different tags, no conflict, but **TAB-separated too**, not space.

Name format resolution: `@<sample>.<mig>[.<g>]:<CB>:<UMI>` — sort's, because UMI last is the fgbio/umi_tools contract and both designs agreed on it. Assemble's `<ci>` moves to the `.<g>` slot.

### C4. Overlapping CLI commands
Union across designs: `checkout, suggest, sheet, refine, stats, plot, correct, assemble, view, subsample, simulate, report, sort`. Thirteen.

**DECISION: five.** `checkout, suggest, refine, assemble, subsample`. See §3. (`plot` folds into `refine`/`assemble` as `--plots`; there is no need for a regenerate-from-TSV command when the TSVs are written in the same run.)

### C5. Where the MIG-size threshold is computed
- refine §5.1: ZTNB mixture → `m_hi` at FDR α, used as the *derate* boundary.
- assemble §8: `--min-reads INT(1; 0=auto from histogram)`.
- Old MIGEC: `Histogram` stage.

**DECISION: refine computes `m_hi` and writes it to `refine.json`; assemble reads it and never re-fits.** `assemble --min-reads` defaults to 1 (keep everything) and is a hard floor only. One model, one place.

### C6. Where the quality derate lives and how the floors combine
- refine §3.2: `p_final = 1 - (1 - p_cons)·gamma_mig`, `gamma_mig = (1-P(error))·P_real(c)`, Dirichlet-α consensus.
- assemble §3: `p_final = p_cons + p_rt + p_umi·δ`, LL-posterior consensus.

**DECISION: assemble's additive form, refine's inputs.** Additive in probability space is correct for independent small error sources and is the one that gets the RT floor right; refine's multiplicative `gamma` form conflates "this MIG might be spurious" with "this base might be wrong". Concretely:
```
p_cons(j) = 1 - P(b*|obs)             # assemble's LL posterior, exp-LUT
p_final(j) = p_cons(j) + p_rt + p_umi·δ
Q(j) = clamp(round(-10·log10 p_final(j)), 2, Qmax),  Qmax = -10·log10(p_rt)
```
`p_umi` = refine's `P(error|·)` (posterior the UMI is an unmerged error child), passed in `.pumi`. Refine's `P_real(c)` from the ZTNB mixture is **not** folded into base quality — it is reported as a per-MIG column and as the `cE` tag. Rationale: `P_real` is about whether the molecule is real, not about whether a base is right; multiplying it into every base Q double-counts and is not what a variant caller wants.

### C7. Two low-coverage-retention specs
refine §3.1 (5 rules, `tau`, `tau_lo`, `m_hi`, `m_lo`) vs assemble `--min-reads 1`.
**DECISION: refine's rules 1–4, with rule 5 (`m_lo` DROP) deleted** — it is unreachable at the default and its existence invites someone to set it. `tau=0.95`, `tau_lo=0.5`. Assemble drops nothing.

### C8. `--max-reads-per-umi` placement (sort's open question)
sort §3 caps at 100k in the IO layer; assemble §7 caps at 10k (`--max-reads-per-mig`, highest-Σqual-first).
**DECISION: one cap, in assemble, after sub-clustering, at 10 000, subsampled highest-Σqual-first.** The IO layer keeps everything (a bucket of 100k reads for one UMI is 20 MB — not a memory problem). This resolves sort's stated bias concern (option (a)) at zero cost. Delete `--max-reads-per-umi` from checkout.

### C9. `--umi-scope` default (refine's open question)
**DECISION: `cell` when CB present, `sample` otherwise; `WARN_SATURATED` fires and correction *continues*.** Firing the warning rather than suppressing correction is right — MIGEC's global gate had no statistical meaning (refine §1.5 shows it reduces to `0.0375·L`) and `A_ind(x) = n·π(x)` already suppresses merging per-barcode in the saturated regime. Surface the warning in `refine.json` and in the CLI epilogue. **Flag to user for confirmation** — this changes results vs MIGEC v1 on every bulk RepSeq sample.

### C10. seqtree's role
checkout says "do not add an index here", refine says "I would cut seqtree", repo says "runtime Python dep only", the project brief says seqtree is *the* fuzzy-search dependency.
**DECISION: seqtree is an optional extra (`pip install migec[whitelist2]`), used only by `--whitelist-max-edits 2`, which is off by default.** All three independent designs reached "don't use it" on their own hot paths; that is a strong signal. Ratify, and tell the user explicitly that v1 does not depend on seqtree.

### C11. Checkout writes FASTQ vs binary
checkout defaults `--out-bin`; sort's `checkout` writes `.mig` buckets; repo's pipeline implies FASTQ interop.
**DECISION: checkout always writes `.mig` buckets. `--out-fastq` is an additional, opt-in output for users who want demultiplexed-only FASTQ and will not run assemble.** Not a default, not both.

### C12. Determinism mechanism
checkout: `--deterministic` reorder buffer at the writer. sort: `src_index` tiebreak in the sort key. assemble: block-sequence reorder buffer.
**DECISION: `src_index` tiebreak (sort) + block-sequence reorder buffer (assemble).** Delete checkout's `--deterministic` flag — `.mig` block order is irrelevant because assemble sorts anyway, and `src_index` makes the result byte-identical for free. One mechanism, no flag.

### C13. License
repo §0 surfaces a real blocker: current LICENSE is MiLaboratory proprietary; seqtree/arda are GPL-3.0-or-later.
**DECISION: with C10 (seqtree not a v1 dependency), the GPL-linkage argument evaporates — the choice is free.** Still must be made by the user before the orphan-master push. Recommend **Apache-2.0 or MIT** for a library meant to be embedded in pipelines, with GPL-3.0-or-later as the alternative if matching arda matters more. **This is a question for the user, not a design decision.**

### C14. `migec assemble` input
sort: `migec assemble out/` (bucket dir). assemble: `migec assemble INPUT...` (`.migb` stream). refine: emits "a rewritten/merged record store".
**DECISION: `migec assemble REFINED_DIR/`** — a directory of `.mig` buckets plus `refine.json` plus `.pumi` side files. Positional dir, not a file list.

---

# 3. THE MINIMAL v1 SURFACE

**Five commands. Six file formats.**

1. `migec checkout --sheet S.tsv -o out/` — extract SAMPLE/CELL/UMI by position or pattern from R1/R2/I1/I2, LLR accept, whitelist exact+1mm correct, hash-partition to `out/<sample>.<bbb>.mig`; writes `checkout.manifest.tsv` + `checkout.json` (per-sample yields, undef buckets, whitelist missing-barcode table, qual calibration). `--out-fastq DIR` optional.
2. `migec suggest --r1 A.fq.gz [--r2 B.fq.gz] [-n 200000]` — per-cycle entropy/composition/mean-Q, run segmentation, conserved-kmer adapter recovery; prints a paste-ready pattern per stream; writes `suggest.json` + `cycles.tsv`.
3. `migec refine out/ -o ref/` — one streaming pass to build the barcode table; error-rate estimation (neighbour MoM + phred + whitelist), UMI PWM/KL/`L_eff`/saturation, ZTNB mixture → `m_hi`, 3L-neighbour correction with the birthday+phred+count posterior, keep-orphan retention; rewrites buckets with corrected keys, writes `.pumi`, all QC TSVs, and `--plots` PNGs. `--no-apply` = stats-only.
4. `migec assemble ref/ -o cons/` — group by (cell,umi) per bucket, modal draft + LL column accumulation, multi-consensus split (`d_split` + ΔBIC), RT/`p_umi`-floored quality, doublet call from chain multiplicity; writes consensus FASTQ + `<sample>.mig.tsv.zst` + `<sample>.cells.tsv` + `assemble.log.tsv` + `assemble.json`.
5. `migec subsample --umis 1000 -o example/` — hash-partition UMI selection, keeps **all** reads of selected UMIs (Python, for notebook/example fixtures).

**File formats:**

| Format | Written by | Purpose |
|---|---|---|
| sample sheet TSV (`sample_id, pattern_r1/r2/i1/i2, preset, whitelist_cb/sb, r1/r2/i1/i2, rc_mask, expected_cells`) | user | input; legacy `barcodes.txt` accepted via `--pattern-dialect migec` |
| `.mig` (magic `MIGB`, zstd-1 blocks, CRC32C, footer, provenance JSON; raw ASCII seq+qual; `u64 cell, u64 umi, u32 src_index`) | checkout, refine | the one intermediate |
| `.pumi` (u32 group ordinal → f32 `p_umi`) | refine | per-MIG spuriousness posterior |
| QC TSVs (`mig_size_hist`, `correction_map.tsv.gz`, `umi_pwm`, `qc_summary`, `barcode_rank`, `whitelist_missing`) | checkout, refine | the "barcode stats" + coverage plots deliverable |
| consensus FASTQ, `@<sample>.<mig>[.<g>]:<CB>:<UMI>` + TAB `RX QX CB BC MI cD cE cI cK` | assemble | the pipeline output, bwa `-C` / minimap2 `-y` / arda ready |
| `*.json` per stage (`checkout.json`, `refine.json`, `assemble.json`) | all | machine-readable summary; `refine.json` carries `m_hi` to assemble |

---

# 4. PHASED MILESTONE PLAN

Highest-risk-first ordering: the **consensus quality model** is the scientific claim and the hardest thing to get right, so it ships in M1 on synthetic data before any IO scaling work.

### M0 — repo skeleton + `.mig` round-trip + simulator (≈3–4 days)
Orphan-master migration per repo §1.2 (needs the license answer first). CMake/pyproject/CI(4 workflows, no Windows)/docs shell. `types.hpp`, `fastq.hpp` (kseq++ + zlib), `mig_record.hpp` + `mig_writer/reader`, `tests/synthetic/_sim.py`.
**Verify:** `ctest` green on ubuntu+macos; write 1M synthetic records → read back byte-identical for all lengths incl. 0; truncate the file at every block boundary → clean `MigecError`, no UB; `pip install -e .` then `import migec._core` in the CI python matrix.

### M1 — `assemble` on pre-grouped input; the consensus + quality model (≈1.5–2 weeks) ← **highest risk/value**
`consensus.hpp` (dedup, modal draft, offset placement, LL accumulation, exp-LUT posterior), `subcluster.hpp` (`d_split`, greedy clustering, ΔBIC), `quality.hpp` (additive floors per C6). Input: simulator FASTQ already grouped by UMI in the header. No sort, no refine, no checkout.
**Verify:** on `_sim.py` at `seq_error=1e-3, rt_error=1e-5, coverage∈{1,2,3,5,10,20,50}`: per-base error vs `truth_consensus.fasta` ≤ 1e-5 at coverage ≥5; **no emitted-Q bucket with n≥1000 anti-conservative** (`ê(Q) ≤ 2·10^-Q/10`), ECE ≤ 0.02 (bench G3); at `collision_rate=0.05`, split recall ≥0.9 / precision ≥0.95 vs `truth_umis.is_collision`; coverage-1 MIGs emitted at the read's own Q, floored at Q50.

### M2 — `checkout` (≈2 weeks)
`pattern.hpp` (grammar + `migec` dialect), `matcher.hpp` (shift-or + LLR + `S_min`), `sheet.hpp`, `whitelist.hpp` (exact + 3L probe + CR posterior), `.mig` bucket writer, undef buckets, manifest + report.
**Verify:** on `ci/migec_spikein` (D1 slice) with the original `barcodes.txt` and `--pattern-dialect migec`, per-sample read counts within 2% of MIGEC v1.2.9 `CheckoutBatch`; on D4 masks `NNNNNNNNNNNNtgact`/`agtcaNNNNNNNNNNNN`, dual-end 24nt UMI extracted, MIG count within 2% of MAGERI 1.1.1's published histogram (bench G9); on a `ci/tenx_hpbmc` slice with a **synthetic 100k whitelist**, exact+1mm correction rate matches a brute-force Python reference exactly; determinism: `-t 1` vs `-t 8` produce identical `.mig` content after sort (G8).

### M3 — `refine` (≈2 weeks)
`barcode_table.hpp`, `composition.hpp`, `error_model.hpp` (3 estimators), `mig_size_model.hpp` (ZTNB EM), `corrector.hpp` (3L hash + posterior + union-find), `cell_calling.hpp` (OrdMag + knee), `stats_report.hpp`, bucket rewrite + `.pumi`, `plots.py`.
**Verify:** on `_sim.py` with known `umi_error`, `err_rate_neighbour` within 20% of the injected rate across `error∈{1e-3,1e-2}`; UMI correction precision ≥0.99 / recall ≥0.95 vs `truth_reads.umi_true`; **G4: ≥95% of no-parent 3–5-read MIGs retained**, emitted at Q≤30, `e_out ≤1e-3`; `m_hi` from the FDR fit within ±1 of the empirical FDR-0.05 point; ZTNB EM converges in <100 iters on all four D-set histograms.

### M4 — end-to-end + `suggest` + `subsample` + notebooks + docs (≈1.5 weeks)
Wire checkout→refine→assemble, `profile.hpp` suggest, subsample, 5 marimo notebooks, all doc pages, `SOURCES.md`, `CHANGELOG.md`, `CLAUDE.md`.
**Verify:** `migec checkout|refine|assemble` on the D1 CI slice end-to-end, output FASTQ consumed by `bwa-meme mem -C` and `arda rnaseq run` without error and with `CB`/`RX` tags present in the BAM; `sphinx-build -W` green; `migec suggest` on the D3 primer data prints a pattern containing the verbatim `GTGGTATCAACGCAGAG`; `subsample --umis 1000` output has `mean_reads_per_umi` exactly equal to the input restricted to those UMIs.

### M5 — benchmarks, gates, HF dataset, release (≈2 weeks, partly cluster-bound)
`2026-migec-benchmark` repo, `isalgo/umi_data` mirror (`ci/` + `truth/` + `whitelists/README.md` only), Calib D12 ARI run, D5 8E5 consensus accuracy, D1 clonotypes vs MIGEC v1, D7 10x concordance, speed table, PyPI publish.
**Verify:** bench gates G1 (ARI ≥0.99, within 0.005 of Calib's published), G2 (`e_out ≤1e-5`, `S_f ≥100` on 8E5), G5 (3/3 D1 clonotypes, `FP ≤` MIGEC v1 at same minCount), G6 (`per_read_cb_concordance ≥0.995`), G7 (≥3× MIGEC v1 wall-clock at 8 threads, RSS ≤4 GB), G10 (≥99% of MIGEC v1's MIGs per D3 primer set). Tag `v2.0.0`, wheel smoke-test matrix green.

**Total ≈ 9–11 weeks.** M0–M1 are the go/no-go: if the quality calibration gate (G3) does not pass on synthetic data, nothing downstream is worth building.

---

# 5. THREE BIGGEST RISKS

### R1. The quality model is the product, and it is validated only against a simulator we also wrote.
`_sim.py` embeds the same error model (`seq_error`, `rt_error`, PCR) that `assemble` assumes. Passing G3 against it proves internal consistency, not calibration. A miscalibrated Q60 base propagates into false variant calls downstream — the costly error.
**Cheapest mitigation:** the 8E5 clonal control (D5 `SRR1763769`, 5.4 GB, single-end, 8nt UMI) is a *real* dataset where every non-polymorphic deviation is an error, i.e. ground truth with no simulator in the loop. Pull that one run in **M1**, not M5 — it is the smallest public dataset in the inventory and it makes G2/G3 real for the cost of one download. Two days moved earlier, whole-project risk removed.

### R2. Format/stage contract drift between three teams building three headers, three record layouts, and three CLI verb sets.
The designs already disagree six ways (§2 C1–C6). Left unresolved, M2 and M3 will each write half a format.
**Cheapest mitigation:** freeze `mig_record.hpp`, `key.hpp`, and the consensus FASTQ header string as a **one-page `docs/formats.rst` + a single doctest round-trip test written in M0**, before any stage code. Any later change must break that test. Cost: half a day in M0.

### R3. Scope regrowth — 13 CLI commands, ~90 flags, and 3 optional output tables were proposed for v1; each one is a doc page, a test, and a permanent support surface.
This is the single most likely cause of the project not shipping.
**Cheapest mitigation:** put the five-command surface (§3) and the cut list (§1) verbatim into `CLAUDE.md` §"Non-negotiables" with a rule: *a new CLI flag requires a failing benchmark that the default cannot pass.* Hard-code every constant in §1 item 20 with the measurement that justifies it (arda's C++ header-comment convention), so adding a flag requires deleting a measurement.

---

**Two questions that must go to the user before M0 executes:** (a) the LICENSE decision — the current file is MiLaboratory proprietary and cannot be inherited (with seqtree dropped as a v1 dep, MIT/Apache-2.0 is now available); (b) confirmation of C9 — `WARN_SATURATED` fires and correction *continues* on bulk RepSeq, which changes results vs MIGEC v1 on every such sample.