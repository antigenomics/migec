# Roadmap

**`checkout`, `refine` and `assemble` all work.** **Implemented:** the `.mig` intermediate format (reader, writer, range partitioning, CRC-checked
blocks, provenance and quality-calibration in the header), FASTQ IO for plain and gzipped input
with strict validation, barcode packing/unpacking and the IUPAC/Phred primitives, the pybind11
module, the read simulator with ground truth, CI (C++ on ubuntu+macos, Python 3.10/3.12 matrix,
pinned ruff), sphinx docs with a zero-warning gate, and the PyPI publish workflow.

**`migec checkout` works**: degenerate barcode patterns in the MIGEC dialect read verbatim from
published barcode tables, quality-aware log-likelihood acceptance, sample assignment with an
ambiguity verdict distinct from "unmatched", UMI extraction from non-contiguous runs, trimming of
adapter/tag/UMI, barcode transfer into SAM-style FASTQ headers, the power-of-two coverage
histogram, per-position base composition with Shannon entropy and information content, collision
entropy and effective barcode length, data-driven UMI error-rate estimation, and count correction
with a sequencing/polymerase/independent-molecule mixture plus the collision-corrected molecule
count. `scripts/spikein_ratio.py` computes the published spike-in validation metric, and
`scripts/compare_calib.py` scores UMI grouping against Calib by adjusted Rand index with splitting
and merging reported separately.

**Throughput and footprint**: 1.18 M reads/s at 16 threads (2 M single-end 115 nt reads, four
patterns), output byte-identical at any `-t`; ~22 bytes per distinct UMI against a hash map's ~48.
The counters are not yet partitioned, which is the open memory item and lands with `.mig` buckets.

Milestones are ordered by risk, not by pipeline order: the consensus quality model is the
scientific claim and is validated before any throughput work.

## M0 — skeleton, format, simulator

- [x] Archive the Groovy implementation on `legacy-v1` / `v1-final`, start master fresh
- [x] `.mig` format frozen in `docs/formats.rst`, with a round-trip and truncation test
- [x] FASTQ reader/writer; a record straddling a buffer refill is covered by a test
- [x] Read simulator with truth files, and tests for the simulator itself
- [x] CMake + scikit-build-core + CI + docs shell
- [x] Push to origin

## X — experiments that must precede the milestones they inform

- [x] **X1 — read-start dispersion within `(CB,UMI)` on one 10x run.** Done, 2026-08-13, on
      `pbmc_1k_v3` at MALAT1/ACTB/B2M — 8,326 multi-read groups. **The co-terminal assumption is
      false**: 7.8% of groups overall, and 0.3% of those with ≥6 reads, so what co-terminality
      there is is a small-group coincidence rather than a property of the chemistry. 92% of groups
      are wider than one 91 nt read. But **72.7% still form a single overlap component** (rising
      to 81.3% at ≥6 reads), so partitioning by overlap first reduces three quarters of the work
      to the ungapped problem MIGEC already solves — and the remaining **27.3%** is why the
      partition is mandatory: a single consensus over them asserts sequence across a gap no read
      covers. Also: only 1.5% of groups hold more than one read at this depth, so for shallow 3'
      GEX the UMI buys counting, not error correction. Written up in `docs/fragmented.rst`,
      script in `scripts/read_start_dispersion.py`.
- [x] **X2 — emitted-quality calibration on a clonal control, stratified by MIG size.** Done,
      2026-08-13, on `SRR1763769` (2.12 M reads, HIV-1 Primer ID). **The floor is of order 1e-4,
      not 1e-6**: 1.54e-4 [1.36e-4, 1.74e-4] at MIGs of ≥80 reads, so no emitted quality above
      **~Q38** is supportable, and the 1e-6 guess is excluded by two orders of magnitude. Matches
      the ~1 in 10,000 the source paper reports (doi:10.1128/JVI.00522-15). Note: The curve is still
      declining at 80 reads, so this is an upper bound: the 9 nt Primer ID puts the library at
      49.6% occupancy and `checkout` flags it `saturated`, so collided MIGs contribute mismatches
      counted here as error. Not the `p_floor + a/c` least-squares fit the plan specified — that
      model is wrong for a majority vote and returned a negative probability on simulated data.
      `docs/quality_floor.rst`, `scripts/quality_floor.py`.
- [x] **X3 — three permutation nulls on one deep real dataset.** Done, 2026-08-13, on
      `SRR1763769` — 125,369 distinct 9 nt barcodes at 47.8% occupancy. (1) **Position
      independence holds to ~1%** of the collision rate. Tested as a distribution — JSD against the
      product measure `q(u) = Π_j p_j(u_j)`, floored by a same-size draw from `q` — not by one
      functional. Invisible on all distinct barcodes (saturation), clear at reads ≥ 2 (z up to 33),
      and **entirely nearest-neighbour**: every adjacent pair positive, every distant pair zero.
      Cause measured: **0.55% of reads carry a barcode one base short**, a coupling that did not
      fire, which frameshifts everything after it. `Π_j m_j` stays; the 1.86x collision excess is
      the read threshold. Note: The first version reported 1.04x and all of it was artefact — N as a
      fifth base, and the plug-in `Σ p̂²` whose bias grows with k.
      (2) **97% of distance-1 pairs are chance** (844,243 observed, 817,358 under a column
      shuffle), and a size-preserving *count* shuffle over the fixed graph finds **~18,000 genuine
      error children**, plateauing from a count ratio of 5 upward. The permutation background puts
      the barcode error at 1.4e-3, 0.70x of the Phred + polymerase prediction, against
      `checkout`'s analytic 8.0e-4 at 0.39x — M3 takes the permuted background. (3) **The split threshold is 8.68, not 2.00**: a curveball
      randomisation preserving both margins of the reads x positions minor-allele matrix puts the
      1% false-positive point 19x above the nominal `p < 0.01`, because a low-quality read carries
      minor bases at many positions at once and mimics a linked subclone. Note: bootstrap 95% CI
      [8.42, 9.14]** over 82,800 randomisations; a tenth as many gave 9.61 and 11.66, so the
      interval is the number. (4) **Shallow libraries (1-3 reads/UMI) are a separate regime** and
      three of these results do not transfer to them — written up rather than quoted.
      `docs/nulls.rst`, `scripts/permutation_nulls.py`.

## M1 — `assemble`, the consensus and quality model

- [x] Grouping on the **whole** barcode — sample + cell + UMI. A UMI repeats across cells and
      samples by design; the sort key is `(cell, umi, src_index)` and the range partition is on
      the cell whenever there is one
- [x] Range partition into `.mig` buckets, one bucket resident: 531 k reads/s, 121 MB at 16
      buckets against 203 MB at one, output identical whatever the bucket count
- [x] Contig assembly (`--contig`): seed placement, union-find overlap components, one consensus
      per component, **never** bridged across a gap. This is one molecule's fragments only —
      full-length receptor assembly and doublet filtering are arda's job
- [x] Column log-likelihood posterior: `LL[j][b] = Σ_i (r==b ? log(1−e) : log(e/3))`
- [x] Sub-clustering by *linkage*, not by count of polymorphic sites. Threshold 8.68 (`-log10 p`,
      two-sided, Bonferroni'd within the MIG) from X3's false-positive curve — **not** the nominal
      2.00, which over-calls by 19x. Note: It implies a minimum group size: the strongest evidence a
      pair of columns can carry is `log10 C(n, n/2)`, so a 50/50 split needs ~34 reads to clear it
- [x] Quality floor **added**, not compared: `Q = −10 log10(p_cons + p_floor)`, default 1e-4 from
      X2, so nothing above ~Q38 is emitted
- [x] The birthday arithmetic re-run on the barcodes assemble saw: `expected_molecules_per_group`
      says how many molecules a group holds when the UMI is short by design, and contig mode warns
      when that makes contigs untrustworthy
- [x] Shallow libraries (1-3 reads/UMI) run, report the coverage histogram, threshold nothing, and
      say that the UMI is buying counting rather than error correction. Benchmarked as the
      memory-hostile shape: 190,595 reads/s at 1.02 reads/UMI, 259 B resident per distinct barcode
- [ ] `--rt-error auto`, fitted per dataset rather than taken from the default
- [ ] R1/R2 overlap merge (as a special case of placement, not a second matcher in checkout)
- Gate: per-base error ≤1e-5 at coverage ≥5 Done: (`tests/synthetic/test_assemble.py`, stratified by
  depth); `ê(Q) ≤ 2·10^(−Q/10)` for every bucket with n≥1000

## M2 — `checkout`

- [x] MIGEC-dialect pattern grammar, read verbatim from published barcode tables
- [x] Quality-aware log-likelihood acceptance; ambiguous distinguished from unmatched
- [x] UMI extraction from non-contiguous runs, trimming, SAM-style header transfer
- [x] Coverage histogram, composition/entropy/information, collision entropy, count correction
- [ ] Bit-parallel matcher (the current scan is O(offsets x pattern) and is not the bottleneck yet)
- [x] **`suggest`** — per-cycle PWM segmented into UMI / constant / payload, paste-ready pattern.
      Done early (2026-08-13) because X2 needed it; recovered SRR1763769's layout unaided.
- [x] **Barcode space and error budget built in** — nominal vs effective space, occupancy, Poisson
      λ, `p_multi`, and the Phred + polymerase error prediction against the distance-1 estimate.
      Logged, warned on, documented, notebooked, tested.
- [x] **Whitelists with a background hypothesis in the posterior; `N` expanded, not discarded.**
      `--cell-whitelist` in `refine`. The background prior is measured from barcodes at distance
      ≥2 from every entry, and it is a prior on *this barcode* -- the off-list read share divided
      by the distinct off-list barcodes -- because the whitelist prior is spread over every entry.
      Snapping scales as `n_parent · e/3`, so it needs a well-used parent *and* a poor base
- [x] Paired-end input; strand normalisation (tag searched in either mate, pair swapped)
- [x] Multi-core, byte-identical output at any thread count; compression on the workers
- [x] Speed and memory reported per run; `tests/benchmark/` regressions
- [x] **Dual-end barcodes** — column 3 of the sheet is MIGEC's slave pattern, on the other mate,
      extending the UMI rather than starting a new one. MAGERI's `NNNNNNNNNNNNtgact` /
      `agtcaNNNNNNNNNNNN` gives a 24 nt UMI; 2000/2000 assigned. Both halves or nothing
- [x] **`--max-offset`**, because a positional chemistry cannot be checked out without it. Also
      fixed: the acceptance bar was charged for the offsets a read *could* hold rather than the
      ones actually scanned, so an anchored 5 nt handle (10 bits) was billed 12.6 bits and refused
- [ ] `.mig` bucket output, which is also what bounds the UMI counters
- [ ] i7×i5 contingency table — the only way index hopping is actually estimable
- Gate: per-sample counts within 2% of MIGEC v1.2.9 on the spike-ins; identical output at 1 and 8
  threads ; >1 M reads/s at 16 threads

## M3 — `refine` (the stage works; cell calling and the FDR threshold are open)

**Measured before building** (`scripts/correction_accuracy.py`, 2026-08-13): the existing
count-ratio correction works from **3.1 reads/UMI upward** (recall ≥0.8, precision ≥0.9) and
collapses below it — 0.02 recall / 0.25 precision at 1.1 reads/UMI. `ε` follows: 0.98x of injected
at 3.1 reads/UMI, 0.20x at 1.1. Molecules are never lost (≥0.99 kept at every depth), so the
failure is missed correction, not destroyed data. That is the whole shallow regime, so M3's error
model has to use the evidence that survives at one read:

| reads/UMI | recall | precision | ε / true |
|---|---|---|---|
| 1.11 | 0.022 | 0.254 | 0.20 |
| 1.51 | 0.235 | 0.567 | 0.58 |
| 2.32 | 0.587 | 0.846 | 0.86 |
| 3.12 | **0.800** | **0.936** | 0.98 |
| 7.12 | 0.959 | 0.977 | 0.96 |

- [x] **Barcode base quality as evidence** (`BarcodeEvidence::position_error`). `QX` was carried
      through checkout and unread. A sequencing miscall in the barcode has *low* Phred at the base
      it changed; an early-PCR child has a high one in every read. Works at one read.
- [x] **Payload agreement as evidence** (`BarcodeEvidence::payload`). A barcode error child is a
      read of the parent's molecule, so its payload matches. Worth `log(1/clonality)`, and the
      clonality is **measured from the data** by sampling random barcode pairs — decisive in a
      diverse repertoire, worth nothing in a clonal library, and the number says which this is. It
      also lifts the count gates, which is what makes a singleton-vs-singleton merge possible at
      all, and it *refuses* merges the count ratio would have made on a disagreeing payload.
- [x] **The error likelihood is a rate, not a conditional.** The zero-truncated Poisson divided out
      `(1 − e^−λ)` — precisely the term saying whether a child should exist — so for a singleton
      child it tended to 1 for every λ and the error rate stopped mattering at exactly the coverage
      where nothing else was available. Untruncated, both sides are expected counts of
      neighbouring barcodes and the comparison is like for like.

  Measured after (`scripts/correction_accuracy.py`), against the achievable ceiling rather than
  against all children — a child whose parent barcode was never sequenced has nothing to merge into:

  | reads/UMI | reachable | recall of those | precision | molecules kept |
  |---|---|---|---|---|
  | 1.11 | 0.204 | 0.108 | 0.818 | **1.000** |
  | 2.32 | 0.904 | 0.816 | 0.830 | 0.987 |
  | 3.12 | 0.975 | 0.914 | 0.926 | 0.991 |
  | 7.12 | 1.000 | 0.979 | 0.997 | 0.999 |
  | 13.30 | 1.000 | 0.983 | 0.999 | 1.000 |

  Note: At ~1 read/UMI **80% of barcode errors are unfixable in principle** — the parent was never
  sequenced — and of the rest migec fixes 11% while destroying no real molecule. Precision is the
  side to err on: a wrong merge deletes a molecule and nothing downstream can tell, a missed
  correction only inflates the count.
- [x] **`migec refine` works**: barcode table, correction, read rewrite with `OX:Z:` preserving
      the original, `<sample>.barcodes.tsv`, coverage histogram after correction. Recovered 20,055
      molecules from 20,000 simulated with ε at 0.96x of injected. Holds the table, never the
      reads; three streaming passes
- [x] **Quality calibration measured against the pattern's own constant bases**, fitted as
      `ê(q) = ε_qi + a·10^(−q/10)` weighted by bases per Q. On `SRR1763769` the slope is 1.04 over
      46.3 M bases, so the reported Phred's *scaling* is right. Never: The intercept (3.9e-3) is **not**
      a sequencing floor -- the standard is a synthesised oligo, and synthesis runs ~1 defect per
      200-500 bases. It is spread evenly over all 23 anchor positions, none polymorphic, and
      matches the independently measured 0.55% one-base-short rate. Reported as a diagnostic of the
      primer and left out of `error()`
- [ ] Three error-rate estimators, sequencing vs quality-independent separation *of the template*
      (the pattern bases can only calibrate the primer, so the polymerase/RT split needs a
      different standard)
- [ ] Correction posterior: birthday prior with Rényi-2 collision entropy, phred, and a
      polymerase mixture component for early-cycle PCR children. The distance-1 background comes
      from X3's column shuffle, not from `C(n,2)·P_coll·shell`
- [x] **MIG-size threshold at a target FDR; keep-orphan retention.** The residual is measured, not
      derived: a surviving barcode that still looks like a child of a surviving neighbour, by count
      *or* by its reads agreeing on the molecule. Note: Count alone reports **zero residual at 1-3
      reads per UMI**, which is where it is worst. On a 1.23 reads/barcode library: 5.25% of 1-read
      molecules, threshold ≥2; on 4.62 reads/barcode: 0%, threshold 1. Never: Reported, never applied —
      every molecule stays in the output
- [ ] **Note: Bucketed correction.** A range partition on the top b bits splits a barcode from its
      neighbour for the top b/2 positions, so correction cannot be bucketed naively. The fix is two
      passes with the key rotated, so every pair shares a bucket in at least one. Until then the
      table is held whole and its size is reported
- [x] **Cell calling (OrdMag + knee)** and the QC tables. Molecules per cell, never reads; the
      knee reported next to the call and a warning when they disagree by more than 3x. On 500 real
      cells over 20,000 ambient barcodes it calls exactly the 500. `<sample>.cells.tsv`
- [x] **QC tables**: `<sample>.rank.tsv` (barcode-rank curve + CDF, log-spaced),
      `<sample>.bins.tsv` (per MIG size: barcodes, reads, fraction merged as error, payload
      entropy), drawn by `notebooks/refine_diagnostics.py`
- Gate: estimated ε within 20% of injected **at 1–3 reads/UMI, not only at 7**; ≥95% of no-parent
  3–5-read MIGs retained (already ≥99% at every depth measured)

## M4 — end to end

- [ ] `suggest`, `sort`, `subsample`, marimo notebooks, full docs
- Gate: output consumed by `bwa-meme mem -C` and `arda rnaseq run` with tags intact

## M5 — benchmarks and release

- [ ] `2026-migec-benchmark` repo, `isalgo/umi_data`, comparisons against MIGEC v1, MAGERI,
      UMI-tools, Calib, fgbio, Cell Ranger, **UMI-VarCal and UMIErrorCorrect** (the two UMI-aware
      callers benchmarked by Maruzani et al. 2024 for low-frequency ctDNA)
- [ ] ctDNA ground truth in the style of Maruzani et al.: COSMIC variants spiked at 0.005–0.075 VAF
      into a real cfDNA background at 200x/450x/850x. Never: Their deposited runs carry **no** UMI
      (aligned BAM submissions, `suggest` finds no pattern, `CMP_LINKAGE_GROUP` empty), so the UMIs
      have to be simulated — and their rule assigns them to reads sharing start *and* end, which is
      the co-terminal assumption X1 falsified
- Gate: grouping ARI ≥0.99; residual error ≤1e-5 on a clonal control; ≥3× MIGEC v1 wall-clock

## Deliberately not doing

- **Alignment and variant calling** — MAGERI's job; the pipeline hands off to arda/minimap2/bwa-meme.
- **Indels** — Illumina rates ~1e-6/base and no dataset to verify against. Substitutions only.
- **Duplex consensus (DCS)** — v2.0 extracts duplex tags and emits single-strand consensuses. No
  error-suppression claim in this repo is based on duplex data until DCS exists.
- **Full-length receptor assembly, doublet calling and contaminating-chain filtering** — arda's
  job. `--contig` assembles *one molecule's* fragments into a contig and stops there.
- **EmptyDrops-style cell rescue** — Cell Ranger's job. Ours is OrdMag plus a knee, and the
  benchmark gate is written against that rather than against a Jaccard we cannot reach.
- **An external merge sort** — range partitioning plus an in-RAM sort per bucket covers it, and
  `nbuckets == 1` is the in-memory case, so there is one code path rather than two.
