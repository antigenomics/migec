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

**Nothing in the pipeline scales with the library any more.** `checkout`'s counters and `refine`'s
barcode table both range-partition themselves past a byte budget, and correction follows both into
the partition in two passes with the key rotated.

**Throughput and footprint**: every stage threads and every stage is byte-identical at any `-t` --
checkout 1.55 M reads/s end to end (1.70 M matching), refine 1.55 M, assemble 2.47 M, all at 16
threads; ~22 bytes per distinct UMI against a hash map's ~48. Verified under the thread sanitizer
as well as by comparison at 1..16 threads.
The UMI counters range-partition themselves past a 1 GB budget, and correction follows them into
the partition in two passes with the key rotated, so nothing in `checkout` scales with the library
any more. `checkout --mig` writes the reads into that same partition, which is the pass `assemble`
otherwise builds for itself: 1.16 s -> 0.98 s end to end on 500 k reads for identical molecules.

**QC**: twenty gnuplot panels over the tables the stages write, including the four figures a user
already knows how to read -- Cell Ranger's barcode rank plot on unique UMIs, the MIG size spectrum
with molecules *and* reads, the rank/Zipf curve, and consensus quality as a box over an exact
`(depth, quality)` grid.

Milestones are ordered by risk, not by pipeline order: the consensus quality model is the
scientific claim and is validated before any throughput work.

## What is left, in the order it should be done

Everything below is either open or half-open. Nothing else on this page is.

| # | item | milestone | why it is next |
|---|---|---|---|
| 1 | The RT-vs-first-cycle-PCR split *inside* the floor | M3 | the sequencing/pre-amplification split now has a standard — the deep-MIG consensus residual, which is what `--pre-amp-error auto` fits (`docs/quality_floor.rst`). What is still unsplit is what MADE the floor, and that needs two chemistries to compare rather than a better estimator: the same template through an RT protocol and through a DNA one. Data, not code |
| 2 | The rest of the published comparisons | M5 | **MIGEC v1, UMI-tools and fgbio are done and scored** (`docs/validation.rst`, `docs/grouping.rst`). Open: MAGERI, Cell Ranger, UMI-VarCal, UMIErrorCorrect. This is what the version number is waiting on, not the code |
| 3 | The other three callers on the ctDNA arms | M5 | the ground truth is **found and scored**: `PRJNA788522` / `PRJNA507366`, certified VAF, real 12 nt inline UMI, and LoFreq is run end to end against it (reliable to 0.25%). What is open is Mutect2, UMI-VarCal and UMIErrorCorrect on the *same* consensus, so the comparison isolates the caller. Also open: re-running with adapter trimming, which is diagnosed but not measured |
| 4 | Bit-parallel matcher | M2 | last, deliberately: the scan is O(offsets x pattern) and is **not** the bottleneck. It goes in when a benchmark says so |

Nothing is blocked on data any more. **Britanova et al ageing** (bulk TCR, shallow -- the real
1-3 reads/UMI dataset) has been run: one HiSeq lane, ten donors, 149,588,907 read pairs at a 16 nt
UMI and 2.44-2.57 reads per UMI, in `assets/shallow_repertoire.tsv` and `docs/validation.rst`. The
ctDNA ground truth was **found** rather than built, by screening SRA read structure instead of
trusting the published claim that none exists (`scripts/sra_fetch.py probe`).

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
      **Q40** is supportable by default, and a blanket 1e-6 is excluded for an RT protocol by two
      orders of magnitude. Matches
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
      X2 and from 10x's own figure, so nothing above Q40 is emitted unless `--rt-error` names a
      higher-fidelity chemistry (`medium` 1e-5, `high` 1e-6)
- [x] The birthday arithmetic re-run on the barcodes assemble saw: `expected_molecules_per_group`
      says how many molecules a group holds when the UMI is short by design, and contig mode warns
      when that makes contigs untrustworthy
- [x] Shallow libraries (1-3 reads/UMI) run, report the coverage histogram, threshold nothing, and
      say that the UMI is buying counting rather than error correction. Benchmarked as the
      memory-hostile shape: 1,179,549 reads/s at 1.02 reads/UMI, 282 B resident per distinct barcode
- [x] **`--rt-error` names the chemistry** rather than guessing: `rt` 1e-4 (default, Q40 -- 10x's
      figure for the V(D)J RT and X2's own measurement), `medium` 1e-5, `high` 1e-6, or the rate.
      Never: it is the ONE-MOLECULE floor; 10x's Q60 needs two UMIs to agree and that is arda's job
- [x] **`--fast`, counting mode**: the modal exact sequence per group with the per-base best
      quality of the reads carrying it. No column model, so no error correction -- for when the
      deliverable is a molecule count. Refused with `--contig`
- [x] **Coverage capped at 10,000 reads per barcode into the consensus** (10x's rule). Never: the cap
      is on the reads consensed, never on the reads counted
- [x] **`--pre-amp-error auto`** (2026-08-14), fitted per dataset rather than taken from a named
      class: X2's estimator run on migec's own consensuses. Molecules at >= 20 reads, the library's
      modal sequence with one vote per molecule, real variation and divergent templates excluded,
      and the residual is the floor. Injected 1e-4 -> 1.20e-4 [7.94e-5, 1.75e-4]; injected 1e-5 ->
      1.58e-5 [8.14e-6, 2.76e-5]; injected 0 -> the bound 9.56e-6. Never: it REFUSES rather than
      guessing -- a diverse library (every position polymorphic) and a shallow one (no molecule
      deep enough) both fall back to the named class and say why, in the report and in
      `assemble.pre_amp_error.tsv`. Costs a second assembly pass, which is why it is opt-in. Note:
      the flag was `--rt-error` and that name is kept as an alias, but only an RNA library has a
      reverse transcription step -- on a DNA library the same floor is library-prep damage plus the
      first PCR cycle. `auto` fits the floor, not what made it
- [x] **R1/R2 overlap merge** (2026-08-14), as a special case of placement and not a second matcher
      in checkout: `assemble --mate2 <R2>` on the FASTQ route, `--merge-mates` on the `.mig` route
      where checkout already stored both mates in the record. Mate 2 is reverse-complemented and
      placed; overlapping mates give one consensus spanning the insert, non-overlapping mates give
      two contigs. Never: the offset is a property of the MOLECULE, so it is voted once per group
      over up to eight pairs rather than placed all-against-all -- the first version cost 11x the
      single-end path (119,820 record-pairs/s against 1,356,819), the vote costs 1,288,686 against
      2,115,912. Never: the mates are matched by POSITION and a file that ends early is refused,
      because pairing off what is left attaches one molecule's mate to another's
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
- [x] **Positional layouts as the primary mode** (2.0.0a2) — `^NNNN` / `^NNNNXNNN`, or half-open
      slices `0:8` / `0:4,5:10` / `cell:0:16,16:26`. `--max-offset` is now automatic and should not
      be passed: a caret, a slice list, a read structure and a pattern with nothing to score all
      anchor at 0. Presets for the eight chemistries with names, each carrying a citable source
- [x] **A run that matches nothing reports the declaration error** instead of three statistics
      computed from reads that never arrived
- [x] **The counter can bound itself** (2026-08-14): `UmiCounts::enable_spill()` range-partitions
      to disk past a byte budget and `for_each()` streams one bucket at a time, reducing on read.
      A spilled counter gives byte-identical histograms, compositions and distinct counts to a
      resident one; `entries()`, `find()` and `merge()` throw rather than answer from the fragment
      that happens to still be in RAM
- [x] **Bucketed correction, and the spill switched on in `checkout`** (2026-08-14). Correction
      follows the counters into the partition: pass 1 over the buckets as they stand owns the
      barcode positions the prefix does not touch, pass 2 over a rotated copy owns exactly the ones
      it hides, and every pair is weighed in one pass and only one. A bucketed run answers with the
      scalars -- `root`/`corrected` are indexed against `entries()`, which is the array being
      bounded -- and every one of them matches the resident answer field for field, on a simulated
      library with injected barcode errors and on a 500,000-read corpus. `umi_budget_bytes`
      defaults to 1 GB for the run; it costs ~2.2x the wall clock when it fires, and nothing when
      it does not
- [x] **`.mig` bucket output** (2026-08-14), `checkout --mig`: the reads are written into the
      same range partition, on the same key, that `assemble` builds in its first pass, so
      `assemble` reads them and skips that pass. One writer per (sample, bucket), owned by one
      thread for the whole run, so `-t` still changes nothing but the clock; the open-file budget
      is for the RUN, so a 96-plex sheet gets a couple of buckets each rather than 96 x 256.
      Opt-in, and FASTQ stays the default: a `.mig` file is an intermediate nothing else reads.
      500 k reads over four samples at `-t 4`: **1.16 s -> 0.98 s** end to end for the identical
      124,878 molecules, and the consensus FASTQ is byte-identical after decompression
- [x] **i7xi5 contingency table** (2026-08-14) — the only way index hopping is actually estimable
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
- [x] **Sequencing vs quality-independent separation *of the template*** (2026-08-14). The pattern
      bases can only calibrate the primer, so the standard is the template itself: at 20+ reads a
      consensus has suppressed sequencing error to nothing, and what it still gets wrong against
      the library's modal sequence is what was in the molecule before amplification.
      `--pre-amp-error auto` (M1) is that measurement. Still open, and it is data rather than an
      estimator: splitting the floor into RT and first-cycle PCR needs the same template through
      both chemistries
- [ ] Correction posterior: birthday prior with Rényi-2 collision entropy, phred, and a
      polymerase mixture component for early-cycle PCR children. The distance-1 background comes
      from X3's column shuffle, not from `C(n,2)·P_coll·shell`
- [x] **MIG-size threshold at a target FDR; keep-orphan retention.** The residual is measured, not
      derived: a surviving barcode that still looks like a child of a surviving neighbour, by count
      *or* by its reads agreeing on the molecule. Note: Count alone reports **zero residual at 1-3
      reads per UMI**, which is where it is worst. On a 1.23 reads/barcode library: 5.25% of 1-read
      molecules, threshold ≥2; on 4.62 reads/barcode: 0%, threshold 1. Never: Reported, never applied —
      every molecule stays in the output
- [x] **Bucketed correction** (2026-08-14). refine's table range-partitions itself past a byte
      budget and correction follows it, in two passes with the key rotated -- a plain partition
      splits a barcode from its neighbour for the top b/2 positions and would bound the memory
      while silently not correcting them. Two things it needed that the counter version did not.
      The table carries the **evidence**: a `BarcodeEvidence` indexed against `entries()` cannot
      survive a partition, and dropping it leaves the bucketed run on the count ratio alone, which
      reports nothing at 1-3 reads/UMI. And the two passes **scan** rather than merge, with one
      global apply afterwards: a barcode can have a plausible parent on each side of the boundary,
      and merging inside a pass takes the first rather than the best -- 2 barcodes in 6,591 landed
      elsewhere than the resident run put them. With both, a partitioned run and a resident one
      agree on every scalar and on every output file byte for byte
      (`tests/synthetic/test_refine_bucketed.py`). It also folded the evidence pass into the table
      pass, so refine now streams the reads twice rather than three times
- [x] **Cell calling (OrdMag + knee)** and the QC tables. Molecules per cell, never reads; the
      knee reported next to the call and a warning when they disagree by more than 3x. On 500 real
      cells over 20,000 ambient barcodes it calls exactly the 500. `<sample>.cells.tsv`
- [x] **QC tables**: `<sample>.rank.tsv` (molecule rank + CDF, log-spaced), `<sample>.bins.tsv`
      (per MIG size: barcodes, reads, fraction merged as error, payload entropy), drawn by
      `notebooks/refine_diagnostics.py`
- [x] **`<sample>.cell_rank.tsv`** (2.1.0) — Cell Ranger's barcode rank plot, cells sorted by
      **distinct UMIs** with the call on the curve. Never reads: one over-amplified molecule would
      put an empty droplet high up, which is the artefact the plot exists to show
- [x] **`<sample>.sizes.tsv`** (2.1.0) — the MIG size spectrum at **exact** sizes, one row per
      distinct depth. Power-of-two bins turn the rank/Zipf curve into four steps
- Gate: estimated ε within 20% of injected **at 1–3 reads/UMI, not only at 7**; ≥95% of no-parent
  3–5-read MIGs retained (already ≥99% at every depth measured)

## M4 — end to end

- [x] `suggest`, `subsample`, marimo notebooks, full docs
- [x] **Gate met, measured** (`docs/downstream.rst`): `minimap2 -ax sr -y` and `bwa mem -C` carry
      `RX`/`CB`/`MI` into a valid sorted BAM on 600/600 records; `arda amplicon` reads the
      consensus directly and its AIRR `sequence_id` **is** the molecule id; `salmon` and `kallisto`
      quantify it plainly. STAR unverified — the brew arm64 build reads 0 reads from any FASTQ
- [x] `plot` — twenty QC panels drawn with gnuplot from the tables the stages already write, and
      `assets/` holds the pipeline figure (graphviz) and the example panels the README shows,
      regenerated by `scripts/example_figures.py`
- [x] **The four familiar figures** (2.1.0) — barcode rank on Cell Ranger's axes, the MIG size
      spectrum on log1p with molecules *and* reads, the rank/Zipf curve, unique UMIs per sample
      barcode. Consensus quality is a **box** over an exact `(depth, quality)` grid, never a
      thinned scatter: quality is discrete and capped, so a cloud draws a flat line whatever the
      bin holds and the thinning removes the only thing it could have added
- [x] **Publication defaults** (2.1.0) — transparent background, one ink colour that reads on
      light and dark, key inside the plot box rather than in a gutter that widens every figure
- [x] **Overrepresented k-mers in `suggest`** — exact counts in a flat 4^8 array, measured against
      the reads' own composition, stitched back into the sequence they came from. Run on a stage's
      OUTPUT it answers "did the trim remove the primer", which nothing else here could
- [x] **`checkout.trimming.tsv`** — the payload length distribution after trimming. A pattern
      matched one base off still matches; this is where that shows
- Never: `sort` is **not** a command. Partitioning happens inside `assemble` and exposing it would
  document a third format with no independent meaning (`project/design-io-interop.md`).

## M5 — benchmarks and release

- [x] **UMI-tools and fgbio, scored** (2026-08-14). `scripts/compare_grouping.py` runs both
      map-first tools end to end -- barcode into the read name / into `RX`, minimap2 onto the
      simulator's own `clones.fa`, then `umi_tools group` / `fgbio GroupReadsByUmi -s adjacency` --
      and scores all three partitions against the simulator's truth with the same adjusted Rand
      index Calib is scored with. The result is what the position is worth: **on ONE reference
      migec wins** (ARI 0.9967 against 0.9864 and 0.9817) and wins on the direction that cannot be
      undone, putting 0.65% of reads into mixed clusters against 3.0% and 3.9% -- 4.6x and 6x fewer
      molecules destroyed. Nothing rescues that case: molecules colliding on one reference hold the
      SAME sequence, so neither the position nor payload sub-clustering separates them. On 200 or
      20,000 distinct references they win by 0.001 ARI, and that gap is a DEPTH THRESHOLD rather
      than a limit: collided molecules there carry different sequences and `assemble` separates
      them from the payload with no aligner, but only past the depth 8.68 fixes for itself
      (`log10 C(n, n/2)` clears it at n ~ 32). Measured -- 0/12 separated at 9.1 reads on the
      barcode, 7/14 at 40.9, then 10/10, 13/13, 10/10 at 82, 161 and 283. migec is 8-48x faster
      including the alignment they cannot skip. `docs/grouping.rst`,
      `assets/grouping_tools.tsv`, `assets/collision_split.tsv`
- [ ] `2026-migec-benchmark` repo, `isalgo/umi_data`, comparisons against MIGEC v1, MAGERI,
      Cell Ranger, **UMI-VarCal and UMIErrorCorrect** (the two UMI-aware callers benchmarked by
      Maruzani et al. 2024 for low-frequency ctDNA)
- [x] ctDNA ground truth **located, not simulated** (2026-08-13): `PRJNA788522` (72 runs, cfDNA
      reference material at 0 / 0.125 / 0.25 / 1% VAF x 5/20/80 ng x 3.3/10/30x, three replicates)
      and `PRJNA507366` (28 runs, six polymerases plus 0.031% / 0.0625% VAF). Both carry a real
      12 nt inline UMI. Never: Maruzani's runs carry **no** UMI (aligned BAM submissions, `suggest`
      finds no pattern, `CMP_LINKAGE_GROUP` empty), so theirs had to be simulated — 9 nt, Phred
      fixed at 37, assigned to reads sharing start *and* end, which is the co-terminal assumption
      X1 falsified. That is a property of the two runs they picked, not of the public record
- [x] **LoFreq against the certified arms (2026-08-13/14).** Full chain on GRCh38 with the panel
      inferred from coverage: `assemble` consensus, `minimap2 -y`, LoFreq, scored per target against
      the certified VAF. **Reliable to 0.25%**, and the finding that mattered is that at
      `--min-reads 1` the 2-colour dark-G artifact is **additive to true positives** -- the 0.25%
      arm read 0.79% (0.25% + a 0.57% artifact floor). `--min-reads 3` takes specificity from 0% to
      100% at no cost to sensitivity. `assets/ctdna_minreads.tsv`, `docs/detection.rst`
- [ ] The other three callers: **Mutect2, UMI-VarCal, UMIErrorCorrect** on the same arms, so the
      comparison is caller-vs-caller on identical consensus input rather than migec-vs-published.
      Note: the LoFreq result above is one caller's, and the artifact it exposed is a property of
      the *consensus input*, not of LoFreq -- which is exactly why the others have to be run
- [ ] Re-run the arms with adapter trimming or `-q 20`. Adapter read-through is the diagnosed cause
      of the chr2 mismapping (87S52M at MAPQ 4-16, TruSeq adapter), but the fix is **diagnosed, not
      measured**: nothing has been re-run with it applied
- [x] **MIGEC 1.2.9 head to head** (2026-08-14), `scripts/compare_migec_v1.py`. Same dialect, same
      sheet, same library, both pipelines end to end, `--min-count` matched -- v1 defaults to 5 and
      we default to 1, and v1 names its output `.t5.` for exactly that reason, so leaving each at
      its own default compares defaults. **9.4-11.9x the wall clock** against a 3x gate, 3.5-3.7x
      less memory, and the molecule count is 0.09% over truth against v1's 13.6%: 99.8-99.99% of
      our consensuses are exactly a template, against 93.9-95.1%. `assets/migec_v1.tsv`
- [x] **The published version is single-sourced** (2026-08-14). `docs/conf.py` was a fourth
      hand-written copy of the version and the only one nothing checked, which is how every
      published page said `2.0.0.dev0` for three releases while the wheel was right. It reads
      `pyproject.toml` with `tomllib` now, so a release bumps THREE places and the docs follow.
      Never: it must not `import migec` to get it -- `migec._core` is in `autodoc_mock_imports`,
      so the package's own version-agreement assertion fires against a Mock and the docs build dies
- Gate: grouping ARI ≥0.99 **met** (0.9967-0.9987); residual error ≤1e-5 on a clonal control;
  ≥3× MIGEC v1 wall-clock **met** (9.4-11.9x)

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
