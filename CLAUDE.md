# CLAUDE.md — how to work in this repo

migec extracts UMI/cell/sample barcodes from sequencing reads, corrects errors in them, and
assembles the reads sharing a barcode into a consensus. C++20 core, one pybind11 module, typer
CLI. This file is *how to change the repo*; `ROADMAP.md` is what is implemented, `SOURCES.md` is
where every dataset came from, `project/` holds the design record, `CHANGELOG.md` is per release.

This is a **complete rewrite**. The Groovy MIGEC 1.2.9 is on branch `legacy-v1` / tag `v1-final`;
MAGERI is archived at `mikessh/mageri`. Their algorithms are the specification, their code is not.

## Layout

| path | what |
|---|---|
| `include/migec/` | public C++ headers, one per subsystem |
| `src/` | one `.cpp` per header, plus `_bindings.cpp` |
| `python/migec/` | the Python package: `cli.py` and the thin stage wrappers |
| `tests/cpp/` | doctest, built only with `-DMIGEC_TESTS=ON`, run by ctest |
| `tests/{unit,synthetic,realworld,benchmark}/` | pytest tiers; `synthetic/_sim.py` is the simulator |
| `docs/` | flat sphinx, `formats.rst` is the inter-stage contract |
| `project/` | the design record: six subsystem designs and two critiques |
| `skills/migec/` | the Claude skill shipped with the repo |
| `scripts/` | validation and analysis scripts (`spikein_ratio.py`) |
| `notebooks/` | marimo examples |
| gitignored | `build/`, `.venv/`, `*.mig`, `scratch/` |

## Build, test, lint, docs

```bash
bash setup.sh                                                  # uv venv + editable + assert _core
cmake -S . -B build -DMIGEC_TESTS=ON && cmake --build build -j && ctest --test-dir build
python -m pytest tests/unit tests/synthetic -q
ruff check python/ tests/                                      # pinned to 0.15.9 in CI
sphinx-build -W --keep-going -b html docs docs/_build/html      # zero warnings required
```

## House rules

- **No emojis anywhere.** Not in code, comments, docs, commit messages or any committed text.
  Write the word: "Never", "Note", "Warning", "Done". Mathematical and typographic characters are
  not emojis and are fine. This is a global rule, not a repo preference.
- **`isalgo/umi_data` holds sequences and metadata only** - `.txt`, `.md`, `.tsv.gz`, `.json`,
  `.fastq.gz`, `.fa.gz`, `.sam`, `.bam`. Never reports, figures, logs or pipeline output; a
  derived results table is output, not data, even when it is a TSV. Those live here, next to the
  script that made them. Write through the `~/hf/umi_data` git+lfs mirror, never the HTTP API.

## Non-negotiables

These exist because the alternative was tried, or because a first-pass design got it wrong and the
correction is written up in `project/review-algorithms.md`.

- **Never: Every model-derived number is reported next to something that checks it.** The birthday
  collision rate has `scripts/collision_check.py` (model-free, from the sequences); the barcode
  error estimate has the Phred + polymerase prediction; the emitted quality has the measured RT
  floor. A formula that nothing tests is a formula that will be wrong silently — all three of
  these were, and the check is what found it.
- **Never: The distance-1 barcode-error estimator fails downward as the space fills** (0.92x of truth
  at 0.3% occupancy, 0.23x at 50%, 0.001x at 93%). `err_unreliable` is set past 5% neighbourhood
  occupancy. Never quote the estimate when it is set.
- **Never: Five pipeline commands: `checkout`, `suggest`, `refine`, `assemble`, `subsample`.** A
  sixth *pipeline* command or a new CLI flag requires a *failing benchmark that the default cannot
  pass*. `info`, `sheet` and `plot` are outside the count: they read no reads and write no pipeline
  output. `sort` is not a command and never was -- partitioning lives in `assemble`
  (`project/design-io-interop.md`).
  Constants live in headers next to the measurement that justifies them, so adding a flag means
  deleting a measurement.
- **Never: Never subsample reads to make an example.** Sort by UMI and take all reads of N UMIs
  (`blake2b(umi) mod 10000 < K`). Ten thousand random reads gives one read per UMI and is useless.
  Taking the *first* N distinct UMIs is also wrong — a UMI with 100 reads appears early ~100× more
  often than a singleton, so first-appearance order oversamples large MIGs and destroys the very
  size distribution the example is showing.
- **Never: Do not derate base quality for a low-coverage or suspected-error-child MIG.** If UMI *x* is
  genuinely an error child of *y*, then *all* of x's reads are clean reads of y's sequence: the
  consensus is right and only the molecule count is wrong. Emit the identity posterior as a tag;
  do not smear it into every base.
- **Never: Never emit a quality above the RT/first-cycle-PCR floor.** An error made before
  amplification is in every read and no consensus removes it. `p_floor = ε_RT + 2·ε_pol`, fitted
  from data as the intercept of `e_out(c) = p_floor + a/c`, never hard-coded.
- **Never: A null that permutes reads column by column is the wrong null.** A low-quality read carries
  a minor base at many positions at once, so column permutation hands every read an average error
  load and calls the resulting co-segregation significant. Preserve *both* margins of the
  reads × positions matrix (curveball). Measured: the nominal `p < 0.01` split threshold calls
  30.62% of MIGs, the both-margins null puts the 1% false-positive point 19x higher, at 8.68.
- **Never: A molecule is sample + cell + UMI, never the UMI alone.** UMIs repeat across cells and
  samples by design; grouping on the UMI merges two molecules and nothing downstream can tell.
  `assemble` sorts on `(cell, umi, src_index)` and partitions on the cell when there is one.
- **Never: Contig assembly means one molecule's fragments, nothing more.** Random priming leaves reads
  that tile a molecule under one barcode; `--contig` places them and emits one consensus per
  overlap component. Assembling a cell's full-length receptor, calling doublets and filtering
  contaminating chains are **arda's** job. Note: It also needs a barcode that is not saturated: two
  fragments of two different molecules on one barcode have no sequence in common, which is exactly
  what two fragments of one look like. Report `expected_molecules_per_group` and warn.
- **Never: 1-3 reads per UMI is the normal case, not the exotic one.** Bulk repertoire profiling and
  shallow 3' GEX both look like this. Never threshold it away -- `--min-reads` defaults to 1 and
  the answer to a barcode error is correction, not a cut. Three things calibrated on a deep library
  stop applying and must be said rather than quoted: the split threshold is inert (needs ~30
  reads), the count-ratio error-child null has no dynamic range, and singleton filtering is
  unaffordable (79% of barcodes, against 56% on a deep library). Benchmarks use the shallow shape
  because distinct barcodes are what everything scales with.
- **Never: A collision rate is estimated with the U-statistic, never the plug-in.** `sum p_hat^2` is
  biased up by `(1 - sum p^2)/n`, and the bias GROWS as the distribution spreads -- so on k-mers it
  grows with k and reads as dependence accumulating with k. It manufactured a 1.007x excess out of
  independent data. Use `sum n_a(n_a-1)/(N(N-1))`.
- **Never: An `N` is not a fifth base.** Fold it to A, as `pack_barcode` does -- the packed key is what
  every stage groups on. Counting it as a letter let `m_j` fall to 0.2466, below the mathematical
  floor of 1/4, and printed an effective length of 9.01 nt for a 9 nt barcode. `m_j >= 1/4` is the
  check; use it.
- **Never: Range partition, never hash partition.** A hash sends a barcode and its 1-mismatch
  neighbours to uncorrelated buckets, which makes correction impossible to apply and splits the
  molecule permanently — and the halves each look like a well-formed MIG, so nothing detects it.
- **Never: No indels anywhere.** Substitutions only. Reinstate when there is a dataset to verify against.
- **Never: Nothing that scales with the input on the serial path, and no hash map keyed by barcode.**
  zlib level 6 does 7 MB/s on random DNA, so compression belongs on the workers (concatenated gzip
  members are a valid stream) at level 1. A `unordered_map<uint64,uint32>` costs ~48 B per distinct
  UMI against 22 for a sorted array, which at NovaSeq scale is 19 GB against 8.8. Both were
  measured; `tests/benchmark/` guards them.
- **Never: `-t` must never change the output.** Chunks are matched in parallel and written in input
  order. A demultiplexer whose output depends on its thread count produces results that cannot be
  compared between runs, and the failure is invisible.
- **Never: One output file per *sample id*, never per barcode-table row.** Rows sharing an id are a
  sample sequenced with more than one tag — the format's own idiom. A file per row opens the same
  path twice and interleaves two `FILE*` into it, which is not a valid gzip stream, and the summary
  still reports success. Anything else keyed by sample (UMI counters, summary rows) groups the same
  way.
- **Never: Nothing may throw out of a worker thread.** An escaping exception is `std::terminate`:
  SIGABRT, no message, no flush. Workers capture, the driver rethrows. Validate at the boundary
  where the error is attributable — a bad pattern is caught when the pattern is compiled, on the
  caller's thread, not when a read is packed on a worker.
- **Never: A reported clock covers the whole operation.** `wall_seconds` stopping at the parallel driver
  hid a serial stats pass worth 4/5 of the run and published a throughput 5× what a user sees. If a
  stage is serial, time it and report it separately — that is the number that says what to fix next.
- **Note: No `log2`/`exp`/`pow` in a per-base loop.** The score depends only on the reported Phred and
  the IUPAC set size, both small integers; it tabulates into 1.2 kB. The transcendental was 90% of
  checkout's runtime before it did.
- **Note: Strand normalisation happens in `checkout`.** The `.mig` flags describe what has *already*
  been applied. A group containing both orientations silently loses half its reads in consensus.
- **Note: Nominal Phred is not the error rate.** On 2-colour instruments there are ~4 distinct Q
  values. `checkout` measures `ê(q)` against the pattern's own constant bases and fits
  `ε_qi + a·10^(−q/10)`; use the calibrated table wherever a likelihood is computed. Measured
  slope on SRR1763769: 1.04 over 46.3 M bases.
- **Never: The calibration intercept is the PRIMER's defect rate, not a sequencing floor.** The
  standard is a synthesised oligo and synthesis runs ~1 defect per 200-500 bases; the fitted
  intercept is 3.9e-3, spread evenly over all 23 anchor positions with none polymorphic, and it
  agrees with the independently measured 0.55% one-base-short rate from failed couplings. Folding
  it into `error()` would add 4e-3 to every base likelihood in the pipeline on the strength of the
  primer's quality. Report it, never apply it.
- **Note: Do not reproduce MIGEC v1's bugs.** Quality must be indexed at the *match offset*, not the
  read start; low-quality mismatches must actually be counted (v1's dangling `else` meant they
  never were).
- Done: The version string is duplicated in `pyproject.toml`, `include/migec/version.hpp` and
  `python/migec/__init__.py`. `__init__.py` asserts the last two agree and CI prints all three.

## Release process

1. Bump the version in all three places above; add a `CHANGELOG.md` section.
2. Branch `release/<version>` → PR → CI and docs green.
3. `gh release create v<version>`.
4. `publish.yml` builds wheels and publishes via OIDC to the `pypi` environment. Never: Do not rename
   that workflow file or the environment: the PyPI trusted publisher is bound to both.

## Open loops

- M0 done. `migec checkout` works: patterns, trimming, header transfer, UMI statistics, count
  correction, **paired input with strand normalisation**, and **multi-core with byte-identical
  output at any `-t`** (1.06 M reads/s end to end, 1.68 M matching, at 16 threads;
  `scripts/benchmark_threads.py` writes the table the figure is drawn from). Whitelists and
  dual-end barcodes are done; `.mig` bucket output from checkout is still open.
- **Note: The UMI counters are not partitioned yet.** ~22 B per distinct UMI is 8.8 GB at NovaSeq
  scale, held in one piece. The fix is the range partition (M2, with `.mig` bucket output), not a
  smaller struct. Until then checkout warns past 1 GB.
- **Note: The per-sample statistics are the serial tail.** Histogram, composition and `correct_umis`
  run once at the end on one thread, ~1.5–2 µs per distinct UMI, so at 16 threads they are a larger
  wall than the gzip reader. `correct_umis` in particular is O(n · 3L · log n) and arguably belongs
  in `refine` (M3) rather than in checkout at all — decide that when `refine` lands, and until then
  it is why `match_seconds` is reported next to `wall_seconds`.
- Grouping accuracy vs Calib is wired up (`scripts/compare_calib.py`, `docs/grouping.rst`) and the
  migec column is asserted in CI. Running the Calib column needs Calib built locally.
- **X1 is done (2026-08-13) and it answered yes: fragmented mode is mandatory.** Reads sharing a
  `(CB,UMI)` are not co-terminal (7.8% overall, 0.3% at ≥6 reads), but 72.7% still form one
  overlap component, so `assemble` partitions by overlap and then runs the ordinary ungapped
  consensus per component. Never extend a component across a gap. `docs/fragmented.rst`.
- **X2 is done (2026-08-13): the RT/PCR floor is of order 1e-4, so cap emitted quality at Q40.**
  Measured 1.54e-4 [1.36e-4, 1.74e-4] on `SRR1763769`, which agrees with 10x's stated 1e-4 for the
  V(D)J RT. Never: it is the ONE-MOLECULE floor and every record we emit is one molecule -- 10x's
  Q60 is for bases covered by >=2 UMIs, and combining molecules is arda's job. `--rt-error` names
  the class instead of guessing: `rt` 1e-4 (default, Q40), `medium` 1e-5, `high` 1e-6, or the rate
  (TSO500 v2 is 7.37e-5). `--rt-error auto`
  still fits per dataset (the floor is a property of the enzyme and cycle count), but the
  **default is 1e-4**. Note: Still an upper bound — that library is 49.6% occupied on a 9 nt barcode
  and `checkout` calls it saturated, so collisions inflate it. `docs/quality_floor.rst`.
- **`migec suggest` is implemented** (2026-08-13), ahead of its M4 slot, because X2 needed it: it
  segments the per-cycle base composition into UMI (all four bases near 1/4), constant and payload
  runs and prints a paste-ready pattern. Recovered the 9 nt + `CAGTTTAACTTTTGGGCCAT` layout of
  SRR1763769 unaided. Note: It stops the pattern at the last *constant* run — composition alone cannot
  tell a UMI from diverse payload, only the anchor can.
- **The barcode-space and error-budget arithmetic is built in**, logged to
  `checkout.barcode_space.tsv` / `checkout.umi_quality.tsv`, warned on, documented in
  `docs/barcode_space.rst`, drawn in `notebooks/barcode_space.py` and tested in
  `tests/synthetic/test_barcode_space.py`.
- **X3 is done (2026-08-13).** Position independence holds (1.0103x over 9 nt; the first pass said
  1.04x and all of it was N-as-a-fifth-base artefact), so `Π_j m_j` stays
  and the 1.86x collision excess is the read threshold, not the barcode. 92% of distance-1 pairs
  are chance at 47.8% occupancy; the column-shuffle background puts barcode error at 3.4e-3 —
  within 1.7x of the Phred + polymerase prediction, where the analytic estimate is 2.6x below it,
  so **M3 takes the permuted background**. The split threshold is **8.68, not 2.00**.
  `docs/nulls.rst`, `scripts/permutation_nulls.py`, `tests/synthetic/test_nulls.py`.
- **M1 is done (2026-08-13): `migec assemble` works.** `(cell, umi)` grouping, range partition
  into `.mig` buckets with one bucket resident, the column log-likelihood posterior, the RT floor
  added (not compared) so nothing above the named floor is emitted (Q40 by default), linkage sub-clustering at X3's 8.68, and
  `--contig` for random-primed reads (seed placement, union-find overlap components, never bridged
  across a gap). 531 k reads/s; 121 MB at 16 buckets against 203 MB at one; asserted in
  `tests/benchmark/test_assemble_speed.py` in a fresh process per configuration, because
  `peak_rss_bytes` is a process high-water mark.
- **Note: The writer buffer budget is split across buckets, not per bucket.** A fixed per-writer block
  made pass 1 cost grow with the bucket count, which is backwards — more buckets exist to use less
  memory, and cutting finer made peak RSS go *up* (238 MB at 16 buckets against 213 at one). It is
  now `clamp(32 MB / buckets, 256 kB, 4 MB)`.
- Note: A tail quantile is not a constant until its Monte Carlo error is smaller than the digits
  quoted.** X3's split threshold read 9.91, then 9.61, then 11.66 on reruns of ~8,000
  randomisations. At 82,800 it is **8.68, bootstrap 95% CI [8.42, 9.14]**, and that interval is
  what the code uses and the docs quote.
- **Never: The count ratio is not evidence below ~3 reads/UMI**, and that is the common regime. Two
  evidence terms that survive at one read are now in `BarcodeEvidence`: the barcode's own **base
  quality** at the differing position, and **payload agreement** with the candidate parent, worth
  `log(1/clonality)` where the clonality is measured from random barcode pairs rather than assumed.
  Payload agreement lifts the count gates (which is what makes a singleton merge possible) and
  payload *disagreement* refuses a merge the counts would have made. Never tune
  `max_child_fraction` to paper over this — it is the wrong evidence, not the wrong threshold.
- **Never: Compare a rate to a rate.** The error likelihood was a zero-truncated Poisson, i.e. a
  probability conditioned on the child existing, weighed against `a_ind * p_size`, which is an
  expected *count*. The truncation divides out `(1 − e^−λ)` — exactly the term that says whether an
  error child should exist — so ZT-Poisson(1, λ) → 1 for every small λ and the error rate stopped
  mattering at precisely the coverage where nothing else was available. Untruncated.
- **Never: Err on precision, never on recall, when merging.** A wrong merge deletes a molecule and
  nothing downstream can tell; a missed correction only inflates the count and is recoverable. The
  current posterior destroys **no** real molecule at any depth measured (≥0.99, 1.000 at the
  extremes) and that is the property to protect.
- **Note: At ~1 read/UMI, 80% of barcode errors are unfixable in principle** — the parent barcode was
  never sequenced, so there is nothing to merge into. Always report recall against that reachable
  ceiling; against all children it understates by 5x and looks like a bug.
- **`migec refine` works (2026-08-13).** Barcode table, correction with all three evidence terms,
  read rewrite with `OX:Z:` preserving the original, barcode table TSV. 20,055 molecules recovered
  from 20,000 simulated, ε at 0.96x of injected. It holds the **table**, never the reads, and
  streams three times.
- **Note: Correction is not bucketable by a plain range partition.** The top b bits of the key decide
  the bucket, so a substitution in the top b/2 positions sends a barcode and its neighbour to
  different buckets and the pair can never be found. Two passes with the key rotated fixes it.
  Until then the table is whole and its size is reported, as checkout does with its counters.
- **Cell calling works (2026-08-13)**: molecules per cell (never reads), OrdMag with the knee
  reported beside it and a warning past 3x disagreement. Exactly 500 of 20,500 on a synthetic
  droplet library. Note: The cell key is the TOP `2*cell_length` bits of the packed barcode, not
  "everything above the UMI" -- those coincide only when cell+UMI fill all 32 bases, and getting it
  wrong silently shatters cells into fragments that still look like cells.
- **Whitelists work (2026-08-13)**, `refine --cell-whitelist`. Never: The posterior needs the
  background hypothesis or every hopped/undeclared barcode is absorbed into its nearest entry with
  posterior 1.0. Note: That prior is per **barcode**, not per library: the whitelist prior is
  `1-background` spread over every entry (~1e-6 each for 737k), so a library-level "1% off-list"
  is four orders of magnitude too big and wins every time. Use the off-list read share divided by
  the distinct off-list barcodes, both measured.
- **Never: Any "is this a child?" test built on the count ratio reports ZERO at 1-3 reads/UMI**, which
  is where the answer matters most. It caught me twice: once in the correction posterior, once in
  the residual-FDR estimator. Both now use payload agreement as well. On a 1.23 reads/barcode
  library the count-only residual is 0 and the full one is 1,294.
- **The MIG-size FDR threshold is reported, never applied** (2026-08-13). Measured residual, not
  derived. 5.25% of 1-read molecules at 1.23 reads/barcode, 0% at 4.62.
- **Dual-end barcodes work (2026-08-13)**: column 3 of the sheet is the slave pattern, on the
  other mate, extending the UMI. Never: Both halves or nothing -- accepting the master alone emits
  12 nt UMIs beside 24 nt ones and every collision estimate is then over two spaces at once.
- **Never: The acceptance bar is charged for the offsets ACTUALLY SCANNED.** It was using the offsets a
  read could hold, so `--max-offset 0` still paid `log2(61/α)` = 12.6 bits and a 5 nt dual-end
  handle (10 bits) was refused on every read. A free scan refusing it is correct -- `TGACT` occurs
  by chance every kilobase -- but an anchored one must not.
- **Never: Never take the first N records of a stage's output.** `assemble` and `refine` write in
  barcode order, so the first N share their leading bases: sampling the first 4,000 consensuses of
  SRR1763769 reported a UMI effective length of **6.45 nt against the true 8.97**, with positions
  0-1 apparently 100% A. Reservoir-sample. This is the same error as subsampling reads instead of
  whole barcodes, in a third costume.
- **Never: Payload agreement is worth nothing on a clonal library**, and anything using it must
  discount by the measured clonality. `correct_umis` does; the residual-FDR estimator did not, and
  reported 97.4% of singletons as error children on an HIV amplicon whose clonality is 0.80.
- **`scripts/diagnose.py`** answers what a run actually is: places the consensus on a reference
  (SRR1763769 = HXB2 2,328-2,595, minus strand, protease/RT), the barcode's own PWM, flowcell
  coordinates and the i7xi5 table when the headers carry them, and the minor-allele spectrum split
  head-vs-tail so read-end decay is not reported as a quasispecies.
- **Positional is the primary mode (2026-08-13, 2.0.0a2).** `^NNNN` / `^NNNNXNNN`, or half-open
  slices `0:8` / `0:4,5:10` / `cell:0:16,16:26`, both in `--bc-pattern`. **Never: `--max-offset` is
  automatic now and must not be passed** — a caret, a slice list, a read structure and a pattern
  with nothing to score all anchor at 0. Passing `-1` reinstates the old refusal, which is still
  correct. Presets in `sheet.PRESETS`: umi, migec, primerid, duplex, 10x, 10x-v2, tso500,
  smarter-umi. Never: every preset carries a citable source and a test that compiles it; a preset
  nobody can run is worse than none, because it looks supported.
- **The downstream contract is measured (`docs/downstream.rst`).** minimap2 `-y`, bwa `-C`, arda,
  salmon, kallisto all run against real `assemble` output; 600/600 records keep `RX`/`CB`/`MI`
  through a sorted BAM, and arda's AIRR `sequence_id` **is** the molecule id. Note: that is why the
  name is `<sample>.<cell>.<umi>` — `dnaio` drops FASTQ comments, so the name must stand alone.
  Never: alevin/bustools/STARsolo must not see a consensus FASTQ; they deduplicate from a raw
  barcode read that no longer exists. STAR unverified here — the brew arm64 build reads 0 reads
  from any FASTQ, including a one-record file, so it says nothing about our output.
- **SRP150352 (UMI RNA-seq, Sci Rep 2018) cannot be reprocessed from SRA** — `ncgr/UMI-analysis`
  moves the UMI into the FASTQ header and SRA rewrites headers. Confirmed on three runs by the
  missing template-switch `GGG`. `suggest` reports it correctly. The `smarter-umi` preset is
  sourced from the pipeline, not fitted to the data. In `SOURCES.md`.
- **`migec plot` works (2026-08-13, 2.0.0a3)**: fifteen gnuplot panels over the TSVs the stages
  already write, ColorBrewer Dark2, A/C/G/T always the same four colours. Never: it computes nothing
  — a figure that cannot be redrawn from a committed table is a figure that will disagree with the
  report. gnuplot is **not** a Python dependency; without it the `.gp` scripts are still written,
  which is also how the tests exercise that path (`run(..., gnuplot="")`). `assets/` holds the
  graphviz pipeline figure and the example panels the README embeds, each next to the table it was
  drawn from.
- **The RT floor is named per protocol, not fitted (2.0.0a3).** `--rt-error rt|medium|high|<rate>`
  = 1e-4 / 1e-5 / 1e-6. Never: **it is the one-molecule floor.** 10x give Q40 to a base covered by one
  UMI and Q60 only to a base covered by two or more, because an RT error is common-mode within a
  molecule and independent between them. Every record we emit is one molecule, so Q40 is the
  default cap and the Q60 case belongs to arda. Every value is cited in `SOURCES.md` (10x's V(D)J
  figure, X2's own 1.54e-4, TSO500 v2's 7.37e-5, McInerney 2014 for the polymerases, Shagin 2017
  for the 5x first-cycle factor).
- **`--fast` is counting mode (2.0.0a3)**, not a faster consensus: modal exact sequence, per-base
  max quality over the reads carrying it, floor still applied. Refused with `--contig`. Never: the
  max is over the reads carrying the WINNING sequence — across variants it would take its highest
  quality from the reads that disagree.
- **Coverage into the consensus is capped at 10,000 reads/barcode (2.0.0a3).** Never: the cap is on
  the reads consensed, never on the reads counted; `cD` stays the true depth or the abundance of
  the most-amplified molecules is silently flattened.
- **Note: the benchmark corpus was assigning every read to one sample of four** (`i % 4` with four
  reads per molecule is always 0). The matcher scored all four patterns, so throughput was sound,
  but the per-sample counters and the memory figure were one sample's. Fixed; the published table
  is re-measured and now comes from `scripts/benchmark_threads.py`, which writes the TSV the
  figure is drawn from.
- Next: M3's remainder (the template's own error split), then the M2 remainder
  (`.mig` bucket output from checkout, i7xi5, bit-parallel matcher), then `--rt-error auto`.
- **Note: Britanova et al aging (bulk TCR, shallow) lives on aldan3** and is the real dataset for the
  1-3 read regime. Not pulled yet. aldan3 compute goes through SLURM, never the frontend.
- The archive is pushed: `legacy-v1` + tag `v1-final`, and master is the rewrite. Recovery point
  `~/backup/migec-local-mirror-2026-08-13.git`. Note: The canonical repo is **antigenomics/migec**;
  `mikessh/migec` is a redirect, and `gh` commands must use the former.
- Note: `gh repo edit --default-branch` and force-pushes may be blocked by the permission classifier;
  hand those to the user rather than working around them.
- `isalgo/umi_data` **exists and is populated** (2026-08-13): the CI fixture and the derived
  result tables, 1.3 MB. Note: It is **public**. Written through the mirror at `~/hf/umi_data`, never
  the HTTP API. Nothing from aldan3 is in it, and raw ENA runs are not either -- `SOURCES.md`
  carries the regenerating command instead of the bytes.
- Note: Do **not** add seqtk as a dependency. It is the right tool for generic FASTQ slicing in a
  benchmark harness, but it cannot subsample by whole UMI (it samples reads), which is the one
  thing we need and the thing that makes a UMI example fixture correct.
