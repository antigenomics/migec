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
- Done: The version string is in FOUR places -- `pyproject.toml`, `include/migec/version.hpp`,
  `python/migec/__init__.py` and `docs/conf.py`. `__init__.py` asserts the middle two agree and CI
  prints all three; nothing checks `docs/conf.py`, which is how it sat at `2.0.0.dev0` for three
  releases while every published page showed it.

## Release process

1. Bump the version in all FOUR places above; add a `CHANGELOG.md` section. Then check the README:
   it is the PyPI long description, so a repo-relative image or link is broken from the moment the
   wheel lands.
2. Branch `release/<version>` → PR → CI and docs green.
3. `gh release create v<version>`.
4. `publish.yml` builds wheels and publishes via OIDC to the `pypi` environment. Never: Do not rename
   that workflow file or the environment: the PyPI trusted publisher is bound to both.

## Open loops

- M0 done. `migec checkout` works: patterns, trimming, header transfer, UMI statistics, count
  correction, **paired input with strand normalisation**, and **multi-core with byte-identical
  output at any `-t`** (1.06 M reads/s end to end, 1.68 M matching, at 16 threads;
  `scripts/benchmark_threads.py` writes the table the figure is drawn from). Whitelists,
  dual-end barcodes and `.mig` bucket output are all done.
- **The UMI counters bound themselves and correction follows them (2026-08-14).** Past
  `umi_budget_bytes` (1 GB per run, divided by the samples; a `checkout.run()` kwarg, never a CLI
  flag) each counter range-partitions into `<out>/.umi_spill`, removed when the summary is written,
  and the histogram, composition, distinct count, distance-1 census and correction all stream a
  bucket at a time. Never: **a partition alone bounds the memory and silently stops correcting** --
  the bucket is the top bits of the key, so an error in the first `(bits+1)/2` positions puts a
  barcode in a different bucket from its parent forever. Correction runs **two passes**: over the
  buckets as they stand owning positions `[pb, L)`, then over a copy with keys rotated left by `pb`
  owning `[L-pb, L)`, so every pair is weighed once and `merged` counts barcodes, not opportunities.
  Never: a bucketed run answers with the SCALARS -- `root`/`corrected` are indexed against
  `entries()`, which is the array being bounded, and `BarcodeEvidence` is refused rather than
  ignored. Verified field for field against the resident run on a simulated library and a 500 k
  corpus. Costs ~2.2x wall when it fires (718 k -> 333 k reads/s), nothing when it does not. Two
  bugs it found: the FIRST flush swapped the buffer in and returned before the budget check, so a
  library arriving in one buffer never partitioned at all; and the driver reported the census's 0.0
  where the resident path reports the 1e-4 floor it actually corrects at.
- **Note: The per-sample statistics are the serial tail.** Histogram, composition and `correct_umis`
  run once at the end on one thread, ~1.5–2 µs per distinct UMI, so at 16 threads they are a larger
  wall than the gzip reader. `correct_umis` in particular is O(n · 3L · log n) and arguably belongs
  in `refine` (M3) rather than in checkout at all — decide that when `refine` lands, and until then
  it is why `match_seconds` is reported next to `wall_seconds`.
- **The whole chain runs on buckets (2026-08-14, 2.3.0): checkout --mig -> refine -> assemble.**
  refine reads buckets and writes them back **re-partitioned on the CORRECTED barcode** -- a
  corrected barcode is a different key and a key decides its bucket, so copying a bucket through
  unchanged stops it being a partition and the reads corrected across a boundary get grouped with
  strangers. Every number matches the FASTQ route (`tests/synthetic/test_mig_chain.py`): reads,
  barcodes, merges, molecules, the error to nine digits, the barcode table byte for byte, the
  consensus byte for byte. Note: the audit trail on this route is `<sample>.barcodes.tsv` -- a
  `.mig` record has no room for the pre-correction barcode the way a comment has `OX:Z:`.
- **`.mig` is at v2, and a v1 file still reads.** v2 carries the BARCODE's own quality, one Phred
  per base, as two fixed-width columns. Never: v1 stored `umi_minq`/`cell_minq`, the MINIMUM over
  the barcode, and the posterior weighs the quality AT THE POSITION THAT DIFFERS -- a minimum says
  every position is as bad as the worst, which overstates the error everywhere and makes merges
  easier, the wrong direction. A v1 file comes back with empty columns and refine falls back to the
  global rate, as it does for a FASTQ with no `QX:Z:`.
- **Never: an adversarial audit of new code is worth more than another test of old code.** The
  `.mig` work passed every test it had and still carried three data-integrity bugs, all reproduced:
  refine wrote its output buckets OVER its input (same names, `"wb"` truncates); refine's bucket
  rewrite ignored `--limit-*` and emitted uncorrected reads beside corrected ones; and the sample
  id -- which in assemble and refine comes from the DATA, the `BC:Z:` tag -- went into every output
  path unvalidated, so `BC:Z:../..` wrote outside the output directory. `validate_sample_id` is the
  guard, in types.hpp, called where the id is attributable.
- **`checkout --mig` writes the partition `assemble` was building (2026-08-14).** Same key, same
  range partition, `<sample>.<bbb>.mig`, one writer per (sample, bucket) owned by one thread for
  the whole run so `-t` still changes nothing but the clock. 500 k reads, four samples, `-t 4`:
  **1.16 s -> 0.98 s** end to end, identical 124,878 molecules, consensus FASTQ byte-identical
  after decompression. Note: the open-file budget is for the RUN -- 256 buckets on one sample, 64
  each on four, two each on a 96-plex, because that sample also holds a 96th of the reads. Never:
  **only buckets assemble WROTE are deleted after use** -- pass 2 removes each bucket as it
  finishes, which ate three of four samples' partitions the first time `--mig` was run end to end.
  Never: two samples' buckets in one `assemble` are refused by name, and `--limit-*` on a
  partitioned input is refused -- a limit is a prefix and a partition has none. Opt-in; FASTQ
  stays the default because a `.mig` file is a migec intermediate nothing else reads.
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
  (TSO500 v2 is 7.37e-5). `--pre-amp-error auto` fits it per dataset (the floor is a property of
  the enzyme and cycle count), but the **default is 1e-4**. Note: Still an upper bound — that library is 49.6% occupied on a 9 nt barcode
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
- **The barcode error rate is measured at every DEPTH, not once for the library.**
  `<sample>.umi_errors.tsv`, one row per exact parent depth: a parent carrying `c` reads offered
  `c*L` barcode bases to be miscalled, so `u(c) = 3L(1-exp(-c eps/3))` from the distinct children
  and `r(c) = c L eps` from the reads in them invert to the same eps. Reported as
  `error_at_depth` (parents seen >= 10 times, where correction is near-complete), `error_phred`,
  and `error_from_children` (all depths, a LOWER bound). On 1e-3 injected: distance-1 9.73e-04,
  children 9.89e-04 = Q30, all depths 9.98e-04.
- **Never: neither error estimator is saturation-free — both are bounded by the merges correction
  MADE.** As a fraction of an injected truth: distance-1 0.97 / 0.96 / 0.76 / 0.45 / 0.001 and
  children 0.99 / 0.95 / 0.88 / 0.62 / 0.00 at occupancy 0.2% / 2.3% / 9.8% / 33% / 100%. The
  children estimate wins everywhere either works; at 100% BOTH are zero, because `correct_umis`
  refuses to merge on a full space and is right to. `saturated` is the flag that says the answer
  is a floor -- do not read the table instead of it. `docs/umi_errors.rst`,
  `tests/synthetic/test_umi_errors.py`.
- **Never: a sparse spectrum is POINTS, never a line.** `mig_size_spectrum` drew "reads in them"
  `with lines` over a table at EXACT sizes, where past the head almost every size holds one
  molecule -- so reads == size and the line drew the `y = x` diagonal as the figure's most
  prominent feature, a tautology that reads as a second mode. Two or three molecules on a size made
  it sawtooth between `size*1` and `size*2`, and it bridged gaps where nothing was observed. Same
  rule for the two `umi_error_*` panels.
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
- **The downstream contract is measured (`docs/downstream.rst`).** minimap2 `-y`, bwa `-C`,
  **minibwa `map -y`**, arda, salmon, kallisto all run against real `assemble` output; 600/600
  records keep `RX`/`CB`/`MI` through a sorted BAM, and arda's AIRR `sequence_id` **is** the
  molecule id. Note: that is why the name is `<sample>.<cell>.<umi>` — `dnaio` drops FASTQ
  comments, so the name must stand alone. Note: minibwa spells the flag `-y` on `map` (minimap2's)
  and `-C` on the legacy `mem` (bwa's); each rejects the other with a non-zero exit, so the tags are
  never dropped quietly. Never: alevin/bustools/STARsolo must not see a consensus FASTQ; they
  deduplicate from a raw barcode read that no longer exists. STAR unverified here — the brew arm64
  build reads 0 reads from any FASTQ, including a one-record file, so it says nothing about our
  output.
- **Map-first vs collapse-first is written up (`docs/downstream.rst`, 2.0.0a4).** fgbio, UMI-tools
  and UMIErrorCorrect (Österlund 2022, `SOURCES.md`) align raw reads and group on *(position, UMI)*;
  we group on *(sample, cell, UMI)* and align once. Note: the position is worth real key bits when
  the barcode is short (TSO500's 5 nt is 1,024 barcodes) and close to zero on a repertoire library,
  where everything maps to the same V genes. Note: what `assemble`'s linkage sub-clustering at 8.68
  recovers is exactly that discriminating power, from the payload, with no aligner. Note: the
  dividing line for downstream tools is **transport vs deduplicate** — a tool that carries `RX`
  composes with migec, a tool that dedups on it replaces a stage of migec.
- **SRP150352 (UMI RNA-seq, Sci Rep 2018) cannot be reprocessed from SRA** — `ncgr/UMI-analysis`
  moves the UMI into the FASTQ header and SRA rewrites headers. Confirmed on three runs by the
  missing template-switch `GGG`. `suggest` reports it correctly. The `smarter-umi` preset is
  sourced from the pipeline, not fitted to the data. In `SOURCES.md`.
- **The i7xi5 contingency table is in (2026-08-14, 2.3.0).** Read off the instrument's header for
  EVERY read, matched or not -- restricting it to assigned reads would hide the population it
  exists to measure. A combination is "declared" when it holds >= 5% of the reads of its own i7 AND
  its own i5; the sheet carries the in-line barcode rather than the index pair, so the declared set
  is inferred, and the gap is wide (hopping 0.1-2% against a declared combination being the bulk of
  its index). Never: a SINGLE-indexed run is `estimable = false`, not a rate of zero -- nothing can
  be off-diagonal when there are no combinations.
- **The four familiar QC figures are in (2.1.0)**: Cell Ranger's barcode rank plot
  (`<sample>.cell_rank.tsv`), the MIG size spectrum on log1p with molecules AND reads, the
  rank/Zipf curve, and unique UMIs per sample barcode. Never: the knee plot's y axis is **unique
  UMIs, never reads** -- one over-amplified molecule otherwise puts an empty droplet high on the
  curve, which is the artefact the plot exists to show. Never: `<sample>.sizes.tsv` is at EXACT
  sizes, never power-of-two bins -- the Zipf curve is its cumulative count and four bins make four
  steps.
- **Never: a thinned scatter is not a distribution.** `consensus_quality` drew `every 17` molecules
  as grey dots. Emitted quality is DISCRETE and capped at the floor, so at any real depth every
  molecule sits on one or two integers: the cloud drew a flat line whether the bin held ten
  molecules or ten million, and the thinning threw away the tails that were the only thing it could
  have shown the line did not. `assemble` accumulates the exact (depth bin, rounded Phred) grid --
  61 counters per bin -- and the panel is a candlestick over real order statistics.
- **Publication defaults on every panel (2.1.0)**: transparent background so one SVG serves a light
  README, a dark README and print; one ink colour `#808080` that reads on both; **the key inside the
  plot box**, because a legend gutter makes every figure wider than its data. Frame 760x520. The
  pipeline figure is `rankdir = TB` with three `rank = same` groups, not a 5:1 strip.
- **`migec plot` works (2026-08-13, 2.0.0a3)**: twenty gnuplot panels over the TSVs the stages
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
- **All three stages thread (2026-08-13, 2.0.0a3-2.1.0).** refine 222k -> 1.01M reads/s, assemble
  203k -> 2.05M (its partition threaded in a4), checkout unchanged at 1.06M. Never: the FIRST fix was not a thread -- zlib level 6
  was 83% of refine's wall clock compressing an intermediate the next stage immediately
  decompresses. Level 1 is now the default for every stage's output and gave 3x on its own.
  Measure before parallelising.
- **Never: the thread HELPER was the bottleneck once everything else was threaded (unreleased).**
  `parallel_for` claimed one item per atomic `fetch_add`, and a sampling profile of assemble put
  **21% of all CPU samples across every thread on that one instruction** -- more than the parse it
  was handing out. Sixteen cores serialise on one cache line. Batched now: ~8 turns per worker,
  capped at 1024, collapsing to 1 when items are few, which is the uneven case (one bucket per
  item) the counter exists for. Never: profile the helper, not only the work. Current numbers,
  1 thread -> 16: checkout **213,880 -> 1,548,835** end to end (1,697,313 matching), refine
  **617,802 -> 1,554,156**, assemble **554,106 -> 2,470,928**; `assets/benchmark_threads.tsv` is
  the committed table and the figure is drawn from it.
- **Two serial scans went with it**: the distance-1 census in `estimate_umi_error` (which IS
  checkout's per-sample statistics tail -- the end-to-end/matching gap fell from 20% to 9%) and
  refine's residual-FDR scan, 0.53 s of a 2.17 s run on one core. Both are read-only `3L` binary
  searches per barcode tallying integers, so per-worker counters summed afterwards keep `-t`
  changing nothing but the clock. `estimate_umi_error` now takes a `threads` argument, defaulted
  to 1 so the C++ test calls are unchanged.
- **Note: a 64 k chunk in assemble's partition is STILL 32% faster and still not taken**
  (3,075,506 reads/s against 2,324,403 on 4 M, re-measured after the batching). It costs 16 MB of
  resident chunk, which at NovaSeq scale makes pass 1 the memory peak. Never: the benchmark tier
  does NOT fail at 64 k on a 500 k corpus -- the objection is scale, which no test here can see,
  so the number lives next to the constant in `assemble.cpp`. Upgrade path is a persistent worker
  pool, not a bigger chunk.
- **Never: `clear()` on a chunk of `FastqOwned` is four allocations a record.** assemble's
  partition knew this; refine's rewrite still did it. Assigned into now, and the rewrite's comment
  buffer and unpacked barcode are worker scratch with the cell/UMI halves as `string_view` rather
  than `substr` -- `rewrite_seconds` 0.51 -> 0.38. Same class of bug in `assemble`: each group's
  UMI was unpacked twice, once for the composition tally and once for the record name.
- **Never: `-t` changes nothing but the wall clock, on every stage.** refine parallelises the
  neighbourhood SCAN (a pure function of the barcode table -- it reads no union-find state) and
  applies the merges serially in the original order, so the result is identical, not equivalent.
  assemble gives each worker a bucket and concatenates in bucket order. Never: the bucket count is a
  fixed floor of 16 and must NOT depend on `--threads` -- if it did, `-t` would choose the gzip
  member boundaries and two runs would differ byte-wise while holding identical records.
- **Note: the race audit is a build, not a claim**: `cmake -S . -B build-tsan -DMIGEC_TESTS=ON
  -DCMAKE_CXX_FLAGS="-fsanitize=thread -g"`. 104 cases, 224,116 assertions, no race -- and the
  instrumentation was proven to fire on a deliberate unsynchronised counter in `parallel_for`
  before the clean run was believed.
- **`--limit-read` on every stage, `--limit-umi` on refine and assemble (2.0.0a3).** Note: checkout
  has no `--limit-umi` and should not — counting distinct UMIs there means matching first, which is
  the work being limited. Never: a limit is not a sample and the report says so -- the first N reads
  are one corner of one flowcell. `subsample` is the sampler.
- **Note: `migec sort` was asked for and measured instead.** A separate sort adds a whole read+write
  to skip a pass that is now 1.45 s of a 1.95 s run at 2.05 M reads/s, and you pay it again on every
  input. It stays unexposed until a benchmark says otherwise.
- **assemble's partition threads too (2026-08-13, 2.0.0a4).** It was 2.07 s of a 2.69 s run against
  a 0.23 s `gzip -dc` floor, so five sixths of it was not the inflate. Now: **1,481,946 ->
  2,051,937 reads/s, and peak RSS 1,479 -> 789 MB.** Never: **ownership, not locking** -- worker w
  owns every bucket with `bucket % threads == w` for the whole run, so a bucket file has exactly one
  writer and no bucket state is shared; ownership decides *who* writes a record, never which file or
  where, which is why the bytes do not move with `-t`. Never: **assign into the chunk, never
  `clear()` it** -- clear destroys four `std::string` per record and the reader ends up in malloc
  instead of inflate; that was half the win.
- **Never: an estimate that nothing checks will be wrong, and this one was.** `choose_bucket_bits`
  assumed a gzipped FASTQ goes resident at **8x** on disk. Measured: **19x** -- a resident record is
  two heap `std::string` with allocator headers and rounded-up buckets, not the 180 bytes of
  payload. Guessing low is the expensive direction: it picks too few buckets and pass 2 holds
  `kBucketConcurrency` of them at once, which is the entire 1,479 MB. It is 20x now.
- **Note: the chunk is 8192 and a bigger one is 22% faster.** 64 k reads gives 2,510,241 reads/s,
  because `parallel_for` starts and joins its threads per call and 4 M reads at 8 k pays ~15,000
  thread creations. Not taken: 16 MB of resident chunk makes **pass 1** the memory peak on a finely
  partitioned shallow library, breaking "a finer partition costs less, not more" --
  `test_shallow_memory_is_still_bounded_by_the_bucket` is what caught it. The upgrade path is a
  persistent worker pool, not a bigger chunk.
- **`subsample` reports what it kept (2.0.0a4)**: median and deepest reads per kept barcode next to
  the mean, plus five example barcodes with their depths. Never: **in key order, not first-seen
  order** -- a 100-read barcode appears early ~100x more often than a singleton, so the head of a
  file samples the deep MIGs and nothing else. Same trap as subsampling reads, one level down.
- **The docs navigate (2.0.0a4).** Twenty pages of flat toctree put every long page title in the
  pydata header; they are seven sections now (Installation, Examples, Layouts, Commands, Downstream,
  Method, Reference) with landing pages saying what each page answers, and every command page has a
  subtitle. `docs/nextflow.rst` is new. Note: nextflow is not installed here, so the modules are
  reviewed against the nf-core spec, not verified by a run, and the docs say so.
- **Never: `containsKey`, not `?:`, for a nextflow boolean.** Groovy's elvis treats `false` as
  absent, so a per-sample `contig: false` against `params.migec_contig = true` silently meant its
  opposite -- the one direction a per-sample override exists to make possible.
- **`--pre-amp-error auto` fits the floor (2026-08-14).** X2's estimator run on assemble's own
  consensuses: molecules at >= 20 reads, the library's modal sequence with one vote per molecule,
  polymorphic positions and divergent templates excluded, and the residual IS the floor. 1e-4 ->
  1.20e-4 [7.94e-5, 1.75e-4]; 1e-5 -> 1.58e-5; 0 -> the bound 9.56e-6. Never: **zero observed
  mismatches is a BOUND, not a rate** -- take the interval's upper end, because a floor of zero is
  a Q-infinity. Never: **it refuses rather than guessing.** A diverse library has no monomorphic
  position, so "disagrees with the modal base" means "is a different molecule" and the answer would
  be the library's diversity; a shallow one never reaches the plateau. Both fall back to the named
  `rt` class and say why, in the report and in `assemble.pre_amp_error.tsv`. Note: it costs a
  SECOND assembly -- the consensus sequences do not depend on the floor, only the emitted quality
  does, so a probe assembly runs and is deleted. Rewriting the qualities of the gzipped output in
  place would have saved the pass and bought a new way to write a broken file.
- **R1/R2 overlap merge is in (2026-08-14).** `assemble --mate2 <R2>` on the FASTQ route,
  `--merge-mates` on the `.mig` route where checkout already stored both mates in the record
  (`seq2`/`qual2` -- assemble had simply been ignoring them). Mate 2 is rc'd once at bucket load
  and placed; overlapping mates give one consensus over the insert, non-overlapping give two
  contigs, and both routes agree byte for byte. Never: **the offset is a property of the MOLECULE,
  not of the pair.** Feeding the pair into `--contig`'s all-against-all `place_reads` cost **11x**
  the single-end path (119,820 record-pairs/s against 1,356,819); `assemble_pairs` votes once per
  group over <= 8 pairs, stopping at two agreeing votes, for 1,288,686 against 2,115,912. Never:
  the mates are matched by POSITION and a short file is refused -- pairing off the remainder
  attaches one molecule's mate to another's and the consensus looks fine.
- **refine's barcode table bounds itself (2026-08-14).** Past `table_budget_bytes` (1 GB, a
  `refine.run()` kwarg, never a CLI flag) the table range-partitions into `<out>/.refine_spill` and
  correction follows it, two passes with the key rotated. Nothing in the pipeline scales with the
  library any more. Never: **the table carries the EVIDENCE, not just the counts** --
  `BarcodeEvidence` as a side array is indexed against `entries()` and cannot survive a partition,
  and dropping it leaves the bucketed run on the count ratio alone, which reports nothing at 1-3
  reads/UMI. `UmiCounts::carry_evidence` partitions the barcode quality (summed, divided by the
  count at hand-over) and the payload draft (first read wins, so a stable sort is load-bearing)
  with the key. Never: **the two passes SCAN, they do not merge** -- a barcode can have a plausible
  parent on each side of the boundary and merging inside a pass takes the first rather than the
  best; both passes propose and ONE global apply decides, in the resident walk's order (child count
  ascending, key descending). Never: **clonality is sampled in KEY order, never by array index** --
  an index rule measures the array, so two partitionings drew different pairs and a borderline
  merge followed, which made the output depend on the memory budget; `PayloadReservoir` keeps every
  k-th barcode in key order and doubles k when full. With all three, resident and partitioned agree
  on every scalar and every output file byte for byte
  (`tests/synthetic/test_refine_bucketed.py`). Two frees came with it: the evidence pass folded
  into the table pass (three read passes -> two) and the rank curve is derived from the size
  spectrum rather than a sorted array of every molecule's count.
- **What the alignment position is worth, measured (2026-08-14).** `scripts/compare_grouping.py`
  runs UMI-tools `group` and fgbio `GroupReadsByUmi` end to end -- barcode into the read name /
  into `RX`, minimap2 onto the simulator's own `clones.fa`, then group -- and scores all three
  partitions against the truth with `compare_calib.adjusted_rand`. Never: **on ONE reference the
  position carries nothing and migec wins** -- ARI 0.9967 against 0.9864 (UMI-tools) and 0.9817
  (fgbio), with 0.65% of reads in clusters that MIX molecules against 3.0% and 3.9%. That is 4.6x
  and 6x fewer molecules destroyed, and destroying one is the error nothing downstream can detect.
  On 200 or 20,000 distinct references they win by 0.001 ARI, which is the barcode collision rate
  and nothing else. 8-48x faster throughout, including the alignment they cannot skip; depth
  (1.2/2.5/5/10 reads per molecule) changes neither ranking. `docs/grouping.rst`,
  `assets/grouping_tools.tsv`. Note: the simulator now writes `clones.fa` because a map-first tool
  cannot run at all without a reference. Note: fgbio needs a **JDK 17+** (JDK 11 gives
  `UnsupportedClassVersionError` from htsjdk, not a version message); UMI-tools needs
  `--no-build-isolation` and a `setuptools` on Python 3.12+.
- **MIGEC 1.2.9 head to head (2026-08-14).** `scripts/compare_migec_v1.py`, same dialect, same
  sheet, same library, both pipelines end to end. **9.4-11.9x the wall clock** against a 3x gate,
  3.5-3.7x less memory, and the molecule count 0.09% over truth against v1's 13.6% -- 99.8-99.99%
  of our consensuses are exactly a template against 93.9-95.1%. The over-count is barcode errors
  `--filter-collisions` missed: its rule is a count ratio, and a count ratio carries nothing below
  ~3 reads/UMI. Never: **`--min-count` goes to BOTH** -- v1 defaults to 5 and we default to 1, and
  v1 names its output `.t5.` for exactly that reason, so leaving each at its own default compares
  defaults and credits us with recovering molecules v1 was told to discard. `assets/migec_v1.tsv`,
  `docs/validation.rst`. Note: the jar is `gh release download 1.2.9`, and it runs on JDK 11.
- **Next, in order. `ROADMAP.md` has the same list with the reasoning; this is the short form.**
  1. **`2026-migec-benchmark`** and the published comparisons. This is what the version number is
     waiting on, not the code.
  2. **Run the callers themselves** against the ctDNA ground truth below. It no longer has to be
     built -- only the call sets are missing.
  3. **Splitting the fitted floor into RT and first-cycle PCR** — not an estimator problem. The
     same template through both chemistries is what it needs, which is data.
  4. **Bit-parallel matcher**, last and deliberately: the scan is not the bottleneck.
- **The ctDNA ground truth was FOUND, not built (2026-08-13).** `PRJNA788522` (72 runs, cfDNA
  reference material at certified 0/0.125/0.25/1% VAF x 5/20/80 ng x 3.3/10/30x, 3 reps) and
  `PRJNA507366` (28 runs, six polymerases plus 0.031%/0.0625% VAF) both kept a **real 12 nt inline
  UMI**, and `migec suggest` recovers it from base composition alone. Never: "no public ctDNA data
  has usable UMIs" is a statement about the two runs Maruzani picked, not about the archive --
  check read structure (`scripts/sra_fetch.py probe`) before believing it. Note: `PRJNA507366`'s
  design is in **`library_name`**, not `sample_alias`; every alias there is
  `SeraCare_Reference_Material`, so the alias alone makes a designed study look undesigned.
- **Never: the molecule total of a multiplex panel is not the count at a site.** A variant sits on
  one amplicon, so the evidence a caller gets is `molecules / amplicons`. The amplicon count is
  measured from consensus prefixes, and its share threshold must sit in the **gap** below the
  smallest real amplicon (5%, where the gap is 7.6% to 1.1%): at 1% the count moved with depth --
  5 on the deepest run, 10 on the shallowest, because a shallow consensus carries more payload
  error -- which deflates the per-amplicon count on exactly the runs where evidence is thinnest.
- **Note: `scripts/sra_fetch.py` is the fetch-on-demand path**, so `isalgo/umi_data` holds only CI
  fixtures and semi-internal data. Never: ENA's ready-made FASTQ looks like the lazy source and is
  33x slower -- **200 kB/s against NCBI S3's 6.7 MB/s** on 8 connections, measured on SRR17220895.
  S3 plus `fasterq-dump` is the default; `--prefer ena` exists for studies that deposited a third
  file the `.sra` folds away.
- **Note: Britanova et al aging (bulk TCR, shallow) lives on aldan3** and is the real dataset for the
  1-3 read regime. Not pulled yet. aldan3 compute goes through SLURM, never the frontend.
- **Note: the ctDNA titration runs on aldan3** (`~/migec-ctdna`, `migec_ctdna.sbatch`, partition
  `medium`, 16 cores / 32 GB). The job builds its own micromamba env because the cluster's system
  python is **3.8** and migec needs >= 3.10, and it installs from PyPI rather than from source,
  which doubles as a check that the published linux wheels work.
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
