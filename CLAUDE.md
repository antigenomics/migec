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

## Non-negotiables

These exist because the alternative was tried, or because a first-pass design got it wrong and the
correction is written up in `project/review-algorithms.md`.

- ⛔ **Every model-derived number is reported next to something that checks it.** The birthday
  collision rate has `scripts/collision_check.py` (model-free, from the sequences); the barcode
  error estimate has the Phred + polymerase prediction; the emitted quality has the measured RT
  floor. A formula that nothing tests is a formula that will be wrong silently — all three of
  these were, and the check is what found it.
- ⛔ **The distance-1 barcode-error estimator fails downward as the space fills** (0.92x of truth
  at 0.3% occupancy, 0.23x at 50%, 0.001x at 93%). `err_unreliable` is set past 5% neighbourhood
  occupancy. Never quote the estimate when it is set.
- ⛔ **Five commands: `checkout`, `suggest`, `refine`, `assemble`, plus `sort`/`subsample`.** A
  sixth command or a new CLI flag requires a *failing benchmark that the default cannot pass*.
  Constants live in headers next to the measurement that justifies them, so adding a flag means
  deleting a measurement.
- ⛔ **Never subsample reads to make an example.** Sort by UMI and take all reads of N UMIs
  (`blake2b(umi) mod 10000 < K`). Ten thousand random reads gives one read per UMI and is useless.
  Taking the *first* N distinct UMIs is also wrong — a UMI with 100 reads appears early ~100× more
  often than a singleton, so first-appearance order oversamples large MIGs and destroys the very
  size distribution the example is showing.
- ⛔ **Do not derate base quality for a low-coverage or suspected-error-child MIG.** If UMI *x* is
  genuinely an error child of *y*, then *all* of x's reads are clean reads of y's sequence: the
  consensus is right and only the molecule count is wrong. Emit the identity posterior as a tag;
  do not smear it into every base.
- ⛔ **Never emit a quality above the RT/first-cycle-PCR floor.** An error made before
  amplification is in every read and no consensus removes it. `p_floor = ε_RT + 2·ε_pol`, fitted
  from data as the intercept of `e_out(c) = p_floor + a/c`, never hard-coded.
- ⛔ **A null that permutes reads column by column is the wrong null.** A low-quality read carries
  a minor base at many positions at once, so column permutation hands every read an average error
  load and calls the resulting co-segregation significant. Preserve *both* margins of the
  reads × positions matrix (curveball). Measured: the nominal `p < 0.01` split threshold calls
  30.62% of MIGs, the both-margins null puts the 1% false-positive point 19x higher, at 8.68.
- ⛔ **A molecule is sample + cell + UMI, never the UMI alone.** UMIs repeat across cells and
  samples by design; grouping on the UMI merges two molecules and nothing downstream can tell.
  `assemble` sorts on `(cell, umi, src_index)` and partitions on the cell when there is one.
- ⛔ **Contig assembly means one molecule's fragments, nothing more.** Random priming leaves reads
  that tile a molecule under one barcode; `--contig` places them and emits one consensus per
  overlap component. Assembling a cell's full-length receptor, calling doublets and filtering
  contaminating chains are **arda's** job. ⚠ It also needs a barcode that is not saturated: two
  fragments of two different molecules on one barcode have no sequence in common, which is exactly
  what two fragments of one look like. Report `expected_molecules_per_group` and warn.
- ⛔ **1-3 reads per UMI is the normal case, not the exotic one.** Bulk repertoire profiling and
  shallow 3' GEX both look like this. Never threshold it away -- `--min-reads` defaults to 1 and
  the answer to a barcode error is correction, not a cut. Three things calibrated on a deep library
  stop applying and must be said rather than quoted: the split threshold is inert (needs ~30
  reads), the count-ratio error-child null has no dynamic range, and singleton filtering is
  unaffordable (79% of barcodes, against 56% on a deep library). Benchmarks use the shallow shape
  because distinct barcodes are what everything scales with.
- ⛔ **A collision rate is estimated with the U-statistic, never the plug-in.** `sum p_hat^2` is
  biased up by `(1 - sum p^2)/n`, and the bias GROWS as the distribution spreads -- so on k-mers it
  grows with k and reads as dependence accumulating with k. It manufactured a 1.007x excess out of
  independent data. Use `sum n_a(n_a-1)/(N(N-1))`.
- ⛔ **An `N` is not a fifth base.** Fold it to A, as `pack_barcode` does -- the packed key is what
  every stage groups on. Counting it as a letter let `m_j` fall to 0.2466, below the mathematical
  floor of 1/4, and printed an effective length of 9.01 nt for a 9 nt barcode. `m_j >= 1/4` is the
  check; use it.
- ⛔ **Range partition, never hash partition.** A hash sends a barcode and its 1-mismatch
  neighbours to uncorrelated buckets, which makes correction impossible to apply and splits the
  molecule permanently — and the halves each look like a well-formed MIG, so nothing detects it.
- ⛔ **No indels anywhere.** Substitutions only. Reinstate when there is a dataset to verify against.
- ⛔ **Nothing that scales with the input on the serial path, and no hash map keyed by barcode.**
  zlib level 6 does 7 MB/s on random DNA, so compression belongs on the workers (concatenated gzip
  members are a valid stream) at level 1. A `unordered_map<uint64,uint32>` costs ~48 B per distinct
  UMI against 22 for a sorted array, which at NovaSeq scale is 19 GB against 8.8. Both were
  measured; `tests/benchmark/` guards them.
- ⛔ **`-t` must never change the output.** Chunks are matched in parallel and written in input
  order. A demultiplexer whose output depends on its thread count produces results that cannot be
  compared between runs, and the failure is invisible.
- ⛔ **One output file per *sample id*, never per barcode-table row.** Rows sharing an id are a
  sample sequenced with more than one tag — the format's own idiom. A file per row opens the same
  path twice and interleaves two `FILE*` into it, which is not a valid gzip stream, and the summary
  still reports success. Anything else keyed by sample (UMI counters, summary rows) groups the same
  way.
- ⛔ **Nothing may throw out of a worker thread.** An escaping exception is `std::terminate`:
  SIGABRT, no message, no flush. Workers capture, the driver rethrows. Validate at the boundary
  where the error is attributable — a bad pattern is caught when the pattern is compiled, on the
  caller's thread, not when a read is packed on a worker.
- ⛔ **A reported clock covers the whole operation.** `wall_seconds` stopping at the parallel driver
  hid a serial stats pass worth 4/5 of the run and published a throughput 5× what a user sees. If a
  stage is serial, time it and report it separately — that is the number that says what to fix next.
- ⚠ **No `log2`/`exp`/`pow` in a per-base loop.** The score depends only on the reported Phred and
  the IUPAC set size, both small integers; it tabulates into 1.2 kB. The transcendental was 90% of
  checkout's runtime before it did.
- ⚠ **Strand normalisation happens in `checkout`.** The `.mig` flags describe what has *already*
  been applied. A group containing both orientations silently loses half its reads in consensus.
- ⚠ **Nominal Phred is not the error rate.** On 2-colour instruments there are ~4 distinct Q
  values. Use the measured calibration table from the `.mig` header wherever a likelihood is
  computed.
- ⚠ **Do not reproduce MIGEC v1's bugs.** Quality must be indexed at the *match offset*, not the
  read start; low-quality mismatches must actually be counted (v1's dangling `else` meant they
  never were).
- ✅ The version string is duplicated in `pyproject.toml`, `include/migec/version.hpp` and
  `python/migec/__init__.py`. `__init__.py` asserts the last two agree and CI prints all three.

## Release process

1. Bump the version in all three places above; add a `CHANGELOG.md` section.
2. Branch `release/<version>` → PR → CI and docs green.
3. `gh release create v<version>`.
4. `publish.yml` builds wheels and publishes via OIDC to the `pypi` environment. ⛔ Do not rename
   that workflow file or the environment: the PyPI trusted publisher is bound to both.

## Open loops

- M0 done. `migec checkout` works: patterns, trimming, header transfer, UMI statistics, count
  correction, **paired input with strand normalisation**, and **multi-core with byte-identical
  output at any `-t`** (1.18 M reads/s at 16 threads). Whitelists, dual-end barcodes and `.mig`
  bucket output are still open.
- ⚠ **The UMI counters are not partitioned yet.** ~22 B per distinct UMI is 8.8 GB at NovaSeq
  scale, held in one piece. The fix is the range partition (M2, with `.mig` bucket output), not a
  smaller struct. Until then checkout warns past 1 GB.
- ⚠ **The per-sample statistics are the serial tail.** Histogram, composition and `correct_umis`
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
- **X2 is done (2026-08-13): the RT/PCR floor is of order 1e-4, so cap emitted quality at ~Q38.**
  Measured 1.54e-4 [1.36e-4, 1.74e-4] on `SRR1763769`; the 1e-6 guess is dead. `--rt-error auto`
  still fits per dataset (the floor is a property of the enzyme and cycle count), but the
  **default is 1e-4**. ⚠ Still an upper bound — that library is 49.6% occupied on a 9 nt barcode
  and `checkout` calls it saturated, so collisions inflate it. `docs/quality_floor.rst`.
- **`migec suggest` is implemented** (2026-08-13), ahead of its M4 slot, because X2 needed it: it
  segments the per-cycle base composition into UMI (all four bases near 1/4), constant and payload
  runs and prints a paste-ready pattern. Recovered the 9 nt + `CAGTTTAACTTTTGGGCCAT` layout of
  SRR1763769 unaided. ⚠ It stops the pattern at the last *constant* run — composition alone cannot
  tell a UMI from diverse payload, only the anchor can.
- **The barcode-space and error-budget arithmetic is built in**, logged to
  `checkout.barcode_space.tsv` / `checkout.umi_quality.tsv`, warned on, documented in
  `docs/barcode_space.rst`, drawn in `notebooks/barcode_space.py` and tested in
  `tests/synthetic/test_barcode_space.py`.
- **X3 is done (2026-08-13).** Position independence holds (1.04x over 9 nt), so `Π_j m_j` stays
  and the 1.86x collision excess is the read threshold, not the barcode. 92% of distance-1 pairs
  are chance at 47.8% occupancy; the column-shuffle background puts barcode error at 3.4e-3 —
  within 1.7x of the Phred + polymerase prediction, where the analytic estimate is 2.6x below it,
  so **M3 takes the permuted background**. The split threshold is **8.68, not 2.00**.
  `docs/nulls.rst`, `scripts/permutation_nulls.py`, `tests/synthetic/test_nulls.py`.
- **M1 is done (2026-08-13): `migec assemble` works.** `(cell, umi)` grouping, range partition
  into `.mig` buckets with one bucket resident, the column log-likelihood posterior, the RT floor
  added (not compared) so nothing above ~Q38 is emitted, linkage sub-clustering at X3's 8.68, and
  `--contig` for random-primed reads (seed placement, union-find overlap components, never bridged
  across a gap). 531 k reads/s; 121 MB at 16 buckets against 203 MB at one; asserted in
  `tests/benchmark/test_assemble_speed.py` in a fresh process per configuration, because
  `peak_rss_bytes` is a process high-water mark.
- ⚠ **The writer buffer budget is split across buckets, not per bucket.** A fixed per-writer block
  made pass 1 cost grow with the bucket count, which is backwards — more buckets exist to use less
  memory, and cutting finer made peak RSS go *up* (238 MB at 16 buckets against 213 at one). It is
  now `clamp(32 MB / buckets, 256 kB, 4 MB)`.
- ⚠ **A tail quantile is not a constant until its Monte Carlo error is smaller than the digits
  quoted.** X3's split threshold read 9.91, then 9.61, then 11.66 on reruns of ~8,000
  randomisations. At 82,800 it is **8.68, bootstrap 95% CI [8.42, 9.14]**, and that interval is
  what the code uses and the docs quote.
- ⛔ **The count ratio is not evidence below ~3 reads/UMI**, and that is the common regime. Two
  evidence terms that survive at one read are now in `BarcodeEvidence`: the barcode's own **base
  quality** at the differing position, and **payload agreement** with the candidate parent, worth
  `log(1/clonality)` where the clonality is measured from random barcode pairs rather than assumed.
  Payload agreement lifts the count gates (which is what makes a singleton merge possible) and
  payload *disagreement* refuses a merge the counts would have made. Never tune
  `max_child_fraction` to paper over this — it is the wrong evidence, not the wrong threshold.
- ⛔ **Compare a rate to a rate.** The error likelihood was a zero-truncated Poisson, i.e. a
  probability conditioned on the child existing, weighed against `a_ind * p_size`, which is an
  expected *count*. The truncation divides out `(1 − e^−λ)` — exactly the term that says whether an
  error child should exist — so ZT-Poisson(1, λ) → 1 for every small λ and the error rate stopped
  mattering at precisely the coverage where nothing else was available. Untruncated.
- ⛔ **Err on precision, never on recall, when merging.** A wrong merge deletes a molecule and
  nothing downstream can tell; a missed correction only inflates the count and is recoverable. The
  current posterior destroys **no** real molecule at any depth measured (≥0.99, 1.000 at the
  extremes) and that is the property to protect.
- ⚠ **At ~1 read/UMI, 80% of barcode errors are unfixable in principle** — the parent barcode was
  never sequenced, so there is nothing to merge into. Always report recall against that reachable
  ceiling; against all children it understates by 5x and looks like a bug.
- Next: **M3** (`refine`), taking X3's permuted distance-1 background; then the M2 remainder
  (whitelists, dual-end barcodes, `.mig` bucket output from checkout, i7xi5).
- ⚠ **Britanova et al aging (bulk TCR, shallow) lives on aldan3** and is the real dataset for the
  1-3 read regime. Not pulled yet. aldan3 compute goes through SLURM, never the frontend.
- The archive is pushed: `legacy-v1` + tag `v1-final`, and master is the rewrite. Recovery point
  `~/backup/migec-local-mirror-2026-08-13.git`. ⚠ The canonical repo is **antigenomics/migec**;
  `mikessh/migec` is a redirect, and `gh` commands must use the former.
- ⚠ `gh repo edit --default-branch` and force-pushes may be blocked by the permission classifier;
  hand those to the user rather than working around them.
- `isalgo/umi_data` does not exist yet; nothing has been uploaded.
- ⚠ Do **not** add seqtk as a dependency. It is the right tool for generic FASTQ slicing in a
  benchmark harness, but it cannot subsample by whole UMI (it samples reads), which is the one
  thing we need and the thing that makes a UMI example fixture correct.
