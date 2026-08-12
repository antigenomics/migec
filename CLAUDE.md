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
- ⛔ **Range partition, never hash partition.** A hash sends a barcode and its 1-mismatch
  neighbours to uncorrelated buckets, which makes correction impossible to apply and splits the
  molecule permanently — and the halves each look like a well-formed MIG, so nothing detects it.
- ⛔ **No indels anywhere.** Substitutions only. Reinstate when there is a dataset to verify against.
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

- M0 done: `.mig` format frozen, FASTQ IO, simulator, CI, docs shell.
- Next, in order: **X1** (read-start dispersion on 10x — decides whether fragmented mode is
  mandatory), then M1 (`assemble` + the quality model, validated on a clonal control).
- Pushes to `origin` are pending a network window: `legacy-v1` and `v1-final` are additive and go
  first, then `master` with `--force-with-lease`. Recovery point:
  `~/backup/migec-local-mirror-2026-08-13.git`. Note the canonical repo is **antigenomics/migec**; `mikessh/migec` is a redirect.
- `isalgo/umi_data` does not exist yet; nothing has been uploaded.
