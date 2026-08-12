# MIGEC v2 — repo skeleton, build, CI, packaging, docs, notebooks, archive procedure

## 0. Verified state (read-only checks just run)

| Fact | Value |
|---|---|
| `migec` local HEAD | `e511a6e` on `master`, clean, in sync with `origin/master` |
| remote branches | `master`, `develop`, `feature/full-length-preproc` |
| tags | 28, latest `1.2.9`; releases 1.2.5–1.2.9 carry jar assets |
| branch protection / rulesets | none (`404 Branch not protected`, `rulesets: []`) — force-push is possible |
| `mikessh/migec` | public, not archived, default `master`, topics umi/antibody/airr/barcode/exome-sequencing/mutation/rna-seq/single-cell |
| `mikessh/mageri` | public, not archived, default `master`, last push 2017-05-08, 5 releases |
| PyPI `migec` | **FREE** (404) |
| PyPI `mageri` | free (404) |
| seqtree wheel | ships `python/seqtree` only — **no C++ headers in the wheel**, no `install()/export()` |

**BLOCKER to surface to the user before any of this:** `/Users/mikesh/vcs/code/migec/LICENSE` is a custom **OOO «MiLaboratory» proprietary agreement** (`Copyright (c) 2014-2024, OOO «MiLaboratory»`, 76 lines), not GPLv3. GitHub reports `licenseInfo.key = "other"`. The rewrite must not silently inherit it, and seqtree/arda are GPL-3.0-or-later. Decide the new license explicitly (recommend GPL-3.0-or-later to match seqtree/arda, since we link seqtree). Keep the old file on `legacy-v1` only.

---

## 1. Archive procedure

### 1.1 Which shape

Two readings of "archived in orphan branch, start from new master":

| Option | `legacy-v1` | `master` | Verdict |
|---|---|---|---|
| A (recommended) | **normal branch** at `e511a6e`, full 10-year history, additive push | **orphan**, new empty tree, one force-push | Recommended |
| B | orphan branch = one squashed snapshot commit | orphan, new | **Reject** — destroys the history that is the whole point of archiving |

Option A already satisfies the intent: the *new master* is the orphan; the archive keeps its history. All old objects stay reachable through `refs/heads/legacy-v1` and the tags, so the single force-push cannot lose anything.

### 1.2 Exact sequence — `mikessh/migec`

```bash
# ---- STEP 0: verify + off-site recovery point (nothing written to origin) ----
cd /Users/mikesh/vcs/code/migec
git fetch --all --tags --prune
git status --porcelain                 # MUST print nothing
git rev-parse master                   # record this sha; expect e511a6e...
git clone --mirror https://github.com/mikessh/migec.git \
    ~/backup/migec-mirror-$(date +%F).git        # full recovery point, incl. all refs

# ---- STEP 1: park the history on a branch. ADDITIVE, no force ----
git branch legacy-v1 master
git push origin legacy-v1              # creates a new ref; cannot destroy anything

# ---- STEP 2: annotated tag on the final v1 commit. ADDITIVE, no force ----
#   name deliberately != the branch name, so `git checkout legacy-v1` is never ambiguous
git tag -a v1-final -m "Final Groovy MIGEC (1.2.9) before the C++20/Python rewrite" master
git push origin v1-final

# ---- STEP 3: point the default branch at the archive BEFORE master changes ----
gh repo edit mikessh/migec --default-branch legacy-v1
gh api repos/mikessh/migec --jq .default_branch     # verify == legacy-v1

# ---- STEP 4: build the new master locally (still nothing pushed) ----
git checkout --orphan master-new       # no parent commit, index still holds old files
git rm -rf .                           # empty the index AND the worktree
#   ... lay down the new tree from §2 here ...
git add -A
git commit -m "MIGEC v2: C++20 core + Python CLI (complete rewrite)"
git log --oneline                      # MUST be exactly 1 commit
git ls-files | head                    # MUST show only new files

# ---- STEP 5: the ONLY destructive step ----
git branch -M master-new master        # -M force-moves the local ref
git push --force-with-lease origin master
#   --force-with-lease (not --force): aborts if origin/master moved since the STEP 0 fetch

# ---- STEP 6: metadata ----
gh repo edit mikessh/migec --default-branch master
gh repo edit mikessh/migec \
  --description "UMI barcode extraction, correction and consensus assembly — C++20 core, Python CLI" \
  --homepage "https://mikessh.github.io/migec/"
gh repo edit mikessh/migec --add-topic consensus-sequence --add-topic cpp --add-topic pybind11
#   existing topics (umi, single-cell, airr, barcode, rna-seq, ...) are kept; --add-topic is additive
```

Push-safety summary:

| Step | Ref | Force? | Recoverable from |
|---|---|---|---|
| 1 | `origin/legacy-v1` | no | — (creation) |
| 2 | `origin/v1-final` | no | — (creation) |
| 5 | `origin/master` | **`--force-with-lease`** | `origin/legacy-v1`, tag `v1-final`, mirror clone, GitHub reflog (90 d) |

Leave `origin/develop` and `origin/feature/full-length-preproc` untouched — free extra recovery points. Existing releases attach to tags, not branches, so the 1.2.x jars survive; JitPack builds from tags and is unaffected. Known breakage: permalinks/PR compare-links into old `master`.

`legacy-v1` should get one commit on top (optional, additive) prepending to `README.md`:
> **This branch is the archived Groovy MIGEC 1.2.9.** Development continues as a C++/Python rewrite on `master`. Java users: `git checkout v1-final`, or the jars on the 1.2.9 release.

### 1.3 Exact sequence — `mikessh/mageri`

```bash
git clone https://github.com/mikessh/mageri.git /tmp/mageri && cd /tmp/mageri
git fetch --all --tags --prune
git clone --mirror https://github.com/mikessh/mageri.git ~/backup/mageri-mirror-$(date +%F).git

git branch legacy-v1 master && git push origin legacy-v1        # additive
git tag -a v1-final -m "Final MAGERI (1.1.1); superseded by MIGEC v2" master
git push origin v1-final                                        # additive
gh repo edit mikessh/mageri --default-branch legacy-v1

git checkout --orphan master-new
git rm -rf .
#   write README.md (stub below) + LICENSE (copied from legacy-v1)
git add -A && git commit -m "Archive MAGERI; redirect to MIGEC v2"
git branch -M master-new master
git push --force-with-lease origin master
gh repo edit mikessh/mageri --default-branch master
gh repo edit mikessh/mageri --description "ARCHIVED — superseded by https://github.com/mikessh/migec"
gh repo archive mikessh/mageri --yes        # LAST; makes the repo read-only (reversible via gh repo unarchive)
```

`README.md` stub (whole file):

```markdown
# MAGERI — archived

MAGERI is no longer developed. Its consensus assembly and UMI quality model live on in
**[MIGEC v2](https://github.com/mikessh/migec)** (C++20 core + Python CLI).

* Alignment and variant calling are **not** carried over; MIGEC v2 emits consensus FASTQ with
  sample/cell/UMI in the read header, for `minimap2`, `bwa-meme`, or
  [arda](https://github.com/antigenomics/arda).
* Original source: branch [`legacy-v1`](../../tree/legacy-v1), tag `v1-final`, releases 0.1–1.1.1.
* Paper: Shugay et al. (2017) PLoS Comput Biol, PMID 28475621 —
  [mageri-paper](https://github.com/mikessh/mageri-paper).
```

**Do not run any of the above now.** Nothing here has been executed.

---

## 2. Names + file tree

| Thing | Value | Note |
|---|---|---|
| PyPI distribution | `migec` | verified free — no `arda-mapper`-style rename needed |
| Python import | `migec` | |
| CLI | `migec` | `[project.scripts] migec = "migec.cli:app"` |
| C++ namespace | `migec` | nested `migec::detail` for internals |
| Static lib target | `migec_core` | |
| Extension module | `migec._core` (CMake target `_core`, `install(... DESTINATION migec)`) | |
| CMake option prefix | `MIGEC_` | |

Python package sits under `python/migec/` (seqtree convention), because `src/` is the C++ sources.

```
migec/
├── CMakeLists.txt
├── pyproject.toml
├── setup.sh                       # uv env + --no-build-isolation + asserts `import migec._core`
├── LICENSE                        # NEW license, decided by user (see §0 blocker)
├── README.md
├── CLAUDE.md  ROADMAP.md  SOURCES.md  CHANGELOG.md
├── .gitignore  .gitattributes
├── include/migec/
│   ├── types.hpp          # ReadView, ReadPair, Umi (2-bit packed), BarcodeHit, phred LUTs
│   ├── fastq.hpp          # FastqReader / FastqWriter (zlib), zero-copy record views
│   ├── barcode.hpp        # BarcodePattern (grammar compiler), BarcodeMatcher, IUPAC LUT
│   ├── sample_sheet.hpp   # SampleSheet, SampleSpec, whitelist loading
│   ├── checkout.hpp       # CheckoutEngine, CheckoutStats
│   ├── umi_stats.hpp      # CoverageHistogram, UmiPwm, error-rate + birthday-collision model
│   ├── correct.hpp        # UmiCorrector (neighbour search), CorrectionTable
│   ├── extsort.hpp        # ExternalSorter, RecordKey{sample,cell,umi}
│   ├── consensus.hpp      # ConsensusAssembler, MigGroup, cqs(), doublet/collision split
│   └── version.hpp
├── src/                   # fastq.cpp barcode.cpp sample_sheet.cpp checkout.cpp
│                          # umi_stats.cpp correct.cpp extsort.cpp consensus.cpp _bindings.cpp
├── python/migec/
│   ├── __init__.py        # __version__, re-exports from _core
│   ├── cli.py             # typer app: checkout, stats, correct, refine, assemble,
│   │                      #            suggest, subsample, simulate, report
│   ├── checkout.py  stats.py  correct.py  refine.py  assemble.py
│   ├── suggest.py         # infer UMI/primer placement from data
│   ├── subsample.py       # UMI-preserving subsample (§8)
│   ├── simulate.py        # ground-truth read simulator (§7)
│   ├── report.py  plots.py  schema.py     # polars schemas for every output table
│   ├── _core.pyi  py.typed
├── tests/
│   ├── cpp/{doctest.h,test_barcode.cpp,test_fastq.cpp,test_consensus.cpp,
│   │        test_extsort.cpp,test_umi_stats.cpp}
│   ├── conftest.py
│   ├── unit/  synthetic/  realworld/  benchmark/
│   └── data/              # committed: <1 MB, 200 read pairs, 3 samples
├── docs/                  # flat, arda-style (§6)
├── notebooks/             # marimo (§8)
├── bench/bench_migec.cpp  bench/tables/  bench/plots/
├── examples/              # barcodes.txt, sample_sheet.tsv, run.sh
├── scripts/fetch_testdata.sh
└── .github/workflows/{ci.yml,docs.yml,publish.yml,testpypi.yml}
```

---

## 3. CMakeLists.txt

```cmake
cmake_minimum_required(VERSION 3.20)
project(migec LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
if(NOT MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()  # MSVC Release already uses /O2; -O3 is not a valid MSVC flag

option(MIGEC_TESTS  "Build C++ tests"        OFF)
option(MIGEC_BENCH  "Build C++ benchmarks"   OFF)
option(MIGEC_PYTHON "Build the pybind11 module" OFF)

find_package(Threads REQUIRED)
find_package(ZLIB REQUIRED)   # gzipped FASTQ; present in manylinux_2_28 and the macOS SDK

add_library(migec_core STATIC
  src/fastq.cpp
  src/barcode.cpp
  src/sample_sheet.cpp
  src/checkout.cpp
  src/umi_stats.cpp
  src/correct.cpp
  src/extsort.cpp
  src/consensus.cpp
)
target_include_directories(migec_core PUBLIC include PRIVATE src)
target_link_libraries(migec_core PUBLIC Threads::Threads ZLIB::ZLIB)

if(MIGEC_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_core src/_bindings.cpp)
  target_link_libraries(_core PRIVATE migec_core)
  install(TARGETS _core DESTINATION migec)
endif()

if(MIGEC_TESTS)
  enable_testing()
  add_executable(migec_tests
    tests/cpp/test_barcode.cpp
    tests/cpp/test_fastq.cpp
    tests/cpp/test_umi_stats.cpp
    tests/cpp/test_extsort.cpp
    tests/cpp/test_consensus.cpp
  )
  target_include_directories(migec_tests PRIVATE tests/cpp src)
  target_link_libraries(migec_tests PRIVATE migec_core)
  add_test(NAME migec_tests COMMAND migec_tests)
endif()

if(MIGEC_BENCH)
  add_executable(migec_bench bench/bench_migec.cpp)
  target_link_libraries(migec_bench PRIVATE migec_core)
endif()
```

### 3.1 seqtree — evaluated, and my recommendation is to NOT link it in C++ for v1

| Option | Verdict |
|---|---|
| link the pip-installed seqtree | **Impossible.** `wheel.packages = ["python/seqtree"]` — the wheel contains no `include/seqtree/*.hpp` and no static lib. |
| git submodule `extern/seqtree` | works, but `scikit-build-core`'s sdist is git-file-list driven; submodule contents in the sdist need care, and every contributor needs `--recurse-submodules`. |
| `FetchContent` pinned tag | works offline-hostile: every cibuildwheel leg clones seqtree; a network blip fails a release. |
| **seqtree as a plain Python runtime dep (`seqtree>=0.6`), like arda does** | **Recommended for v1.** |

Rationale for the recommendation: UMI/barcode correction runs on the **distinct barcode list**, not per read — after the external sort there are ~10⁵–10⁷ distinct UMIs and a *single* `Index.search_batch(queries, params, threads)` call, which already releases the GIL and threads internally. One batch crossing of the Python boundary per stage is unmeasurable against a multi-GB gzip decode. Requiring a C++ linkage buys nothing and costs a fragile build. seqtree ships wheels on exactly the platforms we ship, with no runtime deps of its own, so it is a core dependency (not an extra) — the arda lesson: an optional heavy-stage dep dies after 45 min of work on a bare `ModuleNotFoundError`.

If profiling later shows the boundary matters, this is the whole change (append to CMakeLists):

```cmake
option(MIGEC_WITH_SEQTREE "Link the seqtree C++ core" OFF)
if(MIGEC_WITH_SEQTREE)
  include(FetchContent)
  FetchContent_Declare(seqtree
    GIT_REPOSITORY https://github.com/antigenomics/seqtree.git
    GIT_TAG v0.6.1 GIT_SHALLOW TRUE)
  FetchContent_MakeAvailable(seqtree)          # dev override: -DFETCHCONTENT_SOURCE_DIR_SEQTREE=/Users/mikesh/vcs/code/seqtree
  target_link_libraries(migec_core PRIVATE seqtree_core)
  target_compile_definitions(migec_core PRIVATE MIGEC_HAVE_SEQTREE)
endif()
```

Upstream seqtree patch that would enable `find_package(seqtree CONFIG)` (5 lines, for completeness — **not worth doing for v1**, do it when a second C++ consumer exists):

```cmake
include(GNUInstallDirs)
target_include_directories(seqtree_core PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include> $<INSTALL_INTERFACE:include>)
install(TARGETS seqtree_core EXPORT seqtreeTargets ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(DIRECTORY include/seqtree DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(EXPORT seqtreeTargets NAMESPACE seqtree:: FILE seqtreeConfig.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/seqtree)
```

### 3.2 Compression deps

- **zlib — required**, `find_package(ZLIB REQUIRED)`. Present in manylinux_2_28 and the macOS SDK. No vendoring.
- **libdeflate — deferred.** It is a real 2–3× on gzip decode, which is the dominant cost of `checkout`. But it means a second code path and a vendored dep in every wheel. Keep the reader behind one `migec::GzipReader` in `fastq.hpp` so it is a one-file swap later. ROADMAP item.
- **zstd — deferred.** External-sort spill files are transient; v1 writes them **uncompressed** with `--tmp-dir` and a fixed-size record layout (which is what makes the sort fast). Revisit only if spill I/O measures as the bottleneck at 100M+ reads.
- **doctest** — vendored `tests/cpp/doctest.h`, exactly as seqtree does. No dep.

---

## 4. pyproject.toml

```toml
[build-system]
# pybind11 is UPPER-BOUNDED: PYBIND11_MODULE switched to multi-phase init in 3.0.0, so an
# unbounded requirement builds wheels against whatever PyPI serves that day. Pinned for
# REPRODUCIBILITY (same rationale as arda).
requires = ["scikit-build-core>=0.10", "pybind11>=3.0.2,<4"]
build-backend = "scikit_build_core.build"

[project]
name = "migec"
version = "2.0.0"
description = "UMI barcode extraction, correction and consensus assembly for UMI-tagged sequencing"
readme = "README.md"
requires-python = ">=3.10"
license = { text = "GPL-3.0-or-later" }        # PENDING the §0 license decision
authors = [{ name = "Mikhail Shugay", email = "mikhail.shugay@gmail.com" }]
keywords = ["umi", "barcode", "consensus", "single-cell", "airr", "rep-seq", "duplex-sequencing"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Programming Language :: Python :: 3",
    "Programming Language :: C++",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dependencies = [
    "polars>=1.0",
    "typer>=0.12",
    # CORE, not an extra. Barcode/UMI error correction is one batched `Index.search_batch`
    # per stage (GIL released, internally threaded), so it is not on a hot Python path -- but
    # every non-trivial run needs it. As an extra it would fail AFTER the expensive pass, past
    # `migec --version`, exactly the failure mode arda had to fix.
    "seqtree>=0.6",
]

[project.optional-dependencies]
docs = ["sphinx", "pydata-sphinx-theme"]
test = ["pytest", "pytest-cov"]
notebooks = ["marimo", "altair"]
dev = ["ruff", "pytest", "pybind11", "marimo"]

[project.scripts]
migec = "migec.cli:app"

[project.urls]
Homepage = "https://github.com/mikessh/migec"
Repository = "https://github.com/mikessh/migec"
Documentation = "https://mikessh.github.io/migec/"
Changelog = "https://github.com/mikessh/migec/blob/master/CHANGELOG.md"

[tool.scikit-build]
wheel.packages = ["python/migec"]
cmake.version = ">=3.20"
build-dir = "build/{wheel_tag}"
editable.rebuild = true
editable.verbose = false
sdist.exclude = ["/bench", "/docs", "/.github", "/notebooks", "/build", "/tests/data"]

[tool.scikit-build.cmake.define]
MIGEC_PYTHON = "ON"

[tool.ruff]
line-length = 100
target-version = "py310"

# Rules SELECTED EXPLICITLY, never inherited: ruff's defaults widen between releases, and CI
# installs a pinned ruff for the same reason. E4/E7/E9 + F = pyflakes/pycodestyle core: real
# errors, no style opinions. `I` deliberately not selected.
[tool.ruff.lint]
select = ["E4", "E7", "E9", "F"]

[tool.pytest.ini_options]
testpaths = ["tests"]
addopts = "-ra"
```

---

## 5. CI — four workflows

### `.github/workflows/ci.yml` — three jobs

```yaml
name: CI
on:
  push: { branches: [master, dev] }
  pull_request:
  workflow_dispatch:
concurrency: { group: ci-${{ github.ref }}, cancel-in-progress: true }

jobs:
  lint:                                  # runs first: seconds, not after the C++ build
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: { python-version: "3.12" }
      - uses: astral-sh/setup-uv@v6
      # PINNED. An unpinned linter in a required gate fails on someone else's release note.
      - run: uv pip install --system 'ruff==0.15.9' && ruff check python/

  cpp:
    strategy: { fail-fast: false, matrix: { os: [ubuntu-latest, macos-latest] } }
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - run: cmake -S . -B build -DMIGEC_TESTS=ON -DCMAKE_BUILD_TYPE=Release
      - run: cmake --build build --parallel
      - run: ctest --test-dir build --output-on-failure

  python:
    strategy:
      fail-fast: false
      matrix: { os: [ubuntu-latest, macos-latest], python: ["3.10", "3.12"] }
    runs-on: ${{ matrix.os }}
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with: { python-version: ${{ matrix.python }} }
      - uses: astral-sh/setup-uv@v6
      # build deps installed into the runtime env so editable.rebuild=true can find pybind11
      - run: |
          uv pip install --system scikit-build-core pybind11 ninja
          uv pip install --system -e ".[test]" --no-build-isolation
      # assert the EXTENSION imports, not just the package
      - run: python -c "import migec, migec._core as c; print(migec.__version__, c.__doc__)"
      - run: pytest tests/unit tests/synthetic -q     # synthetic == the simulator; no downloads
```

`tests/realworld` is **not** in CI (needs fetched SRA/10x data); `tests/benchmark` is gated by `RUN_BENCHMARK=1` and never runs on PRs.

### `.github/workflows/docs.yml`

Identical shape to arda's: push to `master` + `workflow_dispatch`, `permissions: {contents: read, pages: write, id-token: write}`, `concurrency: pages`. Build step is the gate:

```yaml
      - run: uv pip install --system -e . && uv pip install --system -r docs/requirements.txt
      - run: sphinx-build -W --keep-going -b html docs docs/_build/html
      - uses: actions/upload-pages-artifact@v3
        with: { path: docs/_build/html }
  deploy:
    needs: build
    environment: { name: github-pages, url: "${{ steps.deployment.outputs.page_url }}" }
    steps: [{ id: deployment, uses: actions/deploy-pages@v4 }]
```

`docs/conf.py` sets `autodoc_mock_imports = ["migec._core"]`; `polars` and `typer` are installed for real (imported at module level by the CLI, not mockable).

### `.github/workflows/publish.yml`

arda's five-job shape (`build-sdist` → `build-wheels` → `test-wheels` → `publish`), minus the reference-asset job, plus:

```yaml
        env:
          CIBW_BUILD: "cp310-* cp311-* cp312-* cp313-*"
          CIBW_SKIP: "*-musllinux_* *-manylinux_i686 *-win32 *-win_amd64"
          CIBW_MANYLINUX_X86_64_IMAGE: manylinux_2_28     # zlib-devel present
          CIBW_ARCHS_MACOS: "x86_64 arm64"
          CIBW_TEST_SKIP: "*"
```

Publish job: `environment: pypi`, `id-token: write`, `pypa/gh-action-pypi-publish@release/v1`, preceded verbatim by arda's version-vs-tag assertion:

```python
version = tomllib.loads(Path("pyproject.toml").read_text())["project"]["version"]
tag = "${{ github.event.release.tag_name }}"; norm = tag[1:] if tag.startswith("v") else tag
if version != norm: raise SystemExit(f"Version mismatch: pyproject.toml={version}, tag={tag}")
```

`test-wheels` smoke matrix mirrors arda's (glob for exactly one wheel per cp tag, `pip install`, then) — but assert more than the import, since `migec --version` succeeding without a working core is exactly the arda failure mode:

```bash
python -c "import migec, migec._core; print(migec.__version__)"
migec simulate --n-molecules 20 --coverage 5 --out sim/ && migec assemble sim/R1.fastq.gz --min-count 2 -o out/
```
No `pip install --upgrade pip` step (arda's Windows-blocked-release incident; moot here anyway).

### `.github/workflows/testpypi.yml`

arda's verbatim: `workflow_dispatch` only, same CIBW env as publish.yml, `repository-url: https://test.pypi.org/legacy/`, `CIBW_TEST_COMMAND: 'python -c "import migec; print(migec.__version__)"'`.

### Differences from arda, and my recommendations

| Question | Recommendation |
|---|---|
| Windows wheels? | **Drop them.** MIGEC v2 is a cluster/pipeline tool driven by gz FASTQ and an on-disk external sort; zlib on Windows needs vcpkg or vendoring, and arda's Windows legs have already blocked one release. `CIBW_SKIP` includes `*-win_amd64`; drop `windows-latest` from all matrices. If someone asks, WSL. |
| manylinux | `manylinux_2_28` explicitly — ships `zlib-devel`, has a modern GCC for C++20 (`manylinux2014`'s GCC 10 is borderline). |
| macOS | build both `x86_64` and `arm64`; `MACOSX_DEPLOYMENT_TARGET=10.15` for C++20. |
| musllinux | skipped (no C++20 toolchain payoff, no user demand). |
| C++ job on Windows | no. |

---

## 6. Docs

Theme: `pydata-sphinx-theme` (arda/seqtree). Flat `docs/` dir, `docs/requirements.txt`, `docs/_static/images/splash.png` carried over from the old `doc/_static/images/`. Gate: `sphinx-build -W --keep-going -b html docs docs/_build/html` in CI **and** in the `docs` skill locally.

Old → new page map:

| Old `doc/*.rst` | New `docs/*.rst` | Change |
|---|---|---|
| `index.rst` | `index.rst` | splash, 5-line quickstart, card grid |
| `install.rst` | `install.rst` | `pip install migec` / `uv`, wheels, build-from-source, no more Java/JitPack |
| `pipeline.rst` | `pipeline.rst` | end-to-end; ends at consensus FASTQ + hand-off |
| `checkout.rst` | `checkout.rst` | barcode grammar, positional vs pattern, dual-end, whitelists, sample sheets |
| `histogram.rst` | `stats.rst` | coverage histogram, over-seq test, threshold formula, UMI PWM/uniformity |
| — | `correct.rst` | **new**: neighbour search, quality+birthday model, keep-orphan-low-coverage rule |
| — | `refine.rst` | **new**: external sort, header transfer, arda/minimap2/bwa-meme-consumable output |
| `assemble.rst` | `assemble.rst` | consensus, multiple groups per UMI, doublets, RT-error quality cap |
| `cdrblast.rst` | — | **dropped** (out of scope) |
| `cdrfinal.rst` | — | **dropped** (out of scope) |
| `post.rst` | `downstream.rst` | merged: hand-off to arda / vdjtools / minimap2 / Cell Ranger |
| `logs.rst` | `logs.rst` | every output table, column by column |
| `report.rst` | `report.rst` | `migec report` (HTML) + the marimo notebooks |
| — | `formats.rst` | **new**: `barcodes.txt`, sample sheet, manifest, correction table, consensus header |
| — | `suggest.rst` | **new**: inferring UMI/primer placement from data |
| — | `simulate.rst` | **new**: the ground-truth simulator + its truth files |
| — | `api.rst` | **new**: autodoc of `migec.*` with `migec._core` mocked |
| — | `migration.rst` | **new**: v1 Groovy → v2 command table; where `legacy-v1` lives |
| — | `faq.rst` | **new**; leads with the UMI-subsampling rule (§8) |

---

## 7. Tests

Five tiers, arda's shape, gated by **constants in `conftest.py`, not markers** (except the benchmark module-level `pytestmark`).

| Tier | Path | Runs in CI | Content |
|---|---|---|---|
| C++ unit | `tests/cpp/` (doctest, ctest) | yes | pattern compiler + IUPAC/quality acceptance, FASTQ round-trip incl. truncated/corrupt input, PWM & CQS arithmetic, external sort ordering + spill boundaries, consensus offset search |
| Python unit | `tests/unit/` | yes | schemas, CLI arg parsing, table columns, `_core` binding types, pure-python helpers |
| Synthetic | `tests/synthetic/` | yes | everything from the **simulator** — correctness against known truth |
| Real-world | `tests/realworld/` | no | committed 200-pair fixture always; fetched SRA/10x/MIGEC-exp1 sets when present |
| Benchmark | `tests/benchmark/` | no | throughput + peak RSS regression vs a committed baseline |

`tests/conftest.py`:

```python
import os, pytest
from pathlib import Path

DATA = Path(__file__).parent / "data"
BIG  = Path(os.environ.get("MIGEC_TEST_DATA", DATA / "big"))

def _core_available() -> bool:
    try:
        import migec._core  # noqa: F401
        return True
    except Exception:
        return False

HAS_CORE      = _core_available()
HAS_SEQTREE   = __import__("importlib.util", fromlist=["util"]).find_spec("seqtree") is not None
HAS_BIGDATA   = BIG.is_dir() and any(BIG.glob("*.fastq.gz"))
RUN_BENCHMARK = os.environ.get("RUN_BENCHMARK") == "1"
RUN_SLOW      = os.environ.get("RUN_SLOW") == "1"

requires_core     = pytest.mark.skipif(not HAS_CORE,    reason="migec._core not built (run ./setup.sh)")
requires_seqtree  = pytest.mark.skipif(not HAS_SEQTREE, reason="seqtree not installed")
requires_bigdata  = pytest.mark.skipif(not HAS_BIGDATA, reason="set $MIGEC_TEST_DATA (scripts/fetch_testdata.sh)")
requires_slow     = pytest.mark.skipif(not RUN_SLOW,    reason="set RUN_SLOW=1")
# tests/benchmark/*.py carry: pytestmark = pytest.mark.skipif(not RUN_BENCHMARK, reason="set RUN_BENCHMARK=1")
```

### 7.1 The simulator — the single most valuable test asset

Ships in the package (`python/migec/simulate.py`), not in `tests/`, so users can reproduce every claim. CLI: `migec simulate`.

```python
@dataclass(frozen=True)
class SimConfig:
    # molecules
    n_molecules: int = 1000
    template: str = "random:300"      # or a FASTA path; molecules are drawn from its records
    n_templates: int = 5              # distinct source sequences when template is random
    # barcodes
    umi_len: int = 12
    master_pattern: str = "NNNNNNNNNNNNagatcggaagagc"   # v2 grammar, same as checkout
    slave_pattern: str | None = None                     # dual-end barcoding
    sample_barcodes: dict[str, str] | None = None        # sample_id -> barcode
    cell_barcodes: list[str] | None = None               # scRNA-seq CB whitelist
    # coverage
    coverage_mean: float = 20.0
    coverage_dist: str = "lognormal"  # "poisson" | "geometric" | "lognormal"
    coverage_cv: float = 1.0
    # errors
    pcr_cycles: int = 25
    pcr_efficiency: float = 0.9
    pcr_error_rate: float = 1e-5      # per base per cycle  (also caps achievable consensus Q)
    rt_error_rate: float = 1e-4       # per base, applied ONCE before amplification
    seq_error_rate: float = 1e-3      # per base
    umi_seq_error_rate: float | None = None     # defaults to seq_error_rate
    # pathologies the pipeline must survive
    umi_collision_rate: float = 0.0   # fraction of UMIs deliberately shared by 2 molecules
    chimera_rate: float = 0.0
    doublet_rate: float = 0.0         # fraction of CBs carrying two "cells"
    # reads
    read_len: int = 150
    paired: bool = True
    fragment_len: int = 350           # < 2*read_len ⇒ overlapping mates (contig assembly path)
    seed: int = 42

def simulate(cfg: SimConfig, outdir: Path) -> SimTruth: ...
```

Emitted into `outdir` (**exact file list**):

| File | Columns / content |
|---|---|
| `R1.fastq.gz`, `R2.fastq.gz` | the reads; header `@sim:<read_id>:<molecule_id>:<umi_true>:<n_umi_err>:<n_seq_err>` so truth survives tools that drop TSVs |
| `truth_molecules.tsv` | `molecule_id`, `sample_id`, `cell_barcode`, `umi_true`, `template_id`, `n_reads`, `n_rt_mutations`, `n_pcr_mutations`, `is_chimera` |
| `truth_reads.tsv` | `read_id`, `molecule_id`, `umi_true`, `umi_observed`, `n_umi_errors`, `n_seq_errors`, `strand`, `offset` |
| `truth_umis.tsv` | `umi_true`, `n_molecules`, `n_reads`, `is_collision` (= `n_molecules > 1`), `has_error_child` |
| `truth_consensus.fasta` | `>` `molecule_id` → the sequence a perfect assembler MUST emit (post-RT/PCR, pre-sequencing-error) |
| `truth_cells.tsv` | `cell_barcode`, `n_molecules`, `is_doublet`, `source_cells` (only when `cell_barcodes` set) |
| `config.json` | full `SimConfig` + seed + `migec.__version__` |
| `MANIFEST.tsv` | `file`, `sha256`, `n_bytes` |

Assertions the synthetic tier derives directly from these:

- UMI correction precision/recall vs `truth_reads.umi_true`; the keep-orphan rule verified as: at `coverage ≤ 5` and no parent, retained fraction = 1.0.
- Collision detection recall/precision vs `truth_umis.is_collision`, swept over `umi_collision_rate ∈ {0, 0.01, 0.05}`.
- Consensus per-base error rate vs `truth_consensus.fasta`; must be `≤ pcr_error_rate * pcr_cycles` and **must not** claim a quality above `-10*log10(rt_error_rate + pcr_error_rate*pcr_cycles)`.
- Round-trip: `checkout` recovers `umi_true` for every read with `n_umi_errors == 0`, exactly.
- Determinism: same `seed` ⇒ identical `MANIFEST.tsv` sha256s (this is also the CI regression oracle).

---

## 8. Marimo notebooks

`notebooks/*.py` (marimo files are plain Python, so they diff and lint), `[notebooks]` extra, listed in `docs/report.rst`, **not** in `sdist.exclude`… actually excluded from the sdist to keep it lean; linked from the docs.

| Notebook | Shows | Data |
|---|---|---|
| `01_checkout_and_umi_qc.py` | pattern grammar interactively; per-sample yield; undef-m/undef-s; UMI base-composition PWM and its deviation from uniform | `examples/` bundle |
| `02_coverage_and_threshold.py` | reads-per-UMI histogram (log2 bins), the over-sequencing test, the auto threshold, slider over `--min-count` | `examples/` bundle |
| `03_correction_and_collisions.py` | neighbour graph, parent/child ratio, birthday-collision probability curve vs UMI length, the correction table, and the low-coverage-orphan retention rule | `examples/` bundle |
| `04_consensus_accuracy.py` | CQS vs coverage; measured error rate vs `truth_consensus`; the quality cap | **`migec simulate` — no download at all** |
| `05_singlecell_10x.py` | CB+UMI extraction, barcode-rank knee plot, doublet calls | 10x subsample |

### 8.1 The subsampling rule — a first-class command

**Never subsample reads.** A random 10 000-read draw from a 20×-oversequenced library yields ~1 read per UMI and is useless: it destroys the exact structure every stage of this pipeline operates on. The correct operation is *sort by UMI, take the first N distinct UMIs, keep **all** their reads*.

Implementation: `python/migec/subsample.py`, heavy pass in `migec._core.subsample_by_umi`, exposed as:

```
migec subsample --r1 R1.fq.gz [--r2 R2.fq.gz] \
                [--pattern-file barcodes.txt | --from-header] \
                --umis 1000 [--min-reads 1] [--max-reads-per-umi 0] \
                [--order first|random-umi] [--seed 42] \
                [--per-sample] -o example/
```

Semantics — two passes, no random read draw anywhere:

1. Pass 1: extract `(sample, cell, umi)` per read pair (via the checkout pattern, or by parsing the `UMI:` header tag with `--from-header`). Build the distinct-UMI list, in first-appearance order.
2. Select N UMIs: `--order first` = the first N distinct (deterministic, no seed needed); `--order random-umi` = a seeded sample of N **UMIs** (never of reads); `--per-sample` applies N per sample id.
3. Pass 2: write **every** read whose UMI is in the selected set, mates kept in lockstep.
4. Exit assertion (the guard against silently regressing to read-subsampling):
   `mean_reads_per_umi(output) == mean_reads_per_umi(input restricted to the selected UMIs)`, exactly.
5. Emits `subsample.json`: `n_umis`, `n_reads`, `reads_per_umi_mean`, `reads_per_umi_median`, `order`, `seed`, source file sha256s — and this becomes a `SOURCES.md` entry for every example bundle.

Every notebook's data cell calls `migec subsample`, and `docs/faq.rst` names `seqtk sample` / `head -n 40000` as the anti-pattern.

---

## 9. The four markdown files

Markers used throughout (arda convention): 🚫 = never do this · ⚠️ = subtle caveat · ✅ = verified good. Each rule written as **rule + the incident that produced it**.

**`CLAUDE.md`** (committed, authoritative)
1. What this is (one paragraph) + explicit scope fence: *ends at consensus FASTQ; no aligner, no variant calling — that was MAGERI and it is out of scope*.
2. Layout: `include/migec` / `src` / `python/migec` / `tests` tiers / `docs` / `notebooks`.
3. Build, test, run: `./setup.sh`, `cmake -DMIGEC_TESTS=ON && ctest`, `pytest tests/unit tests/synthetic`, `RUN_BENCHMARK=1`, docs gate.
4. Domain conventions: barcode grammar (uppercase seed / lowercase fuzzy / `N` = UMI); master vs slave; consensus header format; the phred/IUPAC acceptance rule.
5. 🚫 Non-negotiables: never random-subsample reads (§8); never drop a no-parent low-coverage UMI (3–5 reads) — demote its quality instead; never report a consensus quality above the RT/PCR error floor; never reproduce the two v1 bugs (quality indexed from read start instead of match offset; the dangling-else that made low-quality mismatches uncountable).
6. Data policy: MIGEC Experiment 1 raw reads must not leave aldan3 — derived/summary artifacts only.
7. **Open loops / next steps** (kept current).
8. Pointers: `ROADMAP.md`, `SOURCES.md`, upstream `seqtree`, `arda`.

**`ROADMAP.md`** — v2.0 (checkout → stats → correct → refine → assemble, ship wheels, docs) · v2.1 (`suggest`, 10x CB/doublets, HTML report) · v2.2 (libdeflate reader, zstd spill, optional seqtree C++ linkage, contig assembly on overlapping mates) · Explicitly out of scope (alignment, variant calling, CDR3 extraction — → arda/vdjtools) · Benchmark plan (MAGERI SRA PRJNA352143, MIGEC PRJNA239303, MIGEC exp1, 10x GEX+VDJ; assets → `huggingface.co/datasets/isalgo/umi_data`).

**`SOURCES.md`** — one table, columns: `name | origin (URL / HF repo+path / aldan3 path) | format | fetch-or-regenerate command | provenance (experimental | derived | simulated) | restrictions`. Seeded from `/Users/mikesh/vcs/projects/2026-arda-benchmark/SOURCES.md`: PRJNA239303 (SRR1200517-20, experimental, public), MIGEC exp1 (`/projects/tcr_bcr_rnaseq/data_migec_exp1/...`, experimental, **internal, raw reads must not leave the cluster**), PRJNA352143/SRR1799908 (experimental, public), the 10x sets, the `tests/data` fixture (derived — with the exact `migec subsample` command), every simulator output (simulated — with the exact `migec simulate` command + seed), and the known-corrupt `scratch/spikein/S1_R2_2M.fq` (past record 1 742 617).

**`CHANGELOG.md`** — Keep-a-Changelog, `## [2.0.0] - YYYY-MM-DD` with a leading **Rewrite** note: v2 is a from-scratch C++20/Python reimplementation; Groovy MIGEC 1.2.9 is branch `legacy-v1` / tag `v1-final`; `CdrBlast`/`CdrFinal` are removed, use arda; MAGERI is folded in for consensus/quality only. Version here is the single source of truth asserted against the release tag by `publish.yml`.

---

## 10. What I would cut from v1

| Cut | Why |
|---|---|
| **seqtree C++ linkage (FetchContent/submodule)** | Correction is one batched, GIL-released `search_batch` per stage. A plain `seqtree>=0.6` runtime dep is identical in speed and removes the entire build-fragility class. Keep the 8-line FetchContent block in ROADMAP. |
| **The upstream seqtree install/export patch** | Only pays off with a second C++ consumer. Not now. |
| **Windows wheels** | zlib on Windows needs vcpkg/vendoring; it is a cluster tool; arda's Windows legs already blocked one release. |
| **libdeflate and zstd** | Real wins, but two extra code paths on day one. Hide the reader behind `migec::GzipReader`; revisit after profiling. |
| **musllinux / manylinux_i686 wheels** | No demand, C++20 toolchain pain. |
| **`cdrblast.rst`, `cdrfinal.rst` and any CDR3 logic** | Out of scope — arda/vdjtools. |
| **Compressed external-sort spill** | Transient files; fixed-size uncompressed records are what makes the sort fast. `--tmp-dir` is enough. |
| **A separate `migec-core` PyPI distribution** | One wheel, one name. |
| **`benchmarks` job in `ci.yml` (seqtree has one)** | Real UMI benchmarks need multi-GB data and aldan3; run them out-of-band, keep only the RSS/throughput regression under `RUN_BENCHMARK=1`. |
| **`migec report` HTML in v1** | The five marimo notebooks already cover it; ship `report` in 2.1. |

---

## 11. Two questions for the user before execution

1. **License** — the current `LICENSE` is a MiLaboratory proprietary agreement, not open source. What licence does v2 ship under? (Recommend GPL-3.0-or-later, matching seqtree/arda, which we depend on.)
2. **Does the v1 Groovy `master` need to stay installable during the transition?** If yes, add one additive commit to `legacy-v1` pointing JitPack/README at tag `v1-final` *before* step 5.