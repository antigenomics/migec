# Notebooks

[marimo](https://marimo.io) notebooks. Each is a plain Python file — readable in a diff, runnable
without a notebook server, and with no hidden execution order.

## Running them

Every notebook declares its own dependencies in a [PEP 723](https://peps.python.org/pep-0723/)
header, so `uv` builds the environment for you and nothing has to be installed first:

```bash
uv run marimo edit notebooks/platforms.py       # interactive, in a browser
uv run marimo run  notebooks/platforms.py       # read-only app view
uv run python      notebooks/platforms.py       # just execute it, no UI
```

If you already have an environment with migec in it, `marimo edit notebooks/<name>.py` works the
same way.

Note: the headers list `pyarrow` and `pandas` wherever a chart is drawn. That is not padding —
altair renders a polars frame through pyarrow and then pandas, so declaring only `altair` and
`polars` fails on a clean machine three imports deep, with an error that names neither.

## What each one answers

| notebook | question | data | needs network |
|---|---|---|---|
| `platforms.py` | how do I declare my barcode layout? every preset, two run end to end | `isalgo/umi_data` fixtures | yes, first run |
| `barcode_space.py` | is my barcode long enough? collisions, occupancy, the error budget | simulated | no |
| `refine_diagnostics.py` | where did the molecules go? coverage curve, barcode-rank plot, error spectrum | simulated | no |
| `exome_capture.py` | duplicates or real molecules? why coordinate deduplication undercounts a capture panel | simulated | no |
| `airr_repertoire.py` | how much of a repertoire is PCR? clonotype counts from reads against molecules | simulated | no |
| `ctdna_variants.py` | how many molecules does a variant caller actually get, at a certified allele frequency? | SRA, fetched on demand | yes |

The three simulated ones are self-contained: they build a library whose true molecule and clone
counts are known, so every number they print can be checked against the answer rather than merely
admired.

## Where the data comes from

Two sources, on purpose:

* **`isalgo/umi_data`** — small CI fixtures, downloaded on first run by `huggingface-hub`. These
  are subsampled by *whole barcodes* (all reads of a fraction of the molecules), never by reads;
  taking random reads gives one read per UMI and is useless for a UMI tool.
* **SRA, on demand** — anything with a public accession is fetched by `scripts/sra_fetch.py`
  rather than mirrored. `SOURCES.md` records the accession and the exact command for each.

Nothing here writes into the repository. Notebooks use a temporary directory and clean up after
themselves.

## Editing them

marimo stores notebooks as Python, so ordinary review tools work:

```bash
ruff check notebooks/                 # they are linted with everything else
git diff notebooks/platforms.py       # a real diff, not a JSON blob
```

Cells are functions and their order is derived from their dependencies, not from the order you ran
them in. If a cell returns a value it becomes the cell's output; the last expression is what you
see.
