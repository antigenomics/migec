# Running migec under SLURM

Two templates. Both run **without** SLURM as ordinary bash scripts, which is how they are tested
and how you should smoke-test a layout before queueing a cohort.

```
migec_sample.sbatch    one sample, all three stages, driven by environment variables
migec_array.sbatch     one array task per row of a sample sheet
samples.tsv            an example sheet
```

## One sample

```bash
sbatch --export=ALL,R1=s1_R1.fq.gz,R2=s1_R2.fq.gz,SAMPLE=s1,PRESET=10x-v2,PAYLOAD_MATE=2 \
       migec_sample.sbatch

# the same thing on a laptop, to check the layout before queueing anything
R1=s1_R1.fq.gz R2=s1_R2.fq.gz SAMPLE=s1 PRESET=10x-v2 PAYLOAD_MATE=2 bash migec_sample.sbatch
```

| variable | meaning |
|---|---|
| `R1`, `R2` | input FASTQ; `R2` optional |
| `SAMPLE` | sample id; defaults to the `R1` basename |
| `OUT` | output root; defaults to `migec_out/$SAMPLE` |
| `PRESET` / `BC_PATTERN` / `READ_STRUCTURE` | where the barcode is — set exactly one |
| `PAYLOAD_MATE` | which mate still carries sequence after trimming. `2` for 10x |
| `RT_ERROR` | pre-amplification floor: `rt` (default, Q40), `medium`, `high`, or a rate |

## A cohort

```bash
sbatch --array=1-$(($(wc -l < samples.tsv) - 1)) migec_array.sbatch samples.tsv
```

The sheet is a header plus one row per sample; `layout` is a preset name, a pattern (`^NNNNNNNN`)
or a slice list (`0:12`), and the script tells them apart. Use `-` for an absent R2.

**Note: task 1 is the first data row.** The header is skipped rather than counted, so the array
range is `1-(rows - 1)`. Passing `--array=0-N` instead runs a task that reads the header as a
sample, and it fails in a confusing place rather than at the sheet.

## Sizing the request

The three stages are bounded by different things, which is why they are separate processes:

| stage | scales with | measured | give it |
|---|---|---|---|
| `checkout` | reads | ~1.5 M reads/s at 16 threads | cores |
| `refine` | **distinct barcodes**, not reads | ~1.5 M reads/s; holds the barcode table in memory | memory; `table_bytes` in its JSON sizes the next run |
| `assemble` | one bucket at a time | ~2.5 M reads/s; peak RSS set by the bucket count, not the library | cores |

16 cores and 32 GB covers a typical targeted or single-cell library. For a NovaSeq-scale run, read
`table_bytes` and `peak_rss_bytes` from a pilot rather than guessing — both are reported in every
stage's JSON summary.

**Never: `-t` changes the wall clock and nothing else.** Every stage is byte-identical at any
thread count, so a retry on a different node, an escalating `--requeue`, or a resumed pipeline
cannot produce a result that disagrees with the first attempt. That property is what makes
automatic retries safe here, and it is asserted in `tests/benchmark/`.

## Getting migec onto the node

The cluster's system Python is often too old — migec needs 3.10 or newer. Bring your own
environment rather than fighting the module system:

```bash
# once, on a login node or in the job
uv venv --python 3.12 ~/envs/migec && source ~/envs/migec/bin/activate && uv pip install migec

# or, if uv is not available
python3.12 -m venv ~/envs/migec && ~/envs/migec/bin/pip install migec
```

Wheels are published for Linux x86_64 and macOS arm64 on CPython 3.10-3.13, so this does not
compile anything. Add `source ~/envs/migec/bin/activate` to the top of the sbatch file, or submit
with `--export=ALL,PATH=$HOME/envs/migec/bin:$PATH`.

## Fetching public data first

`scripts/sra_fetch.py` pulls runs from SRA on demand, so a benchmark cohort does not have to be
mirrored anywhere:

```bash
python scripts/sra_fetch.py probe SRR17220921        # read structure, before downloading anything
python scripts/sra_fetch.py get SRR17220921 -o data/ # NCBI S3, ~7 MB/s on 8 connections
```
