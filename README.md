# migec

<p align="center">
  <a href="https://pypi.org/project/migec/"><img alt="PyPI" src="https://img.shields.io/pypi/v/migec"></a>
  <a href="https://github.com/antigenomics/migec/actions/workflows/ci.yml"><img alt="CI" src="https://github.com/antigenomics/migec/actions/workflows/ci.yml/badge.svg"></a>
  <a href="https://antigenomics.github.io/migec/"><img alt="docs" src="https://github.com/antigenomics/migec/actions/workflows/docs.yml/badge.svg"></a>
  <img alt="python" src="https://img.shields.io/badge/python-3.10%2B-blue">
  <img alt="C++" src="https://img.shields.io/badge/C%2B%2B-20-blue">
  <img alt="license" src="https://img.shields.io/badge/license-GPLv3-green">
</p>

**UMI barcode extraction, correction and consensus assembly for barcoded sequencing data.**

A complete C++20 rewrite of [MIGEC](https://doi.org/10.1038/nmeth.2960) (Shugay et al., *Nature
Methods* 2014) and [MAGERI](https://doi.org/10.1371/journal.pcbi.1005480) (Shugay et al., *PLoS
Computational Biology* 2017).

> **Version 2 is under construction.** All three stages work today, with cell barcodes, whitelists,
> dual-end and positional (10x) layouts, mate merging, cell calling, index-hopping estimation, QC
> figures, and `suggest`/`subsample`/`plot`. The published benchmark comparisons are what remain; see
> [`ROADMAP.md`](https://github.com/antigenomics/migec/blob/master/ROADMAP.md). The Groovy MIGEC 1.2.9 is archived on branch
> [`legacy-v1`](https://github.com/antigenomics/migec/tree/legacy-v1) and at tag `v1-final` — Java users want the jars on the
> [1.2.9 release](https://github.com/antigenomics/migec/releases/tag/1.2.9).

## Why

Tag each molecule with a random barcode before amplification and every read carrying that barcode
descends from one original molecule, so collapsing them into a consensus removes essentially all
sequencing error. The difficulty is entirely in the details: barcodes acquire errors of their own
and an error child has to be told from a genuine collision; a molecule seen three times is still
information rather than something to threshold away; and no consensus can repair an error made
before the first amplification cycle, because it is in every read. migec measures that floor from
the data and refuses to claim a quality above it.

<p align="center"><img alt="the migec pipeline" src="https://raw.githubusercontent.com/antigenomics/migec/master/assets/pipeline.svg" width="42%"></p>

## Install

```bash
uv tool install migec      # or: uv pip install migec, or: pip install migec
migec info                 # prints the three version strings; they must agree
```

Wheels for CPython 3.10–3.13 on Linux x86-64 and macOS arm64, so nothing compiles. On a cluster
whose system Python is older than 3.10, bring your own:

```bash
uv venv --python 3.12 ~/envs/migec && source ~/envs/migec/bin/activate && uv pip install migec
```

From source, for development: `git clone https://github.com/antigenomics/migec && cd migec && bash
setup.sh`. See [installation](https://antigenomics.github.io/migec/installation.html).

## Where the barcode is

This is the one thing migec has to be told. Most libraries put the barcode at a fixed offset in one
read, so that is the primary way to say it — a position, no sheet, no anchor:

```bash
migec checkout reads.fq.gz --bc-pattern '^NNNNNNNN'  -o out/    # 8 nt UMI at the read start
migec checkout reads.fq.gz --bc-pattern '0:8'        -o out/    # the same, as a half-open slice
migec checkout reads.fq.gz --bc-pattern '0:4,5:10'   -o out/    # 9 nt UMI split by one spacer base
migec checkout R1.fq.gz R2.fq.gz --bc-pattern 'cell:0:16,16:26' -o out/     # 10x
```

`N` is a UMI base, `X` a cell-barcode base, and slices are half-open and 0-based like Python's. A
leading `^`, and every slice list, **anchors the barcode at the first base**, so `--max-offset`
never has to be passed. Or name the chemistry — `migec sheet --presets` prints all of them with the
source each layout is written down in:

| preset | layout | |
|---|---|---|
| `umi` | `^NNNNNNNN` | generic inline UMI |
| `migec` | `cagtggtatcaacgcagagtNNNNtNNNNtNNNN` | MIGEC 5'-RACE RepSeq |
| `primerid` | `NNNNNNNNNcagtttaacttttgggccatcca` | HIV-1 Primer ID, as used by MAGERI |
| `duplex` | `^NNNNNNNNNNNN.....` on both mates | duplex sequencing |
| `10x` | `^XXXXXXXXXXXXXXXXNNNNNNNNNNNN` | 10x Chromium 3' v3 |
| `10x-v2` | `^XXXXXXXXXXXXXXXXNNNNNNNNNN` | 10x Chromium 3' v2 and 5' |
| `tso500` | `^NNNNN.....` on R1 | Illumina TSO500 ctDNA — read the warning in [layouts](https://antigenomics.github.io/migec/layouts.html) |
| `smarter-umi` | `^NNNNNNNNNN...` | SMARTer template-switching RNA-seq |

```bash
migec checkout R1.fq.gz R2.fq.gz --preset 10x-v2 -o out/          # a named chemistry
migec checkout R1.fq.gz R2.fq.gz --read-structure 5M5S+T -o out/  # fgbio, Picard, samtools, TSO500
migec checkout reads.fq.gz -b barcodes.txt -o out/                # many samples: MIGEC's own table
migec suggest  reads.fq.gz                                        # if you do not know: read it off the data
```

A preset says where the barcode is. It does not say what a consensus is worth, and that matters
more: the same 12 nt UMI serves a repertoire census and an MRD assay, and the right settings are
opposite. `migec sheet --assay ctdna` prints the second half — the `--min-reads` and the
pre-amplification floor the experiment implies, for eight profiles from `airr` to `mrd`
([assays](https://antigenomics.github.io/migec/assays.html), [layouts](https://antigenomics.github.io/migec/layouts.html)).

## The pipeline

```bash
migec checkout  reads.fq.gz -b barcodes.txt -o out/   # find the barcode, cut it out, demultiplex
migec refine    out/S1.fq.gz -o ref/                  # correct the errors IN the barcode
migec assemble  ref/S1.fq.gz -o cons/                 # one consensus per molecule
migec subsample out/S1.fq.gz -o small.fq.gz --keep 1  # a fixture that is still a library
```

Those four and `suggest` are the pipeline; `plot`, `sheet` and `info` read no reads at all
([commands](https://antigenomics.github.io/migec/commands.html)).

Every stage takes `-t/--threads` (one per core by default) and `--limit-read N` / `--limit-umi N`,
which stop the intake early — for getting an answer out of a 400 GB run in a minute, never as a
sample. **`-t` changes the wall clock and nothing else**: the output is byte-identical at any thread
count, which is what makes a retry on different cores comparable. Every run says what it did, what
it cost, and what the barcode was worth:

```
reads       2,000,000
  assigned  2,000,000 (100.0%)
  unmatched 0 (0.0%)
  ambiguous 0 (0.0%)

1.6 s (1,243,801 reads/s) = 1.5 s matching on 8 threads + 0.1 s UMI statistics
peak RSS 136.0 MB of which UMI counters 11.5 MB

sample             reads        UMIs  reads/UMI  UMI len  eff len
S1               500,000     125,000       4.00       12    12.00
```

## What comes out

Ordinary FASTQ, trimmed of adapter, sample tag and UMI. One record is one molecule, and its identity
is carried twice — in the read **name** (`<sample>.<cell>.<umi>`, for tools that drop FASTQ
comments) and in tab-separated **SAM tags** that survive into a BAM:

```
@r0 RX:Z:GCTAAAGACAAT	QX:Z:IIIIIIIIIIII	BC:Z:S1
TACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTT
```

| output | what |
|---|---|
| `<sample>.fq.gz`, `<sample>.consensus.fq.gz` | the reads, then one record per molecule |
| `checkout.summary.tsv`, `.coverage.tsv`, `.umi_composition.tsv` | yields, MIG sizes, per-position base usage |
| `<sample>.barcodes.tsv`, `.umi_errors.tsv`, `.mig.tsv` | every barcode with its parent, the error rate per depth, every molecule |
| `<stage>.json` | all of it, machine-readable |

```bash
migec plot cons/          # twenty QC panels with gnuplot, straight off those TSVs
```

<p align="center"><img alt="barcode rank plot" src="https://raw.githubusercontent.com/antigenomics/migec/master/assets/cell_rank.svg" width="60%"></p>

**How deeply was each molecule sequenced, and is the barcode any good?** Two panels answer the
questions you ask first, both drawn straight off `checkout`'s own tables on a real 10x library
(`sc5p_v2_hs_PBMC_1k` VDJ-T, 3.16 M read pairs, 311,421 molecules):

<p align="center">
<img alt="MIG size distribution" src="https://raw.githubusercontent.com/antigenomics/migec/master/assets/coverage.svg" width="49%">
<img alt="barcode base composition" src="https://raw.githubusercontent.com/antigenomics/migec/master/assets/umi_pwm.PBMC.svg" width="49%">
</p>

Left, **MIG size in doubling bins** — bars are molecules, the line is the reads inside them. Most
molecules are seen once or twice while most *reads* pile onto a few very deep ones, which is the
distribution that decides whether a consensus is worth assembling at all.

Right, **the barcode's own PWM**, every position of cell barcode and UMI against the 1/4 a free
synthesiser mix would give. The two halves read differently on purpose: the 16 cell positions
(marked on the axis) swing between 0.198 and 0.305 because they are drawn from 10x's whitelist,
while the 10 UMI positions sit tight around 1/4. A UMI position that wandered like a cell position
would be a synthesis defect, and `effective_length` in `checkout.summary.tsv` is what it costs you
— 25.85 usable bases of 26 here.

Each of these was run against real output ([downstream](https://antigenomics.github.io/migec/downstream.html)):

```bash
minimap2 -ax sr -y ref.fa cons/S1.consensus.fq.gz | samtools sort -o S1.bam   # RX, CB, MI in the BAM
minibwa map -y -t8 ref.fa cons/S1.consensus.fq.gz | samtools sort -o S1.bam   # `-y`, not bwa's `-C`
bwa mem -C     ref.fa cons/S1.consensus.fq.gz     | samtools sort -o S1.bam
arda amplicon --r1 cons/S1.consensus.fq.gz -p S1      # AIRR sequence_id IS the molecule id
salmon quant -i tx.idx -l A -r cons/S1.consensus.fq.gz -o quant/   # NumReads are molecule counts
```

> **Never** run alevin, bustools or STARsolo on a consensus FASTQ. They read the barcode out of a
> *raw* barcode read and deduplicate themselves; migec already did, and that read no longer exists.

## What it is worth: rare variants, measured

Commercial cfDNA reference material with **certified** allele frequencies, including a 0%-certified
arm that is a true negative by construction. Three replicates per arm, one panel, one aligner,
matched molecule-support thresholds, substitutions only — migec emits no indels by design and 56%
of UMIErrorCorrect's calls are deletions ([post-processing](https://antigenomics.github.io/migec/postprocessing.html), `assets/ctdna_callers.tsv`).

| pipeline | false calls / sample, 0% arm | 0.125% | 0.25% | 1% | VAF at 1% | median depth |
|---|---|---|---|---|---|---|
| **migec + Mutect2** | **0.67** | 0/3 | 1/3 | **3/3** | 0.0103 | 2,811 molecules |
| **migec + LoFreq** | **2.00** | **1/3** | **3/3** | **3/3** | **0.0102** | 2,832 molecules |
| no consensus + LoFreq | 5.67 | 0/3 | 3/3 | 3/3 | 0.0127 | **52,628 reads** |
| UMIErrorCorrect (its own consensus **and** caller) | 7.67 | 1/3 | 3/3 | 2/2 | 0.0094 | 5,010 |

**migec + LoFreq matches the best sensitivity at every arm and reports 3.8× fewer false positives on
the true negative**; migec + Mutect2 reports the fewest of anything measured and pays for it at
0.25%. Against **MAGERI**, the other descendant of MIGEC 1, on a simulated shallow library where
both tools emit the *same* consensuses at the same accuracy: MAGERI reports 142 variants of which
**137 are at positions nothing was injected at**, migec + LoFreq reports 5 and is right about all
five ([validation](https://antigenomics.github.io/migec/validation.html)).

The **no consensus** row is the same reads, trimming, barcode correction, aligner and caller — only
a record is a read rather than a molecule. Collapsing cuts false positives 2.8×, makes the measured
frequency right (1.02× of certified against 1.27×), and **detects more from 38× less depth**: at
0.125% the consensus finds the hotspot in 1 of 3 replicates and a 197,772× read pileup finds it in
none. A read count is not a molecule count.

Two flags decide more than the choice of caller, and both are measured rather than argued:

```bash
gatk Mutect2 --max-reads-per-alignment-start 0 ...   # or it sees 1.5% of your molecules
migec assemble ... --min-reads 3                     # 0% arm: 10.0 -> 2.0 calls per sample
```

## What makes it different

- Barcode correction uses the evidence that survives at **one read per UMI** — the barcode's own
  base quality and payload agreement, not only the count ratio, which reports zero there
  ([refine](https://antigenomics.github.io/migec/refine.html)).
- Emitted quality is capped at the **measured** RT/first-cycle floor, Q40 by default and fitted
  from the data with `--pre-amp-error auto`, never taken from the instrument's
  ([quality floor](https://antigenomics.github.io/migec/quality_floor.html)).
- Every model-derived number has something model-free beside it: collisions, the barcode error rate,
  the split threshold ([nulls](https://antigenomics.github.io/migec/nulls.html), [barcode space](https://antigenomics.github.io/migec/barcode_space.html)).
- Nothing scales with the library — a range partition into buckets and a sorted counter array, 22 B
  per distinct UMI against a hash map's 48 ([performance](https://antigenomics.github.io/migec/performance.html)).
- Twenty QC panels, each a gnuplot script over a TSV a stage already wrote, so a figure can never
  disagree with the report ([plots](https://antigenomics.github.io/migec/plots.html)).
- For rare variants the **molecule count decides, not the caller**, and it is fixed before any
  software runs ([variants](https://antigenomics.github.io/migec/variants.html), [detection](https://antigenomics.github.io/migec/detection.html)).
- BAM, SAM and CRAM are inputs, so a capture, exome or ctDNA kit that puts the UMI in the **index
  read** never needs `checkout` at all ([bring your own UMI](https://antigenomics.github.io/migec/byo_umi.html)).

## Documentation

<https://antigenomics.github.io/migec/>

| | |
|---|---|
| [Installation](https://antigenomics.github.io/migec/installation.html), [Examples](https://antigenomics.github.io/migec/examples.html) | a copy-paste run per platform, and six marimo notebooks |
| [Layouts](https://antigenomics.github.io/migec/layouts.html), [Assays](https://antigenomics.github.io/migec/assays.html) | where the barcode is, and what a consensus is worth once it is found |
| [Commands](https://antigenomics.github.io/migec/commands.html) | all eight, with the number each one decides |
| [Post-processing](https://antigenomics.github.io/migec/postprocessing.html) | the certified-cfDNA benchmark, then [downstream](https://antigenomics.github.io/migec/downstream.html) tool by tool, [variant calling](https://antigenomics.github.io/migec/variants.html) and [detection limits](https://antigenomics.github.io/migec/detection.html) |
| [Method](https://antigenomics.github.io/migec/method.html) | why every default is what it is: barcode space, nulls, the quality floor |
| [Reference](https://antigenomics.github.io/migec/reference.html) | file formats, speed and memory, pipelines, roadmap |
| [`SOURCES.md`](https://github.com/antigenomics/migec/blob/master/SOURCES.md) | every dataset, where it came from, and the command that re-fetches it |

## Citing

Until the v2 paper exists, cite the original methods:

- Shugay M *et al.* Towards error-free profiling of immune repertoires. *Nat Methods* 11:653–655
  (2014). doi:10.1038/nmeth.2960
- Shugay M *et al.* MAGERI: Computational pipeline for molecular-barcoded targeted resequencing.
  *PLoS Comput Biol* 13(5):e1005480 (2017). doi:10.1371/journal.pcbi.1005480
- Turchaninova MA *et al.* High-quality full-length immunoglobulin profiling with unique molecular 
  barcoding. *Nat Protoc* 11(9):1599-616 (2016). doi: 10.1038/nprot.2016.093

## License

GPL-3.0-or-later. 
