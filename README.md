# migec

**UMI barcode extraction, correction and consensus assembly for barcoded sequencing data.**

A complete C++20 rewrite of [MIGEC](https://doi.org/10.1038/nmeth.2960) (Shugay et al., *Nature
Methods* 2014) and [MAGERI](https://doi.org/10.1371/journal.pcbi.1005480) (Shugay et al., *PLoS
Computational Biology* 2017).

> **Version 2 is under construction.** All three stages work today — `checkout`, `refine` and
> `assemble` — with cell barcodes, whitelists, dual-end and positional (10x) layouts, cell calling
> and `suggest`/`subsample`. Index hopping, `.mig` bucket output and the published benchmark
> comparisons are what remain; see [`ROADMAP.md`](ROADMAP.md). The Groovy MIGEC 1.2.9 is archived on
> branch [`legacy-v1`](../../tree/legacy-v1) and at tag `v1-final` — Java users want the jars on the
> [1.2.9 release](../../releases/tag/1.2.9).

## Why

Tag each molecule with a random barcode before amplification and every read carrying that barcode
descends from one original molecule. Collapsing them into a consensus removes essentially all
sequencing error — which is what makes rare-variant detection and error-free repertoire profiling
possible. The difficulty is entirely in the details:

- **Barcodes acquire errors too.** Distinguishing an error-child barcode from a genuine collision
  needs the birthday bound, the base qualities, *and* the fact that a polymerase error in an early
  PCR cycle carries high quality in every read that inherits it. Treating that as a sequencing
  error is the dominant residual mistake in UMI counting.
- **A molecule seen three times is still information.** Cutting at a coverage threshold throws away
  real sequence. migec keeps low-coverage molecules that have no plausible parent and reports the
  uncertainty rather than deleting the data.
- **Consensus cannot fix an error made before amplification.** An RT or first-cycle PCR error is in
  every read. Any quality above that floor is a fiction, so migec measures the floor from the data
  and refuses to claim more.

## Pipeline

```
FASTQ ──checkout──▶ tagged FASTQ ──refine──▶ corrected ──assemble──▶ consensus FASTQ
          │                                 │                              │
     suggest                       barcode table, QC              per-molecule tables
```

Output is ordinary FASTQ with sample, cell barcode and UMI in the read name and in SAM-style tags,
so `bwa-meme`, `minimap2` and [arda](https://github.com/antigenomics/arda) consume it directly.

## Install

```bash
pip install migec
```

Wheels for CPython 3.10–3.13 on Linux x86-64 and macOS arm64. From source: `bash setup.sh`.

## Usage

Barcode tables are MIGEC's, read verbatim — uppercase is matched exactly (IUPAC degeneracy
allowed), lowercase is the fuzzy adapter region, `N` marks a UMI position, and UMI runs need not be
contiguous:

```
S1	aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S2	aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN
```

`X` is a cell-barcode position, captured separately from the UMI. Column 3 is MIGEC's *slave*
pattern — a second pattern on the other mate whose captured positions **extend** the UMI, which is
how a 24 nt dual-end barcode is declared:

```
S1	NNNNNNNNNNNNtgact	agtcaNNNNNNNNNNNN
```

```bash
migec suggest reads.fq.gz                            # where is the barcode? read it off the data
migec sheet barcodes.txt                             # what will each row extract?
migec checkout reads.fq.gz -b barcodes.txt -o out/
migec checkout R1.fq.gz R2.fq.gz -b barcodes.txt -o out/ -t 8
migec refine out/S1.fq.gz -o ref/                    # correct barcode errors
migec assemble ref/S1.fq.gz -o cons/                 # one consensus per molecule
migec subsample out/S1.fq.gz -o small.fq.gz --keep 1 # a fixture that is still a library
```

```
reads       2,000,000
  assigned  2,000,000 (100.0%)
  unmatched 0 (0.0%)
  ambiguous 0 (0.0%)

2.2 s (903,599 reads/s) = 1.5 s matching on 8 threads + 0.7 s UMI statistics, serial
peak RSS 131.0 MB of which UMI counters 21.2 MB

sample             reads        UMIs  reads/UMI  UMI len  eff len
S1               500,000     125,000       4.00       12    12.00
S2               500,000     125,000       4.00       12    12.00
```

Paired input searches both mates for the tag and swaps the pair so R1 always carries it — an
amplicon library sequenced in both orientations otherwise loses half of each MIG at consensus, and
nothing upstream reports it.

Reads come out trimmed of adapter, sample tag and UMI, with the barcode carried in SAM-style tags
that survive `bwa mem -C` and `minimap2 -y` into the BAM:

```
@r0 RX:Z:GCTAAAGACAAT	QX:Z:IIIIIIIIIIII	BC:Z:S1
TACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTT
```

alongside `checkout.summary.tsv`, `checkout.coverage.tsv` (the MIG size histogram) and
`checkout.umi_composition.tsv` (per-position base usage, entropy, information content).

### Positional chemistries, and one pattern instead of a sheet

10x, dual-end and any other layout where the chemistry fixes the barcode's position needs no
barcode table and no anchor:

```bash
migec checkout R1.fq.gz R2.fq.gz --bc-pattern XXXXXXXXXXXXXXXXNNNNNNNNNN --max-offset 0 -o co/
```

`X` is a cell-barcode position, `N` a UMI position — the interface `umi_tools`, `umitools` and
`mgatk` all take. Note: `umi_tools` spells a cell barcode `C`, which is **cytosine** here; pasting one
is refused with the translation rather than compiled into a pattern that matches nothing.

`--max-offset 0` is not a convenience. A 10x barcode has no constant sequence anywhere, so a free
scan cannot place it and correctly refuses; anchored, the placement is the chemistry's and the bar
does not apply. On `sc5p_v2_hs_PBMC_1k` VDJ-T: **100% of 3,155,166 reads assigned**, 221,024
barcodes at 14.28 reads each, 813 cells called.

Note: On such a chemistry `refine` and `assemble` take **R2**: the barcode read is 26 nt and nothing
else, so R1 has no payload at all.

### It corrects the barcodes, with the evidence that survives at one read

```bash
migec refine out/S1.fq.gz -o ref/
```

A barcode one substitution from another is either an error child of it or an independent molecule.
The count ratio separates them on a deep amplicon and is worth **nothing** at 1–3 reads per UMI, so
`refine` also uses the barcode's own base quality at the position that differs — `checkout` already
writes it to `QX` — and **payload agreement**, since an error child is a read of the parent's
molecule. Agreement is worth `log(1/clonality)`, and the clonality is measured rather than assumed:

```
barcodes    23,910 distinct
  merged    3,855 (16.1%) into a parent, 3,889 reads moved
molecules   20,055 after correction          <- 20,000 were simulated

barcode error   2.87e-03 per base            <- 3.0e-03 injected
clonality       0.0100 of random barcode pairs carry the same payload anyway
                -- payload agreement is worth about 100x odds towards the same molecule here
```

Note: At ~1 read per UMI **80% of barcode errors cannot be fixed by anyone** — the parent barcode was
never sequenced. migec corrects a tenth of the rest and destroys no real molecule at any depth
measured, which is the side to err on: a wrong merge deletes a molecule and nothing downstream can
tell, while a missed correction only inflates a count. Every corrected read keeps what it was in an
`OX:Z:` tag. See [docs/refine.rst](docs/refine.rst).

### It collapses each molecule, and caps what it claims

```bash
migec assemble out/S1.fq.gz -o cons/
migec assemble out/S1.fq.gz -o cons/ --contig     # random-primed reads that tile the molecule
```

A molecule is **sample + cell barcode + UMI**, never the UMI alone — the same UMI in two cells is
two molecules, and that is the design rather than a defect. Reads are range partitioned on the
packed key into `.mig` buckets and one bucket is sorted at a time, so nothing scales with the
library: 531,365 reads/s, and 121 MB at 16 buckets against 203 MB at one.

The per-column posterior is `LL[j][b] = Σ_i (r==b ? log(1−e) : log(e/3))`, and then the number that
matters:

```
Q(j) = −10 log10( p_cons(j) + p_floor )
```

The floor is **added, not compared**. An RT or first-cycle-PCR error is in every read and no
consensus removes it, so `assemble` never emits a quality above ~Q38 — the floor measured in
[docs/quality_floor.rst](docs/quality_floor.rst), not the 1e-6 that gets assumed.

`--contig` is for random-primed libraries, where reads sharing a barcode tile the molecule instead
of starting at the same base. They are placed against each other by seed matching, cut into overlap
components, and one consensus is emitted per component — a component is **never** extended across a
gap, because 27.3% of 10x groups hold more than one and a single consensus over those asserts
sequence no read covers. Assembling a cell's full receptor and calling doublets is
[arda](https://github.com/antigenomics/arda)'s job, not this one.

Note: Contig assembly needs a barcode that is not saturated: two fragments of two *different* molecules
sharing a barcode have no sequence in common, which is exactly what two fragments of one look like.
`assemble` runs the same birthday arithmetic on the barcodes it saw and reports how many molecules
a group holds on average — above 1, the warning says so. See [docs/assemble.rst](docs/assemble.rst).

### It tells you where the barcode is

`migec suggest` reads the layout off the reads rather than off the protocol. A UMI cycle is one the
synthesiser mixed — all four bases near 1/4, ~2 bits. A constant cycle is one base near 100%.
Everything else is payload.

```
 cycle      A      C      G      T  1/4 dev     Q  layout
     0  0.271  0.205  0.257  0.267    0.045    33  N  UMI
     ...
     9  0.020  0.971  0.004  0.006    0.721    37  |  constant

segments:
    0-8   umi         9 nt  (mean 1/4 deviation 0.038)
    9-31  constant   23 nt  (mean 1/4 deviation 0.718)  CAGTTTAACTTTTGGGCCATCCA

pattern  NNNNNNNNNcagtttaacttttgggccatcca
```

That is a real HIV Primer ID library (SRR1763769) with nothing supplied but the FASTQ. The pattern
pastes straight into a barcode table, and checking it out assigns 95.0% of reads.

### It tells you whether the barcode was big enough

A 12 nt UMI is 4¹² = 16,777,216 sequences — if the synthesiser delivered exactly 25% of each base.
It never does, so the usable space is the **collision** (Rényi-2) entropy, `1 / Π_j Σ_a p_j(a)²`,
never Shannon: `H₂ ≤ H₁`, so Shannon overstates the space and understates collisions, which is the
direction that silently merges molecules.

From there the birthday problem, in the form that survives a full space:

```
occupied = S·(1 − e^−λ)     molecules = S·λ     P(k>1 | k≥1) = (1 − e^−λ − λe^−λ)/(1 − e^−λ)
```

```
sample              space  occupancy  MIGs >1 mol   molecules   err pred   err est
CTRL              250,902      49.9%        30.6%     173,482    2.0e-03   2.7e-04

warning: CTRL: 31% of MIGs hold more than one molecule (50% of a 250,902 barcode space is
  occupied). Their consensus is a mixture of templates, not a molecule
warning: CTRL: the barcode error estimate (2.7e-04) is not reliable here -- 50% of each barcode's
  1-substitution neighbourhood is itself occupied ...
```

`scripts/collision_check.py` checks that prediction against something model-free — two molecules
sharing a barcode with *different sequences* are visible in the reads — and finds 1.86× more
collisions than predicted. That is *not* the position-independence assumption in `Π_j m_j`, which a
permutation puts at 1.04× — it is the read threshold, since a collided barcode carries two
molecules' reads and is over-represented among the MIGs big enough to show a split.

The barcode error rate is estimated from the distance-1 excess and reported next to what the
reported Phred and the polymerase predict. Note: The estimator has a working range: it recovers 0.92×
of an injected rate at 0.3% occupancy and 0.23× at 50%, always collapsing *downward*, so it is
flagged unreliable past 5% neighbourhood occupancy rather than quietly believed.

Full derivations in [docs/barcode_space.rst](docs/barcode_space.rst); `notebooks/barcode_space.py`
draws it.

### Every derivation has a permutation that checks it

`scripts/permutation_nulls.py` measures three quantities the pipeline otherwise derives, assuming
nothing. On the same HIV library — 125,369 distinct 9 nt barcodes, 47.8% occupancy:

| derived from a model | measured by permutation | |
|---|---|---|
| positions are independent | 1.01× excess, purely nearest-neighbour | holds, to ~1% |
| distance-1 pairs are error children | 97% are chance; ~18,000 are real | permute the background |
| split a MIG at nominal `p < 0.01` | 1% false positives at `-log10 p` = 8.68 | **19× over-call** |

The independence null is a *distribution* — the product measure `q(u) = Π_j p_j(u_j)` — so it is
tested with Jensen-Shannon divergence against a same-size draw from `q`, not with one functional of
each. The residual dependence is entirely between **adjacent** positions, and it has a cause: 0.55%
of reads carry a barcode one base short, a coupling step that did not fire, which frameshifts
everything after it.

The last one is the one that mattered. Reads are not exchangeable — a low-quality read carries a
minor base at many positions at once and looks exactly like a linked subclone — so the null has to
preserve *both* margins of the reads × positions matrix, per-position error count and per-read
error load. The nominal threshold calls 30.62% of MIGs as two molecules; the permutation calls
1.60%. The threshold is a Monte Carlo estimate and its error is quoted: **8.68, bootstrap 95% CI
[8.42, 9.14]** over 82,800 randomisations — a tail quantile from a tenth as many gave 9.61 and
11.66, so the interval is the number, not the point. See [docs/nulls.rst](docs/nulls.rst).

### Most UMIs have 1–3 reads, and that is the normal case

Bulk repertoire profiling and shallow 3' single-cell both put the MIG size histogram's mass at 1–3
reads. migec runs there and says what it can support rather than quoting a number calibrated on a
deep library: the split threshold is inert (a pair of columns can carry at most `log10 C(n, n/2)`,
so it needs ~30 reads), the count-ratio error-child null has no dynamic range, and **nothing is
thresholded away** — `--min-reads` defaults to 1, because a molecule seen once is still a molecule
and the answer to a barcode error is to correct it, not to delete it.

```
    MIG size      groups    share
           1      31,888    79.4%
         2-3       8,176    20.4%
         4-7         112     0.3%

warning: 79.4% of molecules were seen once. A consensus over one read is that read --
  the UMI is buying counting here, not error correction
```

It is also the memory-hostile case, because distinct barcodes are what everything scales with, so
it is what the benchmarks use: 190,595 reads/s at 1.02 reads/UMI, 259 B resident per distinct
barcode, still bounded by the bucket rather than the library.

### Speed and memory are reported, not assumed

`--threads` defaults to one per core and **the output is byte-identical whatever it is set to** —
reads are matched in chunks and written back in input order, so `-t` changes the wall clock and
nothing else.

| threads | reads/s | matching reads/s | peak RSS |
|---|---|---|---|
| 1 | 193,002 | 206,803 | 52 MB |
| 8 | 903,599 | 1,309,576 | 131 MB |
| 16 | **1,055,543** | **1,655,889** | 220 MB |

2 M single-end 129 nt reads, four barcode patterns, 4 reads per molecule, M-series laptop. Two
things had to be true for the matching to scale. zlib compresses random DNA at **7 MB/s** at its
default level 6, so compression runs on the workers (concatenated gzip members are a valid gzip
stream) at level 1 — 137 MB/s for 13% more bytes. And the log-likelihood score tabulates into
1.2 kB, because the `log2` in the inner loop was 90% of runtime.

Both columns are reported, because they scale with different things and only one of them threads:
matching scales with reads, the UMI statistics with *distinct* UMIs at ~1.5–2 µs each, on one
thread. That serial tail is why the end-to-end column flattens at 16 while matching is still
climbing.

The UMI counters are a sorted `(key, count)` array rather than a hash map: **~22 bytes per distinct
UMI against ~48**, which at the 4·10⁸ distinct UMIs of an ordinary NovaSeq run is 8.8 GB against
19 GB. That still does not fit a laptop — the counters are not yet partitioned, so checkout warns
when they pass 1 GB rather than letting you find out from the OOM killer. See
[`docs/performance.rst`](docs/performance.rst).

### Grouping accuracy is measured against Calib

`scripts/compare_calib.py` scores both tools' read partitions against a known truth with the
adjusted Rand index, and reports **splitting and merging separately** — splitting inflates the
molecule count and is recoverable, merging mixes molecules and destroys real variants.

| UMI | UMI error | ARI | reads split | reads merged |
|---|---|---|---|---|
| 12 nt | 0 | **1.0000** | 0.0000 | 0.0000 |
| 12 nt | 5·10⁻³ | 0.9348 | 0.5165 | 0.0004 |
| 6 nt | 0 | 0.8877 | 0.0000 | 0.3982 |

Calib clusters on barcode *and* sequence; migec today groups on the barcode alone, and the gap is
exactly the collision rate — which `eff len` predicts before any clustering runs. A clean 12 nt
barcode needs nothing cleverer; a 6 nt one cannot be rescued by any amount of barcode cleverness,
only by sequence, which is what `assemble` adds.

### `eff len` is the number to look at

A 12 nt UMI is not worth 12 nt unless its bases are uniform. `effective_length` is
$-\sum_j \log_4 \sum_a p_j(a)^2$ — what the barcode is actually worth. A 12 nt UMI with eight
fixed positions has an effective length of 4, a usable space of 256, and will collide constantly.

The distinction matters more than it looks: a sequence logo draws *Shannon* entropy, but the
probability two molecules collide is the *Rényi-2* (collision) entropy. Since H₂ ≤ H₁, using
Shannon overstates the usable space and understates collisions — the direction that silently merges
distinct molecules. Both are reported; only the collision form feeds any decision.

## Pipelines

`integrations/nextflow/migec/` is an nf-core-style local module — `main.nf`, `meta.yml`,
`nextflow.config`, `environment.yml` — that drops into
[nf-core/airrflow](https://nf-co.re/airrflow) or anything else that hands you FASTQ pairs. SLURM is
the pipeline's business, not the module's: it declares `label` and `task.cpus` and nothing more.

Note: Only `checkout` threads. `refine` and `assemble` are single-threaded by construction, so ask for
the cores `checkout` can use and no more.

## Documentation

<https://antigenomics.github.io/migec/> — see [`docs/formats.rst`](docs/formats.rst) for the on-disk
format, and [`ROADMAP.md`](ROADMAP.md) for what is implemented.

## Citing

Until the v2 paper exists, cite the original methods:

- Shugay M *et al.* Towards error-free profiling of immune repertoires. *Nat Methods* 11:653–655
  (2014). doi:10.1038/nmeth.2960
- Shugay M *et al.* MAGERI: Computational pipeline for molecular-barcoded targeted resequencing.
  *PLoS Comput Biol* 13(5):e1005480 (2017). doi:10.1371/journal.pcbi.1005480

## License

GPL-3.0-or-later. The archived v1 code on `legacy-v1` remains under its original MiLaboratory
non-commercial license.
