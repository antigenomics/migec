# migec

**UMI barcode extraction, correction and consensus assembly for barcoded sequencing data.**

A complete C++20 rewrite of [MIGEC](https://doi.org/10.1038/nmeth.2960) (Shugay et al., *Nature
Methods* 2014) and [MAGERI](https://doi.org/10.1371/journal.pcbi.1005480) (Shugay et al., *PLoS
Computational Biology* 2017).

> **Version 2 is under construction.** `checkout` works today — barcode extraction, trimming,
> header transfer and the UMI statistics. `refine` and `assemble` land over the following
> milestones. The Groovy MIGEC 1.2.9 is archived on branch [`legacy-v1`](../../tree/legacy-v1) and
> at tag `v1-final` — Java users want the jars on the [1.2.9 release](../../releases/tag/1.2.9).

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
FASTQ ──checkout──▶ .mig ──refine──▶ .mig + .pumi ──assemble──▶ consensus FASTQ
          │                  │                                        │
     suggest            QC tables, plots                     per-molecule tables
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

```bash
migec sheet barcodes.txt                             # what will each row extract?
migec checkout reads.fq.gz -b barcodes.txt -o out/
migec checkout R1.fq.gz R2.fq.gz -b barcodes.txt -o out/ -t 8
```

```
reads       2,000,000
  assigned  2,000,000 (100.0%)
  unmatched 0 (0.0%)
  ambiguous 0 (0.0%)

2.2 s (917,339 reads/s on 8 threads), peak RSS 204.0 MB of which UMI counters 44.4 MB

sample             reads        UMIs  reads/UMI  UMI len  eff len
S1               500,000     492,751       1.01       12    12.00
S2               500,000     492,606       1.02       12    12.00
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

### Speed and memory are reported, not assumed

`--threads` defaults to one per core and **the output is byte-identical whatever it is set to** —
reads are matched in chunks and written back in input order, so `-t` changes the wall clock and
nothing else.

| threads | reads/s | peak RSS |
|---|---|---|
| 1 | 152,822 | 139 MB |
| 8 | 917,339 | 204 MB |
| 16 | **1,178,648** | 295 MB |

2 M single-end 115 nt reads, four barcode patterns, M-series laptop. Two things had to be true for
that to scale. zlib compresses random DNA at **7 MB/s** at its default level 6, so compression runs
on the workers (concatenated gzip members are a valid gzip stream) at level 1 — 137 MB/s for 13%
more bytes. And the log-likelihood score tabulates into 1.2 kB, because the `log2` in the inner loop
was 90% of runtime.

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
