# migec

**UMI barcode extraction, correction and consensus assembly for barcoded sequencing data.**

A complete C++20 rewrite of [MIGEC](https://doi.org/10.1038/nmeth.2960) (Shugay et al., *Nature
Methods* 2014) and [MAGERI](https://doi.org/10.1371/journal.pcbi.1005480) (Shugay et al., *PLoS
Computational Biology* 2017).

> **Version 2 is under construction.** This build ships the `.mig` intermediate format and the
> FASTQ IO layer; `checkout`, `refine` and `assemble` land over the following milestones. The
> Groovy MIGEC 1.2.9 is archived on branch [`legacy-v1`](../../tree/legacy-v1) and at tag
> `v1-final` — Java users want the jars on the [1.2.9 release](../../releases/tag/1.2.9).

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
