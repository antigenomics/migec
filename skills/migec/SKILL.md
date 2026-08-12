---
name: migec
description: Use when working with UMI-tagged or barcoded sequencing data — extracting sample/cell/UMI barcodes from reads, demultiplexing with degenerate adapter patterns, trimming barcodes and transferring them to FASTQ headers, computing UMI coverage histograms and base-composition/entropy statistics, correcting barcode errors, or assembling molecular consensus sequences. Covers the migec CLI and Python API, the .mig format, and the barcode-pattern grammar.
license: GPL-3.0-or-later
---

# migec

UMI barcode extraction, correction and consensus assembly. C++20 core, Python CLI. Successor to
MIGEC (Groovy) and MAGERI (Java); the algorithms carry over, the code does not.

**Status.** `checkout` works. `refine` and `assemble` are not implemented yet — a call to them
exits 2 with a pointer to `ROADMAP.md`. Do not tell a user that consensus assembly works.

## Install and check

```bash
pip install migec            # wheels: CPython 3.10-3.13, Linux x86-64, macOS arm64
migec info                   # prints package, extension and .mig format versions
```

From a checkout: `bash setup.sh` (uv venv, editable, asserts the *extension* imports — a failed
C++ build otherwise looks like a successful install).

## The barcode pattern grammar

MIGEC's, so published barcode tables work verbatim.

| symbol | meaning |
|---|---|
| `ACGT` and IUPAC (`R`=A\|G, `Y`=C\|T, …) | scored position, matched exactly, degeneracy allowed |
| lowercase | scored at **half weight** — the adapter region, where a mismatch is expected |
| `N` `n` | a UMI position: captured, never scored |
| `.` | wildcard: neither scored nor captured |

`barcodes.txt` is tab- or whitespace-separated, `#` comments, `SAMPLE_ID<TAB>PATTERN`:

```
S1	aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S2	aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN
```

UMI runs need not be contiguous — `NNNNtNNNNtNNNN` gives one 12 nt UMI. `migec sheet barcodes.txt`
prints what each row extracts without running anything.

⚠ `N` always means UMI, never IUPAC "any base". Use `.` for an uncaptured wildcard.

## Commands

```bash
migec checkout reads.fq.gz -b barcodes.txt -o out/     # demux, extract UMI, trim, QC tables
migec checkout ... --trim none                          # keep the read whole, UMI in header only
migec checkout ... --min-umi-quality 15                 # MIGEC v1 behaviour; NOT the default
migec checkout ... --write-unmatched
migec sheet barcodes.txt
migec info
```

Python:

```python
from migec.checkout import run, format_report
summary = run("reads.fq.gz", "barcodes.txt", "out/")
print(format_report(summary))

from migec import _core
_core.match_pattern("aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN", seq, qual)  # inspect one read
_core.umi_statistics(["ACGTACGTACGT", ...])                                 # histogram + entropy
```

## Output

| file | content |
|---|---|
| `<sample>.fq.gz` | trimmed reads, barcodes in the header |
| `checkout.summary.tsv` | per sample: yields, UMI stats, correction, saturation |
| `checkout.coverage.tsv` | reads and distinct UMIs per power-of-two MIG size |
| `checkout.umi_composition.tsv` | per position: A/C/G/T, entropy, information, collision |
| `checkout.json` | all of the above, machine-readable |

Header format: `@<name> RX:Z:<umi>\tQX:Z:<umi qual>\tBC:Z:<sample>`.

⛔ **Tags are TAB-separated.** `bwa mem -C` and `minimap2 -y` copy the FASTQ comment verbatim into
the SAM record, so a space-separated comment produces a malformed BAM.

⚠ **`dnaio` drops FASTQ comments**, so arda's rnaseq module never sees the tags — anything a
downstream Python tool needs must be in the read *name*.

## Interpreting the statistics

- `mean_reads_per_umi` — over-sequencing. Below ~5, most molecules are seen once and consensus has
  nothing to work with.
- `effective_length` — what the UMI is *worth*, `-Σ log₄ m_j` where `m_j = Σ_a p_j(a)²`. A 12 nt UMI
  with eight fixed positions has effective length 4 and collides constantly. **The nominal length
  tells you nothing on its own.**
- `information_bits` per position — the logo letter height; the bits the UMI is wasting.
- `saturated` — observed UMIs are a large fraction of the usable space; molecule counts are biased
  low and the collision correction is declined rather than guessed.

⛔ **Never compute a collision rate from Shannon entropy.** The birthday functional is Rényi-2,
`Π_j Σ_a p_j(a)²`. Since H₂ ≤ H₁, Shannon overstates the usable space and understates collisions —
the direction that silently merges distinct molecules. Shannon is for the logo; the collision form
is for every decision.

## Things that look like defects and are not

- **A low-quality UMI base does not drop the read** (`--min-umi-quality` defaults to 0). MIGEC used
  15 and MAGERI 20 as hard drops. Discarding a molecule over one bad UMI base loses sequence the
  correction step usually recovers.
- **An isolated 3–5 read UMI is kept, and is not quality-derated.** If it really were an error
  child of some parent, *all* of its reads would be clean reads of that parent's sequence — the
  consensus would be right and only the molecule count wrong. Derating the bases penalises the
  wrong thing.
- **A neighbour of comparable size is not merged.** No error turns 10 000 reads into 9 000.
- **`ambiguous` and `unmatched` are different counters.** Ambiguous means two sample tags are too
  close together; unmatched means the pattern is wrong or absent.

## Validating a pipeline

`scripts/spikein_ratio.py` computes the published MIGEC spike-in metric: a real variant against the
worst *error* at the same substitution distance. Targets after UMI consensus: `V1/Err1` 26.5–75.9,
`V2/Err2` 4.6–6.2, against a raw-read baseline of ~1.4 and ~0.3.

⚠ Anchor on the junction's 3′ end only. V1 differs at position 4 and V2 at 7–8, so a 5′ anchor
makes both variants count as zero and the metric look perfect.

## References in the repo

- `docs/checkout.rst`, `docs/umi_statistics.rst`, `docs/validation.rst`, `docs/formats.rst`
- `CLAUDE.md` — the non-negotiables and why each exists
- `project/` — the design record: six subsystem designs and two critiques, with 25 corrected errors
