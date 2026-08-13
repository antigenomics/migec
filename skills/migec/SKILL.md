---
name: migec
description: Use when working with UMI-tagged or barcoded sequencing data — extracting sample/cell/UMI barcodes from reads, demultiplexing with degenerate adapter patterns, trimming barcodes and transferring them to FASTQ headers, computing UMI coverage histograms and base-composition/entropy statistics, correcting barcode errors, or assembling molecular consensus sequences. Covers the migec CLI and Python API, the .mig format, and the barcode-pattern grammar.
license: GPL-3.0-or-later
---

# migec

UMI barcode extraction, correction and consensus assembly. C++20 core, Python CLI. Successor to
MIGEC (Groovy) and MAGERI (Java); the algorithms carry over, the code does not.

**Status.** `checkout`, `suggest`, `refine` and `assemble` work. `subsample` is not implemented yet
— a call to it exits 2 with a pointer to `ROADMAP.md`.

## Install and check

```bash
pip install migec            # wheels: CPython 3.10-3.13, Linux x86-64, macOS arm64
migec info                   # prints package, extension and .mig format versions
```

From a checkout: `bash setup.sh` (uv venv, editable, asserts the *extension* imports — a failed
C++ build otherwise looks like a successful install).

## Telling migec where the barcode is

Four ways, in the order to reach for them. **Positional is the primary mode** — most libraries fix
the barcode at an offset in one read, and that needs no sheet and no anchor flag.

```bash
migec checkout reads.fq.gz --bc-pattern '^NNNNNNNN' -o out/   # 8 nt UMI at the read start
migec checkout reads.fq.gz --bc-pattern '0:8'       -o out/   # the same, as a half-open slice
migec checkout R1.fq.gz R2.fq.gz --bc-pattern 'cell:0:16,16:26' -o out/            # 10x
migec checkout R1.fq.gz R2.fq.gz --preset 10x-v2                -o out/            # named
migec checkout R1.fq.gz R2.fq.gz --read-structure 5M5S+T -o out/                   # fgbio
migec checkout reads.fq.gz -b barcodes.txt -o out/                                 # many samples
```

Slices are half-open and 0-based like Python's: `0:8` is eight bases, the next slice may start at
8, each is a UMI slice unless prefixed `cell:`, and they must not overlap. Gaps become `.`.

**Never: `--max-offset` is now automatic and should not be passed.** A leading `^`, a slice list
and a read structure all anchor the barcode at the first base, and a pattern with nothing to score
is anchored even without them. Passing `--max-offset -1` on such a layout reinstates the old
failure: a free scan has no evidence to choose an offset with and `compile()` refuses.

Presets (`migec sheet --presets` prints each with its source): `umi`, `migec`, `primerid`,
`duplex`, `10x`, `10x-v2`, `tso500`, `smarter-umi`. Never: `duplex` extracts the tags and emits
**single-strand** consensuses; duplex pairing is not implemented, so no duplex error rate may be
quoted from it.

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

`X`/`x` is a **cell barcode** position, captured separately (migec's one extension to the dialect;
`X` is not IUPAC, so published tables keep their exact meaning). Column 3 is the **slave** pattern,
on the other mate, extending the UMI — `S1<TAB>NNNNNNNNNNNNtgact<TAB>agtcaNNNNNNNNNNNN` is a 24 nt
dual-end UMI. Never: Both halves must match or the read is unmatched.

**Note: Why a short handle must be anchored.** Five bases are 10 bits; the acceptance bar is
`log2(offsets × patterns / α)`, which is 12.6 bits over a free scan of a 77 nt read and 6.6 when
anchored. A free scan refusing it is correct — `TGACT` occurs by chance every kilobase. Write the
handle with a `^` (or leave the pattern purely positional) and the anchor is applied for you.

Note: `N` always means UMI, never IUPAC "any base". Use `.` for an uncaptured wildcard.

## Commands

```bash
migec checkout reads.fq.gz -b barcodes.txt -o out/      # demux, extract UMI, trim, QC tables
migec checkout R1.fq.gz R2.fq.gz -b bc.txt -o out/      # paired; tag found in either mate
migec checkout ... -t 8                                 # threads; output identical at any -t
migec checkout ... --trim none                          # keep the read whole, UMI in header only
migec checkout ... --min-umi-quality 15                 # MIGEC v1 behaviour; NOT the default
migec checkout ... --write-unmatched
migec refine out/S1.fq.gz -o ref/                      # correct barcode errors, rewrite RX
migec refine out/S1.fq.gz -o ref/ --min-posterior 0.99 # correct less
migec refine out/S1.fq.gz -o ref/ --no-payload         # what the count ratio alone would do
migec assemble ref/S1.fq.gz -o cons/                   # one consensus per molecule
migec assemble out/S1.fq.gz -o cons/ --contig          # random-primed reads tiling a molecule
migec assemble out/S1.fq.gz -o cons/ --rt-error 3e-5   # a measured floor for THIS chemistry
migec suggest reads.fq.gz -o out/                      # where is the barcode? read it off the data
migec suggest reads.fq.gz --cycles 40 --umi-deviation 0.18
migec sheet barcodes.txt
migec info
```

`suggest` before `checkout` when the layout is unknown or the protocol description is doubtful. It
segments the per-cycle base composition: a UMI cycle is one the synthesiser mixed (all four bases
near 1/4, ~2 bits), a constant cycle is one base near 100%, the rest is payload. It prints a
paste-ready pattern. Note: It stops the pattern at the last *constant* run — a uniform stretch with no
anchor after it is what diverse payload looks like, and claiming it would give a pattern that
matches everywhere. If the UMI is genuinely 3', raise `--cycles`. For paired data with no UMI in R1,
try R2; the note says so.

Python:

```python
from migec.checkout import run, format_report
summary = run("reads.fq.gz", "barcodes.txt", "out/", reads2=None, threads=0)
print(format_report(summary))

from migec.suggest import run as suggest_run, format_report as suggest_report
print(suggest_report(suggest_run("reads.fq.gz", cycles=40)))

from migec import _core
_core.match_pattern("aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN", seq, qual)  # inspect one read
_core.umi_statistics(["ACGTACGTACGT", ...])                                 # histogram + entropy
_core.suggest("reads.fq.gz", cycles=60)                                     # per-cycle PWM + pattern
```

## Output

| file | content |
|---|---|
| `<sample>.fq.gz` | trimmed reads, barcodes in the header (`_R1`/`_R2` when paired) |
| `checkout.summary.tsv` | per sample: yields, UMI stats, correction, saturation |
| `checkout.coverage.tsv` | reads and distinct UMIs per power-of-two MIG size |
| `checkout.umi_composition.tsv` | per position: A/C/G/T, entropy, information, collision |
| `checkout.barcode_space.tsv` | nominal vs effective space, occupancy, λ, molecules, `p_multi`, the error budget |
| `checkout.umi_quality.tsv` | reported Phred histogram over barcode bases |
| `checkout.json` | all of the above, machine-readable |
| `suggest.cycles.tsv` | per cycle: A/C/G/T, entropy, collision, 1/4 deviation, mean Q |
| `suggest.segments.tsv` | the UMI / constant / variable runs and their consensus |

Header format: `@<name> RX:Z:<umi>\tQX:Z:<umi qual>\tBC:Z:<sample>`. Both mates carry the tags.

`summary` also carries `wall_seconds`, `reads_per_second`, `threads`, `peak_rss_bytes` and
`umi_memory_bytes`; `format_report` prints them.

**Never: Tags are TAB-separated.** `bwa mem -C` and `minimap2 -y` copy the FASTQ comment verbatim into
the SAM record, so a space-separated comment produces a malformed BAM.

**Note: `dnaio` drops FASTQ comments**, so arda's rnaseq module never sees the tags — anything a
downstream Python tool needs must be in the read *name*. The consensus name is
`<sample>.<cell>.<umi>` for exactly this reason, and it arrives in arda's AIRR TSV `sequence_id`
unchanged.

### Downstream (measured, `docs/downstream.rst`)

| tool | command | what arrives |
|---|---|---|
| minimap2 | `minimap2 -ax sr -y ref.fa cons.fq.gz` | name + all tags |
| bwa / bwa-mem2 | `bwa mem -C ref.fa cons.fq.gz` | name + all tags |
| arda | `arda amplicon --r1 cons.fq.gz -p out` | name only; `sequence_id` is the molecule id |
| STAR | name truncated at whitespace, comment dropped | name only |
| salmon / kallisto | plain quant, no UMI mode | sequence only |

**Never: Do not run alevin, bustools or STARsolo on a consensus FASTQ.** They read the barcode out
of a *raw* barcode read and deduplicate themselves; migec already did, and the barcode read no
longer exists. One consensus is one molecule, so a plain `salmon`/`kallisto` count already is a
molecule count.

## refine

**Never: Err on precision, never recall.** A wrong merge deletes a molecule and nothing downstream can
tell; a missed correction only inflates a count. `--min-posterior` (0.95) is the knob.

Three pieces of evidence. The **count ratio** is the whole game on a deep amplicon and worth
nothing at 1–3 reads/UMI. The barcode's **base quality** at the differing position separates a
miscall (low Phred there) from an early-PCR child (high Phred in every read). **Payload agreement**
— an error child is a read of the parent's molecule — is worth `log(1/clonality)`, with the
clonality measured from random barcode pairs. It lifts the count gates, which is what makes a
singleton merge possible; disagreement refuses a merge the counts would have made.

**Note: At ~1 read/UMI, ~80% of barcode errors have no observable parent** and cannot be corrected by
any method. Report the molecule count next to the coverage histogram, never alone.

**Note: Merges chain**, so `OX` can differ from `RX` at two positions while every step was one.

| file | content |
|---|---|
| `<sample>.fq.gz` | reads with `RX` corrected and `OX:Z:` = the original |
| `<sample>.barcodes.tsv` | umi, reads, corrected reads, parent |
| `refine.coverage.tsv` | molecules per power-of-two MIG size, after correction |
| `<sample>.cells.tsv` | cell, molecules, called — OrdMag, only with cell barcodes |
| `<sample>.rank.tsv` | the barcode-rank curve and its CDF, log-spaced ranks |
| `<sample>.bins.tsv` | per MIG size: barcodes, reads, merged as error, payload entropy |

**Never: Cells are called on MOLECULES, never reads** — read depth is amplification. OrdMag
(`--expect-cells`, default 3000) makes the call; the knee is reported beside it, and a disagreement
past 3x is warned. EmptyDrops-style rescue is Cell Ranger's job and is not reproduced.

**Note: `fraction_erroneous` by MIG size is the diagnostic to read.** ~94% of singleton barcodes are
error children, ~0.2% at 2-3 reads. A flat curve, or one rising at high counts, means correction is
merging real molecules.

## assemble

**Never: A molecule is sample + cell + UMI, never the UMI alone.** UMIs repeat across cells and samples
by design. The sort key is `(cell, umi, src_index)`; the partition is on the cell when there is one.

**Never: The RT floor is added, not compared**: `Q = −10 log10(p_cons + p_floor)`. An error made before
amplification is in every read, so the two failure modes are independent. Default `--rt-error 1e-4`
(X2), which caps every emitted quality at ~Q38.

**Never: `--contig` assembles one molecule's fragments and stops.** Full-length receptor assembly,
doublet calling and contaminating-chain filtering are **arda's** job, downstream. Contig mode
places reads by exact seeds, cuts them into overlap components, and never extends a component
across a gap.

**Note: Contig assembly needs an unsaturated barcode.** Two fragments of two *different* molecules on
one barcode share no sequence — indistinguishable from two fragments of one. Read
`expected_molecules_per_group` in `assemble.json`: above 1 means a short UMI cannot tag every input
molecule distinctly, which is a design choice, and the contigs are not trustworthy.

**Never: 1-3 reads per UMI is normal, not broken.** Bulk repertoire and shallow 3' GEX both look like
this. `--min-reads` defaults to 1 and nothing is thresholded away. Do not quote these numbers on
such a library: the split threshold is inert (needs ~30 reads in a group), and the count-ratio
error-child null has no dynamic range. What is reportable is the coverage histogram and the fact
that the UMI is buying *counting*, not error correction.

**Note: A tie is not an `N`.** It is resolved by base order and the emitted quality says so (~Q3). An
`N` would throw away the fact that it is one of two.

| file | content |
|---|---|
| `<sample>.consensus.fq.gz` | one record per molecule, barcodes in the name and in RX/CB/BC/MI/cD |
| `<sample>.mig.tsv` | cell, umi, contig, molecule, reads, length, quality, consensus error, linkage |
| `assemble.coverage.tsv` | groups per power-of-two MIG size |

## Interpreting the statistics

- `mean_reads_per_umi` — over-sequencing. Below ~5, most molecules are seen once and consensus has
  nothing to work with.
- `effective_length` — what the UMI is *worth*, `-Σ log₄ m_j` where `m_j = Σ_a p_j(a)²`. A 12 nt UMI
  with eight fixed positions has effective length 4 and collides constantly. **The nominal length
  tells you nothing on its own.**
- `information_bits` per position — the logo letter height; the bits the UMI is wasting.
- `saturated` — observed UMIs are a large fraction of the usable space; molecule counts are biased
  low and the collision correction is declined rather than guessed.

**Never: Never compute a collision rate from Shannon entropy.** The birthday functional is Rényi-2,
`Π_j Σ_a p_j(a)²`. Since H₂ ≤ H₁, Shannon overstates the usable space and understates collisions —
the direction that silently merges distinct molecules. Shannon is for the logo; the collision form
is for every decision.

### The barcode space and the birthday problem (`checkout.barcode_space.tsv`)

- `nominal_space` — `4^L` over the *captured* positions. `NNNNtNNNNtNNNN` captures 12, not 14: the
  `t`s are scored pattern positions, not barcode.
- `effective_space`, `bias_loss` — what the observed base composition supports, and how much the
  synthesiser mix cost. A real oligo "N" is not 25/25/25/25.
- `occupancy`, `lambda`, `molecules` — molecules land independently, so occupancy is Poisson and
  what you see is the *occupied* count: `occupied = S(1 − e^−λ)`, `molecules = S·λ`.
- **`p_multi`** — `P(k>1 | k≥1)`, the fraction of MIGs that are two or more molecules pooled. This
  is the number to read. Their consensus is a mixture of templates, and over-sequencing cannot fix
  it. `checkout` warns past 5%.
- `saturated` — past 90% occupancy the estimate is declined, because `S` is inferred from the
  observed barcodes and the inversion would report "no collisions" for the most collided library
  there can be.

Note: `Π_j m_j` assumes the positions are independent, so it is a *lower* bound on the collision
probability. Measured against a model-free count on real data (`scripts/collision_check.py`),
collisions ran **1.86×** the prediction — but a permutation puts the position-independence part of
that at only **~1.01×** (`scripts/permutation_nulls.py`), so the rest is the read threshold: a
collided barcode carries two molecules' reads and is over-represented among the MIGs big enough for
a split to be visible.

### The error budget (same file)

`err_estimated` comes from the excess of barcode pairs at Hamming distance 1. `err_predicted` is
`⟨10^(−Q/10)⟩ + ε_pol·cycles`.

**Note: The Phred term is the mean of the probabilities, not `10^(−mean Q/10)`.** The function is
convex, so the low-Q tail carries nearly all the error: half at Q40 and half at Q10 is 5%, not the
0.3% "mean Q25" suggests.

**Never: The distance-1 estimator has a working range, and fails downward.** It subtracts the
coincidence expectation from the observed pair count; once a barcode's `3L` neighbours are mostly
real barcodes that is a small difference of two large numbers. Against an injected 3·10⁻³ it
recovered 0.92× at 0.3% occupancy, 0.65× at 16%, 0.23× at 50%, 0.001× at 93%. `err_unreliable` is
set past 5% neighbourhood occupancy — **do not quote the estimate when it is set**, quote the
prediction and say the barcode is too short.

Where it is set, `scripts/permutation_nulls.py` does better: it replaces the derived coincidence
expectation with a **column shuffle** of the observed barcodes (every marginal kept, every error
child destroyed). On the HIV library that background is 92% of all distance-1 pairs, and the
resulting estimate is 3.4·10⁻³ — within 1.7× of the Phred + polymerase prediction, where the
analytic estimate sits 2.6× *below* it.

## Splitting a MIG into two consensuses

**Never: The threshold is a Bonferroni'd `-log10 p` of 8.68, not the nominal 2.00.** Reads are not
exchangeable: a low-quality read carries a minor base at many positions at once, which is
indistinguishable from a linked subclone if you only look at the columns. The false-positive curve
comes from randomising the reads × positions minor-allele matrix while preserving **both** margins
(per-position error count *and* per-read error load — a curveball swap). The nominal threshold
calls 30.62% of MIGs as two molecules; the measured one calls 1.60%. `docs/nulls.rst`.

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
- **`normalised` is not an error count.** It is reads whose tag was found on the other mate or the
  other strand, so the pair was swapped or the read flipped. An amplicon library sequenced both
  ways otherwise loses half of every MIG at consensus, silently.
- **Output gzip is level 1, not 6.** zlib does 7 MB/s on random DNA at level 6 and 137 MB/s at
  level 1 for 13% more bytes. The file is a little larger on purpose.
- **`-t` never changes the output.** If it appears to, that is a bug, not a tuning question.

## Speed and memory

Reported on every run. ~1.18 M reads/s at 16 threads; matching *and* compression run on the
workers, and the serial stage only appends bytes.

**Note: The UMI counters are the allocation that scales with the library** — ~22 bytes per distinct
UMI (a sorted `(key, count)` array; a hash map is ~48). They are **not yet partitioned**, so a
NovaSeq-scale run holds ~8.8 GB in one piece. checkout warns past 1 GB. The fix is the range
partition with `.mig` bucket output (M2), not a smaller struct.

## Comparing grouping accuracy

`scripts/compare_calib.py` scores read partitions against a truth TSV by adjusted Rand index, and
reports **splitting and merging separately** — splitting inflates the molecule count and is
recoverable, merging mixes molecules and destroys real variants.

migec today groups on the barcode alone; Calib clusters on barcode *and* sequence. The gap is the
collision rate: ARI 1.0000 on a clean 12 nt barcode, 0.8877 at 6 nt with 40% of reads merged.
`effective_length` predicts it before any clustering runs.

## Validating a pipeline

`scripts/spikein_ratio.py` computes the published MIGEC spike-in metric: a real variant against the
worst *error* at the same substitution distance. Targets after UMI consensus: `V1/Err1` 26.5–75.9,
`V2/Err2` 4.6–6.2, against a raw-read baseline of ~1.4 and ~0.3.

Note: Anchor on the junction's 3′ end only. V1 differs at position 4 and V2 at 7–8, so a 5′ anchor
makes both variants count as zero and the metric look perfect.

## References in the repo

- `docs/checkout.rst`, `docs/umi_statistics.rst`, `docs/performance.rst`, `docs/grouping.rst`,
  `docs/validation.rst`, `docs/formats.rst`, `docs/nulls.rst`
- `CLAUDE.md` — the non-negotiables and why each exists
- `project/` — the design record: six subsystem designs and two critiques, with 25 corrected errors
