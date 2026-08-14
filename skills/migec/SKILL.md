---
name: migec
description: Use when working with UMI-tagged or barcoded sequencing data — extracting sample/cell/UMI barcodes from reads, demultiplexing with degenerate adapter patterns, trimming barcodes and transferring them to FASTQ headers, computing UMI coverage histograms and base-composition/entropy statistics, correcting barcode errors, or assembling molecular consensus sequences. Covers the migec CLI and Python API, the .mig format, and the barcode-pattern grammar.
license: GPL-3.0-or-later
---

# migec

UMI barcode extraction, correction and consensus assembly. C++20 core, Python CLI. Successor to
MIGEC (Groovy) and MAGERI (Java); the algorithms carry over, the code does not.

**Status.** All five pipeline commands work: `checkout`, `suggest`, `refine`, `assemble`,
`subsample` — plus `sheet`, `info` and `plot`, which read no reads and write no pipeline output.

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

## Assay profiles: what the experiment implies

A preset places the barcode. `migec sheet --assay NAME` says what a consensus is worth, which is
the axis that decides whether the run is usable. `--assay all` prints every one.

| assay | sensitivity | `--min-reads` | also | layout |
|---|---|---|---|---|
| `airr` (`repseq`, `quasispecies`) | counting | 1 | | `migec`, `primerid` |
| `amplicon` (`targeted`, `panel`) | sensitive | 2 | | `umi` |
| `exome` (`capture`) | sensitive | 2 | | `umi` |
| `ctdna` (`cfdna`) | ultrasensitive | 3 | `--pre-amp-error 7.37e-5` | `tso500`, `duplex` |
| `mrd` | ultrasensitive | 3 | | `migec`, `duplex` |
| `rnaseq` | counting | 1 | `--fast` | `smarter-umi`, `umi` |
| `10x-gex` (`10x`, `gex`) | counting | 1 | `--fast` | `10x` |
| `10x-vdj` (`vdj`) | counting | 1 | `--contig` | `10x-v2`, `10x` |

**Never: `amplicon` is not an alias for `airr`.** A targeted panel of a few PCR-amplified regions
is also an amplicon assay and wants the opposite settings. Two profiles; `amplicon` is the targeted
one, because that is what the word means outside immunology.

**Never: there is no reverse transcriptase in a DNA assay**, and four of the eight profiles are
DNA. `--rt-error` is now spelled `--pre-amp-error` (the old name still works). On an RNA library
the floor is an RT miscall plus the first PCR cycle; on a DNA library it is library-prep damage
plus the first PCR cycle — guanine oxidation during acoustic shearing gives 8-oxoG and therefore
`G>T`/`C>A`, cytosine deamination gives `C>T`/`G>A`. Both predate amplification, so every read of
the group carries them and no consensus removes them. Note: that damage signature is **not** the
`-> G` artifact measured here, which is the 2-colour dark-G instrument artifact. A `C>A` excess
instead means oxidative damage, and `--min-reads` will not fix it — enzymatic repair before
ligation will, in the wet lab.

**Never: do not carry a threshold across the axis.** `--min-reads 3` on a shallow repertoire
library discards 79% of the barcodes and nothing downstream can tell a filtered molecule from an
absent one. `--min-reads 1` on a rare-variant assay feeds raw read quality to the caller: measured
on certified cfDNA, the 2-colour dark-G artifact was **additive to true positives**, so the 0.25%
arm read 0.79% — a quantitative answer wrong threefold that looks perfectly well-formed.

Note: `payload_diverse` in `sheet.ASSAYS` records whether the template itself distinguishes
molecules. On AIRR and MRD the rearrangement is near-unique, so `refine`'s payload-agreement term
and `assemble`'s linkage sub-clustering have real evidence; on exome, RNA-seq and a ctDNA panel
every molecule at a locus reads the same and the barcode carries the whole burden. Never: payload
agreement must be discounted by the measured clonality — on a clonal library it is worth nothing.

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
migec checkout ... --mig                                # .mig buckets, not FASTQ: assemble skips
                                                       # its partition pass over them
migec assemble ref/S1.fq.gz -o cons/                   # one consensus per molecule
migec assemble out/S1.000.mig -o cons/                 # ...or the buckets --mig wrote
migec assemble out/S1.fq.gz -o cons/ --contig          # random-primed reads tiling a molecule
migec assemble out/S1.fq.gz -o cons/ --rt-error medium # the floor by chemistry: rt|medium|high
migec assemble out/S1.fq.gz -o cons/ --fast            # counting mode: modal sequence, max quality
migec plot cons/                                       # QC figures from the tables just written
migec suggest cons/S1.consensus.fq.gz -o qc/           # what synthetic sequence is still in there?
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

### Variant calling (`docs/variants.rst`)

Same rule, one level on: a caller that *transports* the barcode composes, one that *deduplicates*
on it replaces a stage.

| caller | after `assemble`? | why |
|---|---|---|
| Mutect2, LoFreq, FreeBayes, VarDict, bcftools | **yes** | they read a BAM and ignore `RX`; depth already is a molecule count |
| UMI-VarCal | **no** | own UMI pileup and consensus -- an alternative to `assemble` |
| UMIErrorCorrect | **no** | aligns, then groups on (position, UMI) -- an alternative pipeline |
| DREAMS-vc, Shearwater | **no** as normally run | per-position error model fitted on *reads*, against a panel of normals |

**Never: do not apply a UMI-aware caller's family-size filter to a consensus.** After `assemble`
every family has size 1 by construction, so `--min-family-size 3` discards the entire library and
reports zero variants without erroring.

Which caller matters less than the molecule count it is given:

```
molecules at a site = input DNA / 3.3 pg x strands recovered x efficiency
variant molecules   = molecules at a site x VAF
```

**Note: depth buys molecules until it does not.** Deeper sequencing recovers more of the molecules
that are in the tube -- measured 6,310 / 10,299 / 16,809 per amplicon at 3.3 / 10 / 30x on the same
20 ng library -- but the ceiling is the number of input molecules, and past it further reads only
raise reads-per-molecule. Both depth and input are worth spending on until the molecule count stops
rising. The supporting count is then a Poisson draw: an expectation of 3 means a third of
replicates see fewer than 3, whatever the caller.

**Never: for a multiplex panel divide by the panel size.** A variant sits on one amplicon, so the
library total overstates the evidence by exactly that factor. Measured: 5 amplicons for
`PRJNA788522`, 3 for `PRJNA507366`. At 0.125% VAF this is what separates a certain call from a coin
flip -- 20 ng gives 15.0 variant molecules (P(>=3) = 1.00), 5 ng gives 3.4 (P = 0.65).

Defaults, from Maruzani et al. 2024 (BMC Genomics, doi:10.1186/s12864-024-10737-w): **LoFreq** for
a balanced call set, **Mutect2** for sensitivity (and run `FilterMutectCalls`, which that
comparison deliberately did not), **UMI-VarCal** for the fewest false positives on UMI data,
**bcftools** never below ~8% VAF.

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
| `<sample>.sizes.tsv` | the MIG size spectrum at EXACT sizes: molecules and their reads |
| `<sample>.umi_errors.tsv` | per parent depth: error children, their reads, the rate implied |

**Note: the barcode error rate is reported twice, and that is the point.** `estimated_error` inverts
the excess of distance-1 NEIGHBOURS, of which a barcode has only `3L`, so it saturates and fails
downward as the space fills. `error_at_depth` divides the reads in a parent's error children by the
`c*L` barcode bases that parent had to miscall; reads have no ceiling, so it does not saturate.
`error_phred` is the same number as a Phred, for comparison with the barcode's own reported Q.
Against a known injected rate, as a fraction of truth: 0.99 / 0.97 at 0.2% occupancy, 0.88 / 0.76
at 9.8%, 0.62 / 0.45 at 33%.

**Never: BOTH are bounded by the merges correction made**, so on a FULL barcode space both fall to
zero -- `correct_umis` refuses to merge there, rightly, because a distance-1 neighbour is more
likely a real molecule than a child. Read the `saturated` flag; it is what says the answer is a
floor. And read the table at DEPTH: a child whose parent was never sequenced cannot be counted,
which at 1-3 reads/UMI is 80% of them, so `error_from_children` over all depths is a lower bound.

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
(X2, and 10x's own figure for the V(D)J RT), which caps every emitted quality at Q40. `--rt-error`
takes a class rather than a guess: `rt` 1e-4 (default), `medium` 1e-5, `high` 1e-6, or the rate.

**Never: The floor is the ONE-MOLECULE floor.** 10x give Q40 to a base covered by a single UMI and
Q60 only to one covered by two or more, because an RT error is common-mode within a molecule and
independent between them. Every record migec emits is one molecule, so it never claims the
two-UMI number; combining molecules is arda's job.

**Note: `--fast` is counting mode**, not a faster consensus: the modal exact sequence per group, each
base carrying the best quality any read of that sequence reported, no column model and no
sub-clustering. Use it when the deliverable is molecule counts. Refused with `--contig`, whose
tiling reads share no exact sequence to vote on. `support` in the per-molecule table says how many
of the group's reads carried what was emitted.

**Note: Coverage into the consensus is capped at 10,000 reads per barcode** (10x's rule). Never: the cap
is on the reads CONSENSED, never on the reads COUNTED — `cD` stays the molecule's true depth.

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

Reported on every run, and every stage is byte-identical at any `-t`. Measured on this machine,
1 thread against 16:

| stage | 1 thread | 16 threads | scales with |
|---|---|---|---|
| `checkout` | 213,880 | 1,548,835 | reads |
| `refine` | 617,802 | 1,554,156 | distinct barcodes |
| `assemble` | 554,106 | 2,470,928 | reads, then the largest bucket |

reads/s. `assets/benchmark_threads.tsv` is the committed checkout table and the figure is drawn
from it; regenerate with `python scripts/benchmark_threads.py --reads 2000000 -o assets/`.

**Note: profile the thread helper, not only the work.** With all three stages threaded, the largest
single cost in the pipeline was `parallel_for` handing out one item per atomic `fetch_add` — 21% of
all CPU samples across every thread, on one instruction, because sixteen cores were serialising on
one cache line. Items are claimed in batches now. Anything that spreads per-read work should be
checked against the same trap.

**Note: `checkout --mig` writes the partition `assemble` builds.** `<sample>.<bbb>.mig` buckets
instead of per-sample FASTQ, on the same key and the same range partition, so
`migec assemble out/S1.000.mig -o asm` skips its first pass -- 1.16 s -> 0.98 s end to end on 500 k
reads over four samples, identical molecules. One bucket file names the whole partition. Opt-in:
FASTQ is the default because a `.mig` file is a migec intermediate nothing else reads. Never: a
directory of two samples' buckets is refused (a UMI repeats across samples), and `--limit-*` on a
partitioned input is refused (a limit is a prefix; a partition has none). `refine` cannot read
buckets yet, so `--mig` goes checkout -> assemble.

**Note: The UMI counters scale with the library, and bound themselves past a budget** — ~22 bytes
per distinct UMI (a sorted `(key, count)` array; a hash map is ~48), which is ~8.8 GB at NovaSeq
scale. Past `umi_budget_bytes` (a `checkout.run()` kwarg, 1 GB per run, `0` disables) they
range-partition into `<out_dir>/.umi_spill`, which is removed once the summary is written, and every
statistic over them streams a bucket at a time. `umi_spilled` in the summary says whether it fired;
it costs ~2.2x the wall clock when it does. Never: **correction runs twice, the second pass on keys
rotated by the width of the partitioned prefix.** A single pass would bound the memory and silently
stop correcting every error that landed in the prefix — a third of them — while reporting the
smaller merge count as if the library were cleaner.

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

## Limit of detection (`docs/detection.rst`, `scripts/detection_limit.py`)

Two numbers answer "how low can this assay go", and the caller is neither: **N** molecules covering
the site, and **p** the per-MOLECULE error floor. Every assay is in one of two regimes:

| | molecule-limited | floor-limited |
|---|---|---|
| fix | more input DNA, or track more sites | lower floor: proofreading enzyme, or duplex |
| does **not** help | deeper sequencing, a better caller | deeper sequencing, **more input DNA**, a better caller |

**Never: the crossover is `VAF = p/3`** -- the frequency at which a true variant molecule is as rare
as the chemistry's own false ones. At the `rt` floor of 1e-4 that is **3.3e-5**, and no amount of
input reaches below it. A 50 ng / 30-site MRD panel has molecules for 6.9e-6 and a single-strand
floor at 3.3e-5: the molecules promise what the chemistry cannot deliver.

```bash
python scripts/detection_limit.py --input-ng 20 --sites 5
python scripts/detection_limit.py --input-ng 50 --sites 30 --rt-error duplex
```

**Note: MRD pools evidence across sites.** Tracking 30 patient-specific variants is 30x the
molecules that can carry a signal and 30x lower a reachable VAF, from the same blood draw. This is
migec's original application -- the leukaemic clone's IGH rearrangement -- with the clonotype half
belonging to **arda**, whose AIRR `duplicate_count` is then a molecule count.

**Never: measured on a real panel, a library total is not on-target evidence.** Off-target product
took **5-7% at 80 ng, 24% at 20 ng and 47-58% at 5 ng** -- its absolute count barely moves while
on-target scales with input, so precisely when input is scarce the total is most misleading. And
coverage is not uniform: the weakest target held 0.31-0.61x of the on-target mean. Count molecules
**per target** (distinct `MI` in the region) and quote the weakest, not the mean.

## Running a cohort

Two, both in `docs/nextflow.rst`:

- **`integrations/slurm/`** — `migec_sample.sbatch` (env-driven, one sample) and
  `migec_array.sbatch` (one array task per sheet row). Both run as **ordinary bash without SLURM**,
  which is how they are tested and how a layout should be checked before queueing.
  Note: array task 1 is the first *data* row -- the header is skipped, not counted, so the range is
  `1-(rows - 1)`.
- **`integrations/nextflow/`** — `--mode consensus | ctdna | airr`. The `ctdna` mode aligns and
  calls variants, `airr` calls clonotypes with arda. The align module exits non-zero if `MI:Z:`
  does not reach the BAM, because that failure is otherwise silent. Note: nextflow is not installed
  on this machine, so those modules are spec-reviewed rather than run.

Sizing: `checkout` wants cores (scales with reads), `refine` wants memory (scales with **distinct
barcodes**; `table_bytes` in its JSON sizes the next run), `assemble` wants cores (peak memory is
set by the bucket count, not the library).

## Fetching public data

`scripts/sra_fetch.py` — `probe` (read structure from metadata, no download), `url`, `peek`
(first N spots, enough for `migec suggest`), `get`. Note: NCBI's S3 mirror is the default because
it is **33x faster than ENA's ready-made FASTQ** here (6.7 MB/s on 8 connections against 200 kB/s);
`--prefer ena` exists for studies that deposited a file the `.sra` folds away. Never: some runs
list their *original submitted* file ahead of the `.sra`, so select on the file `type` rather than
taking the first entry.

## References in the repo

- `docs/checkout.rst`, `docs/umi_statistics.rst`, `docs/performance.rst`, `docs/grouping.rst`,
  `docs/validation.rst`, `docs/formats.rst`, `docs/nulls.rst`, `docs/variants.rst`,
  `docs/nextflow.rst`
- `notebooks/` — six marimo notebooks, four of which need no network; `notebooks/README.md`
- `CLAUDE.md` — the non-negotiables and why each exists
- `project/` — the design record: six subsystem designs and two critiques, with 25 corrected errors
