# Changelog

Hand-written and prose-heavy: each entry says what changed and, where it matters, which failure it
prevents. Releases before 2.0.0 are the Groovy MIGEC and are described by their git tags on the
`legacy-v1` branch.

## Unreleased

### MIGEC 1.2.9, head to head

`scripts/compare_migec_v1.py` runs the Groovy implementation this repo replaced and migec 2 over
the same simulated library, same barcode dialect, same sheet, both pipelines end to end.
**9.4-11.9x the wall clock** against an M5 gate of 3x, 3.5-3.7x less memory, and the molecule count
is 0.09% over truth against v1's 13.6% -- 99.8-99.99% of our consensuses are exactly a template,
against 93.9-95.1%. The over-count is barcode errors `--filter-collisions` did not catch, because
its rule is a count ratio and a count ratio carries nothing below ~3 reads per UMI.

`--min-count` is given to both: v1 defaults to 5 and migec 2 to 1, so leaving each at its own
default would have credited us with recovering molecules v1 was told to throw away.
`docs/validation.rst`, `assets/migec_v1.tsv`.

### UMI-tools and fgbio, measured

`scripts/compare_grouping.py` runs both map-first tools end to end and scores all three partitions
against the simulator's truth with the same adjusted Rand index Calib is scored with. What it
measures is what the alignment position is worth, because that is the only difference: they group
on *(position, UMI)*, we group on *(sample, cell, UMI)* and align once, afterwards.

On **one reference** — a single amplicon, a clonal control, a targeted ctDNA panel — the position
carries nothing, and migec wins: ARI 0.9967 against 0.9864 (UMI-tools) and 0.9817 (fgbio), with
0.65% of reads in clusters that mix molecules against 3.0% and 3.9%. That is 4.6x and 6x fewer
molecules destroyed, and destroying one is the error nothing downstream can detect. Nothing
rescues that case: two molecules that collide on one reference hold the same sequence, so neither
the mapping position nor payload sub-clustering can tell them apart, and the barcode is the only
evidence there is.

On 200 or 20,000 distinct references they win by 0.001 ARI. That is where the comparison stops
rather than a limit of the method: collided molecules there carry *different* sequences, which is
what makes them separable, and `assemble`'s linkage sub-clustering separates them from the payload
with no aligner. The column is scored at `refine`, one stage earlier. migec is 8-48x faster
throughout, including the alignment they cannot skip.

`docs/grouping.rst` has the table, `assets/grouping_tools.tsv` the numbers, `SOURCES.md` how to
install both (fgbio needs a JDK 17+; UMI-tools needs `--no-build-isolation` on Python 3.12+). The
simulator now writes `clones.fa`, the reference a map-first tool has to align to.

### `refine`'s barcode table bounds itself

The last thing in the pipeline that scaled with the library. Past 1 GB the table range-partitions
to disk, correction follows it into the partition, and every table refine writes is streamed one
bucket at a time. Nothing changes but the wall clock: a partitioned run and a resident one agree on
every scalar and on every output file byte for byte.

Two things this needed that `checkout`'s counters did not. The table carries the **evidence** —
the barcode's own quality at each position and its payload draft — because a side array indexed
against the entry list cannot survive a partition, and dropping it would leave the bucketed run
correcting on the count ratio alone, which reports nothing at 1–3 reads per UMI. And the two passes
**scan**; they do not merge. A barcode can have a plausible parent on each side of the partition
boundary, and merging inside a pass takes the first candidate rather than the best: 2 barcodes in
6,591 landed on a different parent than the resident run gave them, and every table downstream
moved with them. Both passes now propose and one global apply decides, in the same order the
resident walk uses.

The budget is a `refine.run()` argument, never a CLI flag. Two smaller things came with it: the
evidence pass folded into the table pass, so refine streams the reads twice rather than three
times, and the barcode-rank curve is read off the MIG size spectrum rather than off a sorted array
of every molecule's count — the last allocation in the stage that was one entry per barcode.

### `--mate2`: both ends of the fragment are one molecule

A read pair carries one barcode because it is one molecule, so `assemble` will take the other mate
and place it:

```bash
migec assemble out/S1_R1.fq.gz --mate2 out/S1_R2.fq.gz -o cons/
migec assemble out/S1.000.mig --merge-mates -o cons/     # the pair is already in the record
```

Mate 2 is reverse-complemented and placed against mate 1. Overlapping mates give one consensus
spanning the insert; mates that do not reach each other give two contigs, because the bases between
them are covered by no read. That is `--contig`'s rule rather than a second one, and the two output
routes agree byte for byte.

The offset is a property of the molecule, not of the pair, so it is voted on once per group over
up to eight pairs instead of placing every read against every other. The first version did the
latter and cost 11x the single-end path (119,820 record-pairs/s against 1,356,819); the vote costs
1,288,686 against 2,115,912, which is the price of the second read's bases and nothing else.

The mates are matched by position, which is what checkout's `_R1`/`_R2` output guarantees. A file
that ends before the other is refused: pairing off what is left attaches one molecule's mate to
another's, and the consensus that comes out looks perfectly well formed.

### `--pre-amp-error auto` fits the floor instead of naming it

The cap on every emitted quality is the error that was already in the molecule before
amplification, and until now it was a named bracket: `rt` 1e-4, `medium` 1e-5, `high` 1e-6. It can
be measured on the dataset in front of you, and `auto` does it -- X2's estimator, run on migec's
own consensuses rather than on a bespoke one:

```bash
migec assemble out/S1.fq.gz -o cons/ --pre-amp-error auto
```

Molecules at 20 reads or more have suppressed sequencing error to nothing, so what they still get
wrong against the library's modal sequence is what was there before the first cycle. Injected
1e-4 comes back as 1.20e-4 [7.94e-5, 1.75e-4] and injected 1e-5 as 1.58e-5 [8.14e-6, 2.76e-5]; with
no injected error at all the answer is the interval's upper end, 9.56e-6, because a floor of zero
is a quality of infinity.

It refuses more often than it answers, and that is the point. A diverse library has no monomorphic
position to score, so "disagrees with the modal base" means "is a different molecule" and the
number that would come out is the library's diversity; a shallow library has no molecule deep
enough for the curve to have flattened. Both refuse, fall back to the named `rt` class, and say
why -- in the report and in the new `assemble.pre_amp_error.tsv`, which carries every input the
fit rested on next to the fit.

It costs a second assembly pass, so it is opt-in: the consensus sequences do not depend on the
floor, only the emitted quality does, so a probe assembly runs, is read back and is deleted.

## 2.3.0 — 2026-08-14

**The intermediate is a partition now, end to end**, and the `.mig` format is at v2 (a reader takes
v1 and v2, so an intermediate written by 2.2.1 still reads). Nothing in `checkout` scales with the
library any more, `refine` and `assemble` both take buckets, and an adversarial audit of all of it
found three bugs that destroyed or invented data -- fixed, each with the test that reproduces it.

Upgrading: nothing changes unless you pass `--mig`. FASTQ is still the default at every stage, the
CLI has one new flag, and every number a 2.2.1 run produced is the number a 2.3.0 run produces.

### `checkout --mig` writes the partition `assemble` was building for itself

`assemble`'s first pass range-partitions the reads on the barcode into `.mig` buckets, and that
pass is most of its wall clock. `checkout --mig` writes those buckets directly -- same key, same
range partition, `<sample>.<bbb>.mig` -- and `assemble` reads them instead of building them:

```bash
migec checkout reads.fq.gz -b barcodes.txt -o out --mig
migec assemble out/S1.000.mig -o asm      # one bucket names the whole partition
```

Measured on 500,000 reads over four samples at `-t 4`: **1.16 s to 0.98 s** end to end for the
identical 124,878 molecules, and the consensus FASTQ is byte-identical after decompression (the
gzip framing differs, because bucket boundaries are member boundaries and there are a different
number of buckets). `checkout` pays 0.06 s of that and `assemble` saves 0.25 s.

Opt-in, and FASTQ stays the default: a `.mig` file is a migec intermediate that nothing else
reads, while every aligner and every pipeline in `docs/downstream.rst` speaks FASTQ.

Never: **`-t` still changes nothing but the clock.** One writer per (sample, bucket), owned by one
thread for the whole run, and every worker walks the chunks in input order -- so ownership decides
who writes a record and never which file it lands in or where.

Note: the open-file budget is for the RUN, not per sample. 256 buckets each on a 96-plex sheet
would be 24,576 open writers; each sample of such a sheet also holds a 96th of the reads, so a
couple of buckets each is proportionate rather than a compromise.

Two things are refused rather than guessed at: buckets from two samples in one `assemble` (a UMI
repeats across samples by design, so grouping them together invents molecules), and
`--limit-read`/`--limit-umi` on an already-partitioned input (a limit is a prefix of the input, and
a partition has no prefix -- the first records of bucket 0 are one corner of the barcode space).

Never: **only buckets `assemble` wrote are deleted after they are consumed.** Pass 2 removes each
bucket as it finishes with it, which is right for its own temporaries and catastrophic for
checkout's output -- the first run through `--mig` ate three of the four samples' partitions, and
a second `assemble` over an eaten one would report a smaller library with nothing to say why.

### The i7 x i5 contingency table

`checkout` reads the index pair out of the instrument's own read header -- the last field of
`1:N:0:ATCACG+CGTGAT` -- for every read, matched or not, and writes
`checkout.index_pairs.tsv`: reads per observed combination, each one's share of its own i7 and i5,
and whether it looks ordered. Index hopping is the one contamination a per-sample yield cannot see,
because a hopped read lands in another real sample and looks like one of its reads; the header is
the only evidence, and it was already in the file.

Never: **a single-indexed run is not estimable, and that is not zero.** With one index there are no
combinations, so nothing can be off-diagonal.

Note: it matters most where it is smallest. At 0.1% hopping a 1% variant in a deeply sequenced
sample contaminates its neighbour at 1e-5, which is the level a rare-variant caller is asked to
believe.

### refine reads and writes buckets, and the format carries the barcode's own quality

`checkout --mig` -> `refine` -> `assemble` runs on buckets end to end. refine's output is
**re-partitioned on the corrected barcode**: a corrected barcode is a different key and a key
decides its bucket, so copying a bucket through unchanged would stop it being a partition, and the
reads whose barcode was corrected across a boundary would be grouped with strangers by the next
stage. Every number matches the FASTQ route -- reads, barcodes, merges, molecules, the estimated
error to nine digits, the barcode table byte for byte, and the consensus FASTQ byte for byte.

`.mig` **v2** adds the BARCODE's own quality, one Phred per base, as two fixed-width columns. v1
stored the minimum over the barcode, which is not what the correction posterior wants: it weighs
the reported quality AT THE POSITION THAT DIFFERS, and a minimum says every position is as bad as
the worst one -- overstating the error everywhere, which makes merges easier, which is the wrong
direction. A v1 file still reads and falls back to the global rate, exactly as a FASTQ with no
`QX:Z:` tag does.

Note: the audit trail on the bucket route is `<sample>.barcodes.tsv`, every barcode with its
parent. A `.mig` record has no room for the pre-correction barcode the way a FASTQ comment has
`OX:Z:`, and two `u64` per READ to carry what is one row per BARCODE is the wrong trade.

### Three data-integrity bugs, found by attacking the new code

All three reproduced before they were fixed, and all three are in the `.mig` work above.

- **refine wrote its output buckets over its input.** Same names on both sides, and `MigWriter`
  opens `"wb"`, which truncates. Pointed at its own input directory the run destroyed the reads it
  was halfway through reading. Refused now, on the resolved path, before a writer is opened.
- **refine's `.mig` rewrite ignored `--limit-read`/`--limit-umi`.** It wrote reads the table was
  never built from -- uncorrected, beside corrected ones, with nothing saying which was which.
- **The sample id went into every output path unvalidated**, and in `assemble` and `refine` it comes
  from the DATA (the `BC:Z:` tag of the first read), so a crafted tag wrote outside the output
  directory. `validate_sample_id` refuses path separators, `..`, a leading `-` and control
  characters.

### ...and six more from the same audit

A truncated spill file was read into an undersized buffer (a 15-byte heap overrun; trigger is a run
killed mid-spill). A `.mig` header's `bucket_bits`/`bucket_index` went unvalidated into a shift and
a vector index. A block's `n_records` was `memcpy`'d before it was bounded by the payload. `assemble`
left its whole temp partition behind on the error path and leaked the per-bucket table handle on a
throw. `fclose` was unchecked in four writers, so a disk that filled during the final flush was
reported as a successful run. A failed thread SPAWN unwound over joinable threads, which is
`std::terminate`.

Two usability bugs went with them: a barcode **sheet** could not use the positional spellings that
work through `--bc-pattern`, and an **empty input** printed statistics computed from no reads and
exited 0.

### Tests and docs

`tests/realworld/` exists now -- it was an empty `__init__.py` -- and runs the two public CI
fixtures through the stages they can support. `tests/unit/test_error_paths.py` covers empty,
truncated and malformed input into every stage. Paired-end flows past `checkout` for the first time
in a test, and the four-stage chain is asserted as one invariant. `.github/workflows/benchmark.yml`
runs the benchmark and realworld tiers nightly and on dispatch, never on the PR path.

README.md went from 809 lines to 206 -- install, where the barcode is, the five commands, what comes
out, and a documentation table -- with every fact moved into a docs page rather than deleted, and
`docs/assays.rst` is a new home for the per-assay settings table that was buried in a 762-line ctDNA
page.

### The UMI counters bound themselves, and correction follows them

This is roadmap items 1 and 2, which turned out to be one item. `checkout`'s UMI counters were the
last allocation in the pipeline that grew with the library rather than with the chunk: ~22 bytes per
distinct barcode held in one piece, 8.8 GB at NovaSeq scale, with a warning past 1 GB standing in
for a fix.

They now range-partition to disk past `umi_budget_bytes` (1 GB for the whole run, divided by the
samples, `0` to keep everything resident). The buckets go in `<out_dir>/.umi_spill` and are removed
once the summary has been written; `umi_spilled` in the summary says whether it happened. The
histogram, the composition, the distinct count, the distance-1 census and the correction all stream
one bucket at a time. Not a CLI flag: 1 GB is where a counter becomes a problem on any machine, and
nothing below it changes.

Never: **a range partition on its own would have bounded the memory and silently stopped
correcting.** The bucket is the top bits of the key, so a barcode whose error landed in the first
few positions sits in a different bucket from its parent and the two can never meet — and a third
of all barcode errors are in those positions. Correction therefore runs twice: once over the
buckets as they stand, owning the positions the prefix does not touch, and once over a copy whose
keys are rotated left by the width of the prefix, owning exactly the ones it hides. Every pair is
weighed in one pass and only one, so `umis_merged` stays a count of barcodes rather than of
opportunities to look at them. Checked field for field against the resident answer, on a simulated
library with injected barcode errors and on a 500,000-read corpus: identical, including the
estimated error rate.

Never: a bucketed run reports the **scalars**. `root` and `corrected` are indexed against
`entries()`, which is the array being bounded, so they come back empty rather than wrong, and
`BarcodeEvidence` — indexed the same way — is refused rather than ignored. The stage that rewrites
reads is `refine`, and `refine` does not spill.

Two bugs found on the way, both by tests written to state the property rather than the number:

- The **first** flush of a counter swaps the append buffer in wholesale and returned before the
  budget check, so a counter whose whole library arrived in one buffer's worth never partitioned at
  all — the bound switched itself off on exactly the small case, and reported a resident answer
  that happened to be right.
- The bucketed driver reported the census's `0.0` for a clean library where the resident path
  reports the 1e-4 floor it actually corrects at. The rate quoted has to be the rate used; it feeds
  the error budget, which divides by it.

Cost: about 2.2x the wall clock when it fires (718,000 reads/s against 333,000 on the benchmark
corpus), because the partition is read back four times. Nothing when it does not fire — 217,000
reads/s single-threaded and 23.7 B per distinct UMI, both unchanged.

## 2.2.1 — 2026-08-14

Documentation, pipelines and tooling. Note: **the wheel is unchanged from 2.2.0** — no command,
flag, output file or number that migec itself produces is different. Everything below is in the
repository: the SLURM templates, the Nextflow downstream layer, two more notebooks, a new docs
page, and the two scripts behind it.

### Pipelines you can actually run

`integrations/slurm/` is new: `migec_sample.sbatch` (one sample, environment-driven) and
`migec_array.sbatch` (one array task per sample-sheet row), with an example sheet and a README that
sizes the request per stage. Both run as **ordinary bash without SLURM** — every SLURM variable has
a fallback — which is how they are tested here and how a layout should be checked before a cohort
is queued. Note: array task 1 is the first *data* row; the header is skipped rather than counted,
so the range is `1-(rows - 1)`, and an `--array=0-N` runs a task that reads the header as a sample.

`integrations/nextflow/` gained a runnable entry point (`--mode consensus | ctdna | airr`) and a
downstream layer: `modules/downstream/{align,callvariants,arda}` plus `subworkflows/migec_ctdna`
and `subworkflows/migec_airr`. The align module **checks its own output** — if no `MI:Z:` tag
reaches the BAM it exits non-zero and names the flag each aligner needs, because that failure is
otherwise silent and surfaces much later as an inexplicably untagged BAM.

### Notebooks

`notebooks/README.md`, plus `exome_capture.py` (why coordinate deduplication undercounts a capture
panel) and `airr_repertoire.py` (clonotype counts from reads against molecules). Both are
simulated, so both run offline against a known answer. Never: the PEP 723 headers were
**under-declared** — altair renders a polars frame through pyarrow and then pandas, so a clean
machine failed three imports deep with an error naming neither. Every header is now what the
notebook actually needs, and all of them were executed to prove it.

### Which variant caller, and the count that decides it

`docs/variants.rst` answers the question people arrive with — *rare variants out of ctDNA, what do
I run after `assemble`?* — and extends `docs/downstream.rst`'s transport-vs-deduplicate rule one
level: `Mutect2`, `LoFreq`, `FreeBayes`, `VarDict` and `bcftools` read a BAM and ignore `RX`, so
they **compose**; `UMI-VarCal`, `UMIErrorCorrect`, `DREAMS-vc` and `Shearwater` group and consense
themselves, so they **replace** a stage. Never: do not apply a UMI-aware caller's family-size
filter to a consensus — after `assemble` every family has size 1 by construction, so
`--min-family-size 3` discards the entire library and reports zero variants without erroring.

The caller ranking is quoted from Maruzani et al. 2024 (BMC Genomics,
[doi:10.1186/s12864-024-10737-w](https://doi.org/10.1186/s12864-024-10737-w)), which is the only
independent comparison of the six. The page says which half of it transfers and why.

### The ctDNA ground truth was found, not built

`ROADMAP.md` item 7 said this had to be constructed because Maruzani's deposited runs carry no UMI.
That is true of those two runs and not of the public record. Screening SRA read structure instead
turned up **`PRJNA788522`** (72 runs, cell-free DNA reference material at certified 0 / 0.125 /
0.25 / 1% VAF, crossed with 5/20/80 ng input and 3.3/10/30x depth, three replicates) and
**`PRJNA507366`** (28 runs, six polymerases on the same material plus 0.031% and 0.0625% VAF).
Both kept a real 12 nt inline UMI, which `migec suggest` recovers from base composition alone.

Note: `PRJNA507366`'s design is in `library_name`, not `sample_alias` — every alias there reads
`SeraCare_Reference_Material`, so reading the alias alone would have called the study undesigned.

### `scripts/sra_fetch.py` — fetch on demand instead of mirroring

`probe` reports reads-per-spot and their lengths from metadata alone, which is the cheap "did the
UMI survive deposition" test; `url` gives the direct S3/HTTPS links with md5 and size; `peek`
streams the first N spots for `migec suggest`; `get` downloads and converts. Never: the
obvious-looking source is the slow one — ENA serves ready-made FASTQ and skips `fasterq-dump`, but
caps at one connection. Measured on `SRR17220895`: **ENA 200 kB/s against NCBI S3 at 6.7 MB/s over
8 connections**, so S3 is the default and `--prefer ena` is there for the studies that deposited a
third file the `.sra` object folds away.

This is what keeps `isalgo/umi_data` small: anything with a public accession is regenerated by a
one-line command recorded in `SOURCES.md`, and only the CI fixtures and the semi-internal data that
has nowhere else to live are mirrored.

### `scripts/ctdna_titration.py` and `notebooks/ctdna_variants.py`

The three stages over the titration, reporting molecules per site against the certified frequency.
Never: **the molecule total of a multiplex panel is not the count at any one site** — a variant sits
on one amplicon, so the total overstates the evidence by exactly the panel size. The amplicon count
is measured from consensus prefixes rather than assumed, and its share threshold has to sit in the
gap below the smallest real amplicon: at 1% the count moved with depth (5 amplicons on the deepest
run, 10 on the shallowest) because a shallow consensus carries more payload error, which deflated
the per-amplicon count on exactly the runs where the evidence was thinnest. The script warns when
the count is not constant across runs, which is the check that the threshold is in the gap.

## 2.2.0 — 2026-08-13

Throughput and one new measurement. No breaking change: every 2.1.0 command, flag and output file
behaves as it did, and the additions are a new table, four new JSON fields and two new plot panels.

### The barcode error rate, measured at every depth

`refine` writes `<sample>.umi_errors.tsv`: one row per exact parent depth, carrying the distinct
error children that parent spawned, the reads in them, and the error rate each implies. A parent
carrying `c` reads offered `c*L` barcode bases to be miscalled, so the same eps falls out two ways:

    distinct children  u(c) = 3L (1 - exp(-c eps / 3))   -> eps = -(3/c) ln(1 - u/3L)
    reads in children  r(c) = c L eps                    -> eps = r / (c L)

A barcode has only `3L` neighbours one substitution away, so the first **saturates**; reads have no
ceiling, so the second does not. Where the two part company on the figure is where this library's
barcode neighbourhood filled up, read off the data instead of predicted. Past saturation the first
is left blank rather than reported as a small number: inverting a full neighbourhood returns "no
errors" for the most error-ridden library there can be.

New report lines and JSON fields `error_at_depth`, `error_phred`, `error_from_children`,
`error_depth`. The number to quote is `error_at_depth` — restricted to parents seen >= 10 times,
where correction is near-complete — and `error_phred` is the same thing as a Phred so it can be put
beside the barcode's own reported Q. On a diverse library sequenced 25 deep with 1e-3 injected, the
distance-1 excess gives 9.73e-04, the children give 9.89e-04 (Q30) and all depths give 9.98e-04.
Three routes to one number agreeing within 3% is what makes any of them believable.

Never: **both estimators are bounded by the merges correction actually made**, so neither is
saturation-free. Against a known injected rate, as a fraction of truth:

| occupancy | 0.2% | 2.3% | 9.8% | 33% | 100% |
|---|---|---|---|---|---|
| distance-1 excess | 0.97 | 0.96 | 0.76 | 0.45 | 0.001 |
| from the children | **0.99** | **0.95** | **0.88** | **0.62** | 0.00 |

The children estimate is the better of the two wherever either works, and at 100% they both go to
zero for the same reason: on a full barcode space `correct_umis` refuses to merge, correctly,
because a distance-1 neighbour there is more likely a real molecule than a child. The `saturated`
flag is what says the answer is a floor — read it, do not read this table instead of it. And read
the table at depth: a child whose parent was never sequenced cannot be counted, which at 1-3
reads/UMI is 80% of them. `tests/synthetic/test_umi_errors.py` holds all of it, `docs/umi_errors.rst`
works through it.

Two new `migec plot` panels over that table, `umi_error_children` and `umi_error_rate`, with the
`3L` ceiling and the reported estimate drawn as reference lines. Note: `neighbours` and `estimate`
are constant down their columns on purpose — the panels draw them, and a figure that needs a value
its own table does not carry will one day disagree with the report.

### `reads in them` is points, never a line

The `mig_size_spectrum` panel drew its reads series `with lines` over a table that has one row per
**exact** size. Past the head of the distribution almost every size holds exactly one molecule, so
reads == size, and the line rendered the `y = x` diagonal as the most prominent feature of the
figure — a tautology drawn as a second mode. Where two or three molecules shared a size the same
line sawtoothed between `size*1` and `size*2`, which is integer quantisation drawn as signal, and it
bridged gaps in the support where no size was observed at all. It is points now, and the new
barcode-error panels are points for the same reason.

### The thread helper was the bottleneck

With all three stages threaded, a sampling profile of
`assemble` put **21% of all CPU samples, across every thread, on a single instruction**: the atomic
`fetch_add` in `parallel_for` that handed out one item at a time. When an item is one read's tag
scan, sixteen cores serialise on one cache line and the counter costs more than the work it
distributes. Items are claimed in batches now, sized so each worker takes ~8 turns and capped so a
ten-million-item scan still hands out ten thousand batches. Never: the batch collapses to 1 when
there are few items, which is the uneven case -- one bucket per item, one holding ten molecules and
the next ten million -- that the counter exists for.

Two serial blocks fell out with it, both read-only `3L`-binary-search scans of the barcode table:
the **distance-1 census** in `estimate_umi_error` (checkout's per-sample statistics tail) and
refine's **residual-FDR scan**, which was 0.53 s of a 2.17 s run on one core after everything around
it had been parallelised. Both tally integers into per-worker counters summed afterwards, so `-t`
still changes nothing but the clock -- verified byte-for-byte at `-t` 1/3/16 on a clean library and
on one with injected barcode errors (residual 1,294 either way), and under the thread sanitizer
(104 cases, 224,116 assertions, no race).

Two allocation fixes on paths every read takes: refine's rewrite chunk is now **assigned into
rather than cleared** -- the same trade assemble's partition already made, since `clear()` destroys
four `std::string` per record -- and the rewrite's comment buffer and unpacked barcode are worker
scratch reused read after read, with the cell and UMI halves taken as `string_view` instead of
`substr`. That is four allocations a read removed; `rewrite_seconds` fell 0.51 s to 0.38 s.
`assemble` also unpacked each group's UMI twice, once for the composition tally and once for the
record name.

| stage | 1 thread | 16 threads | was, at 16 |
|---|---|---|---|
| `checkout` | 213,880 | **1,548,835** | 1,056,472 |
| `refine` | 617,802 | **1,554,156** | 1,012,368 |
| `assemble` | 554,106 | **2,470,928** | 2,051,937 |

reads/s, 500 k-read sample except checkout at 2 M; `assets/benchmark_threads.tsv` is regenerated
and the figure redrawn from it. checkout's gap between end-to-end and matching throughput fell from
20% to 9%.

Note: re-measured while it was in hand, a **64 k chunk in assemble's partition is still 32% faster**
(3,075,506 reads/s against 2,324,403 on 4 M reads) and is still not taken: it costs 16 MB of
resident chunk, which at NovaSeq scale makes pass 1 the memory peak and breaks the property that a
finer partition costs less rather than more. The benchmark tier does not fail at 64 k on a 500 k
corpus -- the objection is one of scale -- so the number is written down next to the constant.

## 2.1.0 — 2026-08-13

First non-alpha of the rewrite. The three stages, the eight commands and the on-disk formats have
been stable across a3 and a4; what this release adds is the QC layer catching up with them.

**Four figures you already know how to read.** They come off new tables rather than out of a new
computation, so each is still redrawable from a committed TSV:

- **The barcode rank plot**, on [Cell Ranger's](https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-ab-barcode-rank-plot)
  axes, because it is the figure every user of a droplet protocol has seen: barcodes sorted by
  content, log-log, knee where cells stop. Never: the y axis is **unique UMIs, never reads** — one
  over-amplified molecule would otherwise put an empty droplet high on the curve, which is the
  exact artefact the plot exists to show. The call is drawn on the curve, not described in a
  caption. New table `<sample>.cell_rank.tsv`.
- **The MIG size spectrum**, molecules *and* the reads they account for, against `log(1 + size)`.
  Both series, because they peak in different places the moment a library is over-sequenced: most
  molecules are shallow, most reads are in the deep ones, and a figure with only one of them says
  the opposite of a figure with only the other. `log1p`, so a molecule seen once has a place on
  the axis.
- **The rank/Zipf curve**, molecule size against rank on log-log. Never: this is why the new
  `<sample>.sizes.tsv` is written at **exact** sizes and not power-of-two bins — four bins make
  four steps, and a straight line cannot be told from a bent one. It costs one row per distinct
  depth, a few thousand on a real library, not one row per molecule.
- **Unique UMIs and reads per sample barcode**, off `checkout.summary.tsv`. The multiplexed
  analogue of the same question.

**Consensus quality is a box, not a thinned scatter.** Emitted quality is discrete and capped at
the RT floor, so at any real depth every molecule sits on one or two integers: a cloud of dots
draws that as a flat line whether the bin holds ten molecules or ten million, and the `every 17`
that kept the SVG small threw away the tails that were the only thing the cloud could have shown.
`assemble` now accumulates the exact joint distribution of (depth bin, rounded Phred) — both are
small integers, so it is 61 counters per bin per bucket — and `assemble.quality_by_depth.tsv`
carries real order statistics over every molecule.

**Publication defaults on every panel.** Transparent background, so one SVG serves a light README,
a dark README and print. One ink colour (`#808080`), which reads on both. **The legend is inside
the plot box**: a key in the margin makes every figure wider than its data and is the first thing a
journal asks you to move. Frame is 760x520 rather than 900x560 wide.

**The pipeline figure is page-shaped.** `rankdir = TB` and three `rank = same` groups pinning each
side tool level with the data it reads; it was a 5:1 strip that filled the README column with air.
Transparent, and `minibwa` is in the downstream box.

Twenty panels now, from sixteen. `docs/formats.rst` documents every QC table and its columns.

## 2.0.0a4 — 2026-08-13

**The partition threads, and the memory estimate behind it was wrong.** With the consensus already
parallel, `assemble`'s *partition* was 2.07 s of a 2.69 s run — 77% of the stage, on one thread.
`gzip -dc` on the same file takes 0.23 s, so five sixths of it was not the inflate: it was the tag
scan, the barcode packing, the record serialisation and the level-1 deflate of each bucket block.
All four run on the workers now, by **ownership rather than locking** — worker *w* owns every
bucket with `bucket % threads == w` for the whole run, so a bucket file has exactly one writer and
no bucket state is shared.

| 4 M reads, `-t 16` | before | after |
|---|---|---|
| wall clock | 2.70 s | **1.95 s** |
| partition | 2.06 s | **1.45 s** |
| reads/s | 1,481,946 | **2,051,937** |
| peak RSS | 1,479 MB | **789 MB** |

Half of that was the reader rather than the threading: **assign into the chunk rather than clearing
it**, because `clear()` destroys four `std::string` per record and the reader ends up spending its
time in malloc instead of inflate. One chunk is held rather than one per worker, so it costs ~2 MB
at any `-t`.

Note: a **bigger** chunk is 22% faster again — 2,510,241 reads/s at 64 k reads, because
`parallel_for` starts and joins its threads per call and 4 M reads at 8 k a chunk pays ~15,000
thread creations. It is not taken: 16 MB of resident chunk makes the *partition* the memory peak on
a finely partitioned shallow library, which breaks the property that a finer partition costs less
rather than more, and `test_shallow_memory_is_still_bounded_by_the_bucket` catches it. The upgrade
path is a persistent worker pool, not a bigger chunk.

**Never: an estimate that nothing checks will be wrong.** The constant deciding how finely to cut
the input said a gzipped FASTQ goes resident at **8x** its on-disk size. Measured, it is **19x** —
a resident record is two heap `std::string` with their allocator headers and rounded-up buckets,
plus three 8-byte keys, not the 180 bytes of payload. Guessing low is the expensive direction,
because it picks too few buckets and pass 2 holds sixteen of them at once. That single wrong
number is where the 1,479 MB came from.

**`subsample` says what it did.** The report now carries the median and the deepest reads per kept
barcode next to the mean, and five kept barcodes with their depths. Never: **in key order, not
first-seen order** — a barcode with 100 reads appears early about 100x more often than a singleton,
so the head of a file is a sample of the deep MIGs and of nothing else, which is the same trap
`subsample` exists to avoid, one level down.

**minibwa is in the downstream contract**, run and counted like the rest: 600/600 records keep
`RX`/`CB`/`MI` through a sorted BAM. Note: the comment flag is **`-y` on `minibwa map`** (the
minimap2 spelling) and **`-C` on the legacy `minibwa mem`** (bwa's), and each rejects the other's
flag with a non-zero exit rather than dropping the tags quietly.

**Map first, or collapse first?** `docs/downstream.rst` now works the question through: what the
chromosomal position buys as extra key bits (and when — a 5 nt TSO500 UMI is 1,024 barcodes), what
it costs (N x the alignment, a mismapping becoming a grouping error, needing a reference at all),
and why `assemble`'s linkage sub-clustering recovers most of it from the payload without an
aligner. With a table separating tools that *transport* a UMI from tools that *deduplicate* on one:
the first compose with migec, the second replace a stage of it.

**The docs navigate.** Twenty pages of flat toctree put every long page title in the header; they
are grouped into seven sections now — Installation, Examples, Layouts, Commands, Downstream,
Method, Reference — with landing pages that say what each page answers. Every command page carries
a subtitle (`assemble -- one consensus per molecule`), and `docs/nextflow.rst` is new.

**Nextflow.** Never: `containsKey`, not `?:`. Groovy's elvis treats `false` as absent, so a
per-sample `contig: false` against a params default of `true` silently meant its opposite — the one
direction a per-sample override exists to make possible. `--limit-read`/`--limit-umi` now reach
`assemble`, and the container tag is one param rather than four literals. Note: nextflow is not
installed on the machine this was measured on, so the modules are reviewed against the nf-core
spec, not verified by a run, and the docs say so.

## 2.0.0a3 — 2026-08-13

**Every stage threads now, and the first fix was not a thread.** `refine` and `assemble` were the
pipeline's bottleneck at ~200 k reads/s against checkout's 1.06 M. The measurement said why: zlib
at its default level 6 spent **1.78 s of refine's 2.14 s run** compressing an intermediate that the
next stage decompresses immediately. Level 1 costs 21% more bytes and gave 3x on its own. Then the
parallelism:

| stage | before | after (16 threads) |
|---|---|---|
| `checkout` | 1,056,472 | 1,056,472 |
| `refine` | 222,017 | **1,012,368** |
| `assemble` | 202,977 | **1,434,573** |

`refine --threads` splits the neighbourhood scan, which is a pure function of the barcode table --
it reads no union-find state -- so it parallelises exactly, and the merges it finds are applied
serially afterwards in the original smallest-first order. The result is identical, not merely
equivalent. `assemble --threads` gives each worker its own bucket; the buckets are independent by
construction because the partition is on the barcode.

**Never: `-t` still changes nothing but the wall clock, on all three stages.** That is why the
bucket count is a constant floor of 16 rather than a function of `--threads` -- if `-t` chose how
finely the input was cut, it would choose the gzip member boundaries too, and two runs would
produce byte-different files holding identical records. Asserted three ways: in C++ per stage
(`tests/cpp/test_parallel_stages.cpp`), at the CLI over a full three-stage chain
(`tests/synthetic/test_thread_invariance.py`), and under the **thread sanitizer**, which reports
no data race across 104 test cases and 224,116 assertions -- with the instrumentation proven to
fire on a deliberate race in the same helper.

**`--limit-read N` and `--limit-umi N`** on every stage: stop after N reads, or after N distinct
barcodes. For getting an answer out of a 400 GB run in a minute. Never: a limited run is not a
sample, and says so in its own report -- the first N reads of a FASTQ are one corner of one
flowcell. `subsample` remains the sampler.

**The nextflow integration is rebuilt as three modules and a subworkflow**
(`integrations/nextflow/`), nf-core layout, with `meta.yml` and a `stub:` for each. One process
meant one resource label for three stages with different shapes and no resume between them: a
failed assemble re-ran the whole demultiplex. `refine` now carries `process_high_memory` because
its memory is set by distinct barcodes and by nothing else, and every stage passes `task.cpus`,
which is safe precisely because the output does not depend on it.

**The pre-amplification error floor is named, not guessed.** `--rt-error` takes a fidelity class:

```bash
migec assemble ... --rt-error rt        # 1e-4, caps at Q40 -- anything with an RT step (default)
migec assemble ... --rt-error medium    # 1e-5, caps at Q50 -- no RT, an ordinary polymerase
migec assemble ... --rt-error high      # 1e-6, caps at Q60 -- no RT, a proofreading polymerase
migec assemble ... --rt-error 7.37e-5   # or the rate, e.g. TruSight Oncology 500 v2
```

Which class applies is a property of the protocol, so it is declared rather than fitted, and every
value is cited in `SOURCES.md`. The default is 1e-4 because that is 10x's stated figure for the
V(D)J RT reaction *and* what X2 measured here independently (1.54e-4 on SRR1763769). Never: it is
the **one-molecule** floor and every record migec emits is one molecule. 10x assign Q40 to bases
covered by a single UMI and Q60 only to bases covered by two or more — an RT error is common-mode
within a molecule and independent between them — and combining molecules is arda's job. A
per-molecule record claiming Q50 is claiming two-UMI confidence on one-UMI evidence.

**`migec assemble --fast`: counting mode.** The group's most frequent *exact* sequence, with each
base carrying the best quality any read of that sequence reported for it. No column model, so no
per-base error correction and no sub-clustering — and the RT floor still caps what it claims. For
expression and clonotype abundance, where the deliverable is a molecule count. Measured against
the full path at 8 reads and 5e-3 per base: the column posterior clears essentially every
sequencing error and the majority string keeps what it carried. Refused with `--contig`, whose
tiling reads share no exact sequence to take a majority over. The per-molecule table gains a
`support` column: how many of the molecule's reads carried what was emitted.

**Coverage into the consensus is capped at 10,000 reads per barcode**, which is 10x's rule and
their reasoning — past that the column posterior has long since saturated while the group still
costs time and memory. Never: the cap applies to the reads that are *consensed*, never to the reads
that are *counted*. `cD` and the table's `reads` column stay the molecule's true depth, because
capping a count would flatten the abundance of exactly the most-amplified molecules.

**`migec plot`: sixteen QC panels, drawn with gnuplot from the tables the stages already wrote.**
UMI PWM and information content, barcode quality and its calibration, coverage, trimming, barcode
space, the `suggest` cycle trace, overrepresented k-mers, the cell rank curve, consensus quality,
error and layout, and thread scaling. It reads no reads and computes nothing, so a figure can be
redrawn from the table beside it long after the FASTQ is gone and can never disagree with the
number in the report. gnuplot is not a Python dependency: without it the `.gp` scripts are still
written. `migec plot` joins `info` and `sheet` outside the five-pipeline-command rule.

**`migec suggest` now reports overrepresented k-mers**, and stitches them back into the sequence
they came from. Run it on a stage's *output* to find what the trim left behind: an 8-mer occurs by
chance every ~65 kb, so a surviving primer appears as a run of k-mers each shifted one base from
the last. Counted exactly in a flat 4^8 array — no hash map, nothing to size — and measured
against the reads' own base composition, never a flat 1/4, so the table is not just a description
of GC content.

**`checkout` writes `checkout.trimming.tsv`**, the payload length distribution after trimming. A
pattern matched one base off still matches and still trims; it just leaves every payload one base
short, which no counter of matched reads can show. The mean payload length is now a column in the
report.

**Fixes**

- `checkout`'s report divided by zero on an empty input file, so a run with no reads ended in a
  traceback rather than a summary.
- The benchmark corpus assigned every read to one sample of four (`i % 4` with four reads per
  molecule is always 0). The matcher still scored all four patterns, so the throughput figures
  stand, but the per-sample counters and the memory figure were a single sample's. Fixed in
  `tests/benchmark/` and in the new `scripts/benchmark_threads.py`, and the published numbers are
  re-measured: 1,056,472 reads/s end to end and 1,684,654 matching at 16 threads, 217 MB.
- A sample that received no reads reported an effective barcode length of `inf`. Infinity is not a
  length; it is the absence of one, and it was being written into a numeric TSV column.
- The consensus log-likelihood tables were rebuilt per molecule — 122 transcendentals per group
  against the ~3 per emitted base the posterior actually needs.
- `polars` was a runtime dependency of the wheel and nothing in the package imported it; it is now
  in the `notebooks` extra with `marimo` and `matplotlib`. The pipeline's only runtime dependency
  is `typer`.
- Stale documentation corrected: `docs/roadmap.rst` still called `assemble` and `refine` planned;
  `docs/checkout.rst` said dual-end barcodes were unimplemented; `docs/fragmented.rst` documented
  `--mode {amplicon,fragmented}`, which shipped as `--contig`; the README quoted the retracted
  1.04x position-independence excess (it is 1.0103x) and linked the docs badge at a host the docs
  are not published on.
- `_not_yet()` and its claim in the skill that `subsample` exits 2 are gone; it has worked since
  M4.

## 2.0.0a2 — 2026-08-13

**Positional is the primary mode.** Most libraries fix the barcode at an offset in one read, and
saying so no longer takes a sample sheet or a flag:

```bash
migec checkout reads.fq.gz --bc-pattern '^NNNNNNNN' -o out/      # a caret anchors it
migec checkout reads.fq.gz --bc-pattern '0:8'       -o out/      # or a half-open slice
migec checkout R1.fq.gz R2.fq.gz --bc-pattern 'cell:0:16,16:26' -o out/
migec checkout R1.fq.gz R2.fq.gz --preset 10x-v2 -o out/
```

Slices are half-open and 0-based like Python's, each a UMI slice unless prefixed `cell:`, and gaps
between them become skipped bases — which is what a spacer is.

**`--max-offset` is now automatic and should not be passed.** A caret, a slice list, a read
structure and a pattern with nothing to score all anchor at the first base. This was the sharpest
edge in the interface: 10x and TSO500 needed `--max-offset 0` typed by hand, and without it every
read was refused with an error about anchoring — correct, but only after the run had failed.
Passing `--max-offset -1` explicitly still reinstates the refusal, and that is still the right
answer: a free scan over an unanchored pattern has no evidence to choose an offset with.

**Presets** — `umi`, `migec`, `primerid`, `duplex`, `10x`, `10x-v2`, `tso500`, `smarter-umi`.
`migec sheet --presets` prints each with what it is and where the layout is written down; a wrong
name lists all of them. Every preset is a published chemistry with a citable source, so a wrong one
is falsifiable rather than folklore, and each is compiled by a test — a preset nobody can run is
worse than no preset, because it looks supported. Two carry warnings that are part of the preset:

* Never: `duplex` extracts the tags and emits **single-strand** consensuses. Duplex pairing is not
  implemented and no duplex error rate may be quoted from it.
* Never: **TSO500's UMI is 5 nt, on R1 only** — the read structure is `5M5S+T +T`, checked against
  the pipeline this repo was pointed at, which corrects a first pass that had it on both mates.
  1,024 barcodes does not identify a molecule on a real ctDNA panel, and TSO500's own pipeline does
  not claim it does: it groups on the UMI **and the mapping position** (`fgbio GroupReadsByUmi`,
  after alignment). migec groups on the barcode, before any alignment exists, so it will report the
  space as saturated — on this chemistry that warning is the correct answer, not a threshold to
  raise. Extract and tag here; group position-aware, downstream.

**A run that matches nothing now says so, and says nothing else.** Zero assigned reads was rendered
indistinguishably from a successful run of an empty file: the table printed, and three warnings
computed from the reads that never arrived followed it — "base composition costs 100% of the
barcode space", a barcode error of `0.0e+00`, "under-sequenced at 0.0 reads/UMI" — none of which is
a fact about the library. It now reports the declaration error and suppresses the rest. Found by
the `smarter-umi` preset, whose first draft **scored** the template-switch `GGG`: three matching
bases are 6.0 bits against an anchored acceptance bar of 6.64, so it refused every read. The preset
skips those bases instead, as the source pipeline does, and the warning names that trap.

**The downstream contract is measured, not asserted** — new `docs/downstream.rst`. Against a real
`assemble` output on this machine: `minimap2 -ax sr -y` and `bwa mem -C` carry `RX`, `CB` and `MI`
into a valid sorted BAM on 600/600 records; `arda amplicon` reads it directly and its AIRR
`sequence_id` **is** the molecule id, because `dnaio` drops FASTQ comments and the read name was
built to be self-sufficient for exactly that reason; `salmon` and `kallisto` quantify it plainly.
Never: do not run alevin, bustools or STARsolo on a consensus FASTQ — they deduplicate from a raw
barcode read that no longer exists, and one consensus is already one molecule. STAR could not be
confirmed here: the Homebrew arm64 build reports zero input reads for any FASTQ, including a
one-record file, so the failure is the build's.

**`smarter-umi`, and a dataset that cannot be reprocessed.** SRP150352 (Sci Rep 2018,
doi:10.1038/s41598-018-31064-7) is the reference UMI RNA-seq library for the layout, but its
pipeline moves the UMI into the FASTQ *header* and SRA rewrites headers, so the archived copy has
no UMI anywhere — confirmed on three runs by the absence of the template-switch `GGG` where it
would have to be. `migec suggest` reports it unprompted rather than fitting a barcode to payload.
The preset is therefore sourced from `ncgr/UMI-analysis` itself, and `SOURCES.md` records both.

**Translating from zUMIs**, which is what NASC-seq2 drives. Never: zUMIs ranges are 1-based and
inclusive, migec slices are 0-based and half-open, so `UMI(12-19)` is `11:19` and not `12:19` —
subtract one from the start and leave the stop alone. `BC(1-6,20-26)` becomes
`cell:0:6,cell:19:26`. NASC-seq2 is an alternative to `checkout` rather than something downstream
of `assemble`; feeding it a consensus FASTQ would collapse twice.

Also: README badges; `docs/layouts.rst` collecting all four ways to declare a layout; the Nextflow
module takes `migec_preset` and no longer forces an offset.

## 2.0.0a1 — 2026-08-13

First published build of the rewrite. All three stages work and are validated against real data on
four layouts: a bulk amplicon with an anchor (HIV-1 Primer ID, `SRR1763769`), a 10x droplet VDJ
library (`sc5p_v2_hs_PBMC_1k`), dual-end barcodes (MAGERI's design), and TSO500 read structures.

An alpha because the roadmap's remaining items are real: index hopping from the i7 x i5 table,
`.mig` bucket output from checkout (which is what would bound the UMI counters at NovaSeq scale),
and the published benchmark comparisons. Install with `pip install --pre migec`.

Everything below was developed before this tag.

## 2.0.0a1 and earlier

### Cell barcodes, whitelists, cell calling, dual-end barcodes, and what the Phred is worth

**`X`/`x` is a cell-barcode position** in the pattern grammar — the one extension to MIGEC's
dialect, chosen because `X` is not a IUPAC symbol so no published table can contain one. `checkout`
writes `CB:Z:`/`CY:Z:`; `refine` keys on cell+UMI concatenated, because the same UMI in two cells is
two molecules and a UMI-only table would correct them against each other.

**`refine --cell-whitelist`** snaps a barcode to the list that was actually synthesised. The
load-bearing part is the competing hypothesis: without "this barcode is not on the list and was
read correctly", every hopped or undeclared barcode is absorbed into its nearest entry at posterior
1.0. Its prior is measured from barcodes at distance ≥2 from every entry — and it is a prior on
*this barcode*, the off-list read share divided by the distinct off-list barcodes.

**Cell calling** with OrdMag, the knee reported beside it and a warning when they disagree past 3×.
Exactly 500 of 20,500 on a synthetic droplet library. EmptyDrops-style rescue is Cell Ranger's job
and is deliberately not reproduced.

**Dual-end barcodes**: column 3 of the sheet is MIGEC's slave pattern, on the other mate, extending
the UMI. MAGERI's `NNNNNNNNNNNNtgact` / `agtcaNNNNNNNNNNNN` gives a 24 nt UMI. Both halves must
match or the read is unmatched.

**`--max-offset`**, and the bug it exposed: the acceptance bar is a Bonferroni bound over offsets,
and it was charging for the offsets a read *could* hold rather than the ones actually scanned. A
5 nt dual-end handle is 10 bits and was billed 12.6, so an anchored scan refused every read of a
design that is perfectly well determined. A *free* scan refusing it is correct — `TGACT` occurs by
chance every kilobase.

**What the reported Phred is worth**, measured against the pattern's own constant bases: slope 1.04
over 46.3 M bases on `SRR1763769`. Never: The fit's intercept (3.9e-3) is the **synthesised primer's**
defect rate, not a sequencing floor — it is even across all 23 anchor positions with none
polymorphic, and agrees with the independently measured 0.55% one-base-short rate. Reported, never
applied.

**The MIG-size FDR threshold**, measured rather than derived, and reported rather than applied.
Note: A count-ratio criterion reports zero residual at 1–3 reads/UMI, which is where it is worst;
payload agreement is what still works there. 5.25% of 1-read molecules at 1.23 reads/barcode
against 0% at 4.62.

### `migec refine`: the stage that decides how many molecules there were

Reads checkout's tagged FASTQ, folds error-child barcodes into their parents using all three
evidence terms, and rewrites the reads with the corrected barcode in `RX` and **the original
preserved in `OX`** — a correction nobody can audit is a correction nobody can check. Writes
`<sample>.barcodes.tsv` (umi, reads, corrected reads, parent) and the post-correction coverage
histogram.

On 20,000 simulated molecules at 6 reads/UMI with 3·10⁻³ barcode error: 23,910 distinct barcodes
in, **20,055 molecules out**, per-base error estimated at 2.87·10⁻³. The report prints the measured
clonality next to what it buys, because payload agreement is worth `log(1/clonality)` and that is a
property of the library rather than of the method.

It holds the **barcode table**, never the reads — `(key, count)` plus the per-position mean error
and a 32-base payload draft — and streams the reads three times instead. 2.2 MB of table for
110,349 reads.

**Note: Correction is not bucketable by a plain range partition.** The top *b* bits of the key decide
the bucket, so a substitution in the top *b*/2 positions sends a barcode and its neighbour to
different buckets and the pair can never be found. Two passes with the key rotated fixes it; until
then the table is held whole and its size is reported, exactly as `checkout` does with its counters.

Also fixed: `merged_reads` double-counted chained merges. When x folds into y and y later folds
into z, `corrected[y]` already carried x's reads by then. It is now each barcode's own count, since
a read changes barcode once.

### UMI correction learns to work at 1–3 reads per UMI

`correct_umis` ran only inside checkout and was unreachable from Python, so its accuracy had never
been scored against a known truth. It is now bound, and `scripts/correction_accuracy.py` scores it
across depth against the simulator's truth files.

The first score said the count ratio — the only evidence it used — is the whole game on a deeply
sequenced amplicon and carries **nothing** below ~3 reads/UMI: a parent with 2 reads and a child
with 1 is not an asymmetry, and two singletons are not one either. Two count gates made a
singleton-vs-singleton merge impossible by construction. Three changes:

- **The barcode's own base quality** at the position that differs. A sequencing miscall carries a
  low Phred exactly there; an early-PCR child carries a high one in every read. `QX` was already
  travelling through checkout unread. Works at one read.
- **Payload agreement** with the candidate parent, because a barcode error child is a read of the
  *parent's molecule*. It is worth `log(1/clonality)`, and the clonality is measured by sampling
  random barcode pairs — decisive in a diverse repertoire, worth nothing in a clonal library, and
  the reported number says which this one is. Agreement lifts the count gates; disagreement refuses
  a merge the count ratio would have made.
- **The error likelihood is a rate, not a conditional.** It was a zero-truncated Poisson weighed
  against an expected *count*. The truncation divides out `(1 − e^−λ)`, which is precisely the term
  saying whether an error child should exist — so `ZT-Poisson(1, λ) → 1` for every small λ, and the
  error rate stopped mattering at exactly the coverage where nothing else was available.

Scored against the achievable ceiling, because a child whose parent barcode was never sequenced has
nothing to merge into and correctly stays put:

| reads/UMI | reachable | recall of those | precision | molecules kept |
|---|---|---|---|---|
| 1.11 | 0.204 | 0.108 | 0.818 | **1.000** |
| 2.32 | 0.904 | 0.816 | 0.830 | 0.987 |
| 3.12 | 0.975 | 0.914 | 0.926 | 0.991 |
| 7.12 | 1.000 | 0.979 | 0.997 | 0.999 |
| 13.30 | 1.000 | 0.983 | 0.999 | 1.000 |

Note: At ~1 read/UMI **80% of barcode errors are unfixable in principle**. Of the rest migec fixes 11%
and destroys no real molecule at any depth — which is the side to err on, because a wrong merge
deletes a molecule and nothing downstream can tell, while a missed correction only inflates a count.

### `migec assemble`: one consensus per molecule

**A molecule is sample + cell barcode + UMI.** Never the UMI alone — the same UMI turning up in two
cells or two samples is the design, not an error, and 4^12 random tags reused across ten thousand
cells is exactly what a droplet protocol does. The sort key is `(cell, umi, src_index)` and the
range partition is on the cell whenever there is one, which also makes a per-cell scope contiguous
on disk.

**Nothing scales with the library.** Reads are range partitioned into `.mig` buckets and one bucket
is sorted in RAM at a time. 531,365 reads/s; 121 MB at 16 buckets against 203 MB at one, asserted
in `tests/benchmark/test_assemble_speed.py` — each configuration measured in its own process,
because `peak_rss_bytes` is a process high-water mark and two runs in one interpreter cannot be
compared. Note: The writer buffer budget is split *across* buckets rather than being per-bucket: a
fixed per-writer block made cutting the input finer cost *more* memory (238 MB at 16 buckets
against 213 at one), which is backwards.

**The quality floor is added, not compared.** `Q(j) = −10 log10(p_cons(j) + p_floor)`. An RT or
first-cycle-PCR error is in every read of the molecule and no consensus removes it, so the two
failure modes are independent and the emitted quality carries both. Default 1e-4 from X2, so
nothing above ~Q38 is emitted.

**Splitting a group uses X3's measured threshold**, 8.68, two-sided and Bonferroni'd within the
group. Note: It implies a minimum group size: the strongest evidence a pair of columns can carry is
`log10 C(n, n/2)`, so a 50/50 split needs about 34 reads before it can clear 8.68 at all. Below
that the data cannot separate a subclone from two bad reads at a 1% false-positive rate.

**`--contig` for random-primed libraries.** Reads sharing a barcode tile the molecule rather than
starting at the same base (X1: true for 92% of 10x groups). They are placed against each other by
exact seed matching, cut into overlap components by a union-find that carries each read's offset,
and one consensus is emitted per component — never bridged across a gap, because 27.3% of groups
hold more than one and a single consensus over those asserts sequence no read covers. This is one
molecule's fragments and nothing more: full-length receptor assembly, doublet calling and
contaminating-chain filtering are arda's job. Note: It also needs a barcode that is not saturated, so
`assemble` re-runs the birthday arithmetic on the barcodes it saw and reports
`expected_molecules_per_group`.

### X3: three derivations replaced by three permutations

Each of these numbers came out of an argument that assumed something the data had never been asked
about. `scripts/permutation_nulls.py` measures all three on `SRR1763769` — 125,369 distinct 9 nt
barcodes at 47.8% occupancy — and `docs/nulls.rst` has the tables.

**Position independence holds to ~1%.** The null is a *distribution* — the product measure
`q(u) = Π_j p_j(u_j)` — so it is tested by Jensen-Shannon divergence against a same-size draw from
`q` (a column shuffle), which measures the sparsity floor instead of assuming it. Invisible on all
distinct barcodes (at 47.5% occupancy the observed set is nearly a complete enumeration, which is
uniform by construction); clear once singletons are dropped, z up to 33; and **entirely
nearest-neighbour** — every adjacent position pair positive, every distant pair zero.

The cause is measurable: **0.55% of reads carry a barcode one base short**, a coupling step that
did not fire, and a frameshift is exactly a nearest-neighbour correlation, largest next to the
anchor. `Π_j Σ_a p_j(a)²` stays; the 1.86× collision excess is the read threshold.

Note: The first version of this null reported 1.04× and **all of it was artefact**. Two defects, both
now guarded: an `N` counted as a fifth base let `m_j` fall to 0.2466 — below the mathematical floor
of 1/4 — and printed an effective length of 9.01 nt for a 9 nt barcode; and the plug-in `Σ p̂²` is
biased up by `(1 − Σp²)/n`, a bias that *grows* as the distribution spreads and so reads as
dependence accumulating with k. `N` now folds to `A` as `pack_barcode` stores it, and the collision
is the U-statistic `Σ nₐ(nₐ−1)/(N(N−1))`.

**92% of distance-1 barcode pairs are coincidence.** 839,218 observed against 773,684 under a
column shuffle, which keeps every marginal and destroys the error children. Shuffling the read
*counts* over the fixed distance-1 graph — the graph, the composition and the count distribution
all unchanged, only which count sits on which node — finds **~19,400 genuine parent/child pairs**,
plateauing from a count ratio of 5 upward at z ≈ 48. The permuted background puts barcode error at
3.4e-3, within 1.7× of the Phred + polymerase prediction, where the analytic estimate is 2.6×
*below* it. M3's error model takes the permuted background.

**The split threshold is 8.68, not 2.00.** Reads are not exchangeable: a low-quality read carries a
minor base at many positions at once and is indistinguishable from a linked subclone if you only
look at the columns. Randomising the reads × positions minor-allele matrix while preserving *both*
margins (curveball) puts the 1% false-positive point at a Bonferroni'd `-log10 p` of 8.68. The
nominal `p < 0.01` the derivation gives calls 30.62% of MIGs against 1.60% — a 19× over-call, every
one of which would have become a spurious extra molecule.

Note: A tail quantile is not a constant until its Monte Carlo error is smaller than the digits
quoted.** This read 9.91, then 9.61, then 11.66 on reruns of ~8,000 randomisations. At 82,800 it is
**8.68, bootstrap 95% CI [8.42, 9.14]**, and the interval is what the docs quote.

### Shallow libraries are a first-class case

Bulk repertoire profiling and shallow 3' single-cell both put the MIG size histogram's mass at 1–3
reads per UMI. Nothing is thresholded away there — `--min-reads` defaults to 1, because a molecule
seen once is still a molecule and the answer to a barcode error is correction, not a cut — and the
report says that the UMI is buying counting rather than error correction. Three results calibrated
on a deep library are documented as *not* transferring: the split threshold is inert (a column pair
carries at most `log10 C(n, n/2)`, so it needs ~30 reads), the count-ratio error-child null has no
dynamic range at 1–3 reads, and singleton filtering costs 79% of barcodes rather than 56%.

It is also the memory-hostile shape, since distinct barcodes are what everything in `assemble`
scales with, so the benchmarks use it: 190,595 reads/s at 1.02 reads/UMI, 259 B resident per
distinct barcode, still bounded by the bucket rather than the library.

### `migec suggest`: read the barcode layout off the reads

A UMI cycle is one the synthesiser mixed — all four bases near 1/4, ~2 bits. A constant cycle is
one base near 100%. `suggest` segments the per-cycle composition on that and prints a paste-ready
pattern. On `SRR1763769`, with nothing supplied but the FASTQ, it recovered a 9 nt Primer ID
followed by `CAGTTTAACTTTTGGGCCAT`, and checking that pattern out assigned 95.0% of reads.

Note: The pattern stops at the last *constant* run. Composition alone cannot tell a UMI from diverse
payload — both are four flat lines at 25% — and what separates them is that a barcode is anchored
and payload is not. A uniform run with nothing constant after it is reported in the note and left
out, because claiming it would produce a pattern that matches everywhere.

### The barcode space and the error budget, reported on every run

Both were arithmetic anyone could do from the existing output and nobody did, and both change what
the numbers mean.

**Barcode space.** `4^L` over the *captured* positions (`NNNNtNNNNtNNNN` captures 12, not 14) is
nominal; a real oligo mix is not 25/25/25/25, so the usable space is `1 / Π_j Σ_a p_j(a)²` and
`bias_loss` is the shortfall. From there the birthday problem in the form that survives a full
space: molecules land independently, so occupancy is Poisson, `occupied = S(1 − e^−λ)` pins λ, and
**`p_multi` = P(k>1 | k≥1)** is the fraction of MIGs that are two or more molecules pooled. That is
the number that matters — their consensus is a mixture of templates and over-sequencing cannot fix
it. `checkout` warns past 5%.

**Error budget.** The distance-1 estimate is now printed next to what predicts it: `⟨10^(−Q/10)⟩`
over the barcode bases plus `ε_pol × cycles`. Note: The Phred term is the mean of the *probabilities*,
not `10^(−mean Q/10)` — the function is convex, so half at Q40 and half at Q10 is 5%, not the 0.3%
"mean Q25" suggests.

Written to `checkout.barcode_space.tsv` and `checkout.umi_quality.tsv`, warned on in the report,
derived in `docs/barcode_space.rst`, drawn in `notebooks/barcode_space.py`, tested in
`tests/synthetic/test_barcode_space.py`.

### Two errors the checks found

**Never: The UMI error estimator was 3× low, everywhere.** Its expected-children term used ε where it
had to use ε/3: a sequencing miscall has to land on one specific alternative base out of three, and
using ε makes the expectation 3× too large and the solved rate 3× too small. Against an injected
3·10⁻³ it returned 0.31× at every occupancy from 0.3% up. It now returns 0.92×.

**Never: ...and it fails downward as the barcode space fills.** The estimator subtracts the coincidence
expectation from the observed distance-1 pair count; once most of a barcode's 3L neighbours are
themselves real barcodes, that is a small difference of two large numbers. Measured against the
same injected rate: 0.92× at 0.3% occupancy, 0.65× at 16%, 0.23× at 50%, 0.001× at 93%. The
collapse is always downward, so a crowded library under-reports its own barcode error and
under-corrects. `err_unreliable` is set past 5% neighbourhood occupancy.

**The birthday prediction checks out to within a factor of two, in the direction it should.**
`scripts/collision_check.py` measures collisions model-free — two molecules sharing a barcode with
different sequences are visible in the reads — and finds 1.86× the prediction on the HIV library.
`Π_j m_j` assumes the positions are independent and is a *lower* bound on the collision
probability, so more collisions than predicted is what a real synthesiser should produce.

### Audit fixes

A read of the checkout path start to finish, each finding reproduced before it was fixed.

**Never: Two barcode rows declaring the same sample destroyed the output.** A MIGEC barcode table
writes a sample sequenced with more than one tag as several rows sharing the id — which
`sheet.py` documents — and checkout opened one file per *row*, so both rows `fopen`ed the same
path and interleaved two `FILE*` into it. The result was not a truncated FASTQ but a file that is
not a gzip stream at all, while the summary reported a clean run. Rows are now grouped by id: one
output file, one UMI counter, one summary row, and rows disagreeing about the UMI length are an
error rather than a counter holding two lengths at once.

**Never: An exception on a worker thread aborted the process.** It propagated out of the thread
function and hit `std::terminate` — SIGABRT with no message and no output flushed. Reachable
today: `BarcodePattern::compile` accepted a pattern capturing more than 32 UMI bases, which threw
from `pack_barcode` on a worker. Workers now capture and the driver rethrows on the caller's
thread, and the length is checked where the pattern is compiled, so the error names the row that
caused it.

**Note: The offset prune defeated the placement margin.** An offset was abandoned once it could no
longer reach the incumbent best — but a runner-up only has to reach `best − min_margin` to make
the placement ambiguous. Those offsets were dropped silently and the margin came back as the full
score, so a read with two placements 2.6 bits apart was reported as an unambiguous match at
whichever came first. The bar now leaves the margin's worth of room.

**Note: The UMI error-rate estimator assumed a uniform base composition.** The distance-1 shell is
`P_coll · Σ_j (1−m_j)/m_j`; the code used `3L · P_coll`, its `m_j = ¼` special case. On a skewed
UMI that overstates the independent term and so *underestimates* the error rate — the direction
that leaves errors uncorrected, and the same mistake as reaching for Shannon over Rényi.

Also: a sample present in few chunks no longer accumulates one empty gzip member per chunk (but a
sample with no reads at all still gets exactly one, because a zero-byte file is not a gzip stream
and `gzip -t` rejects it); the acceptance threshold is derived per pattern rather than from the
first row's length; `estimate_umi_error` no longer holds a reference into the entry array across a
call that may flush it; and `TrimMode::kPatternOnly`, which was a second name for `kPattern`, is
gone.

### The reported throughput was the matcher's, not checkout's

`wall_seconds` stopped at the demultiplexing driver, leaving the per-sample statistics — coverage
histogram, composition, count correction — outside the clock. They are serial and cost ~1.5–2 µs
per *distinct* UMI, so on 400 k reads at one read per molecule the reported figure was
2,192,882 reads/s against 438,456 actually elapsed.

`wall_seconds` now covers the whole call and `match_seconds` reports the driver separately, because
the two scale with different things and only one of them threads. Every published number has been
re-measured; the table in `docs/performance.rst` gains a matching column, and the benchmark asserts
that the stats stage stays inside the clock.

### Paired input, and strand normalisation

`migec checkout R1.fq.gz R2.fq.gz` looks for the tag in R1 first and, only if R1 came up empty, in
R2 — so the extra work falls on reads that would otherwise be discarded. When the tag turns up
there the pair is swapped, so the output R1 always carries it; single-end input gets the same
fallback against the reverse complement.

Note: This is not a convenience. Amplicon libraries are sequenced in both orientations, and a MIG
holding both orientations of one molecule loses half its reads at consensus while nothing upstream
reports it. The flipped count is in the report and in `normalised`.

The mate is passed through whole — trimming it needs a second tag, and dual-end barcodes are not
implemented yet — but both mates carry RX/QX/BC, because a tool that sees only one of them cannot
group the pair.

### Multi-core, at 1.18 M reads/s, with the output independent of `-t`

`--threads` defaults to one per core. **The output is byte-identical whatever it is set to**: reads
are matched in fixed-size chunks and the chunks are written back in input order. A demultiplexer
whose output depended on its thread count would produce results that could not be compared between
runs, and the failure would be invisible.

2 M single-end 115 nt reads over four patterns went from 53.7 s to 1.7 s. Three things were wrong,
each measured:

- **`log2` in the per-base scoring loop was 90% of runtime.** The score depends only on the
  reported Phred and the size of the IUPAC set, both small integers, so it tabulates into 1.2 kB.
- **zlib at level 6 compresses random DNA at 7 MB/s.** Read payload is close to incompressible, so
  compression on the serial writer capped throughput no matter how many threads were matching. Each
  worker now gzips its own chunk and the writer appends bytes; concatenated gzip members are a valid
  gzip stream (RFC 1952 §2.2), so the output is an ordinary `.fq.gz`.
- **The default compression level is now 1, not 6** — 137 MB/s for 13% more bytes.

An offset that cannot reach the acceptance threshold is now also abandoned mid-scan rather than
scored to the end.

### UMI counters that fit: a sorted array, not a hash map

`UmiCounts` is a bounded append buffer folded into a sorted `(key, count)` array. **~22 bytes per
distinct UMI, measured, against ~48 for `unordered_map<uint64_t, uint32_t>`** once nodes, the
cached hash and the bucket array are counted — at the 4·10⁸ distinct UMIs of an ordinary NovaSeq
run, 8.8 GB against 19 GB. Sorted order is also what the range partition and the 1-substitution
neighbourhood search both want. `CorrectionResult` is indexed in parallel with the entry array for
the same reason, and the buffer grows with the data rather than costing a fixed ceiling per sample.

Note: 8.8 GB still does not fit a laptop, and the counters are **not yet partitioned**. The fix is the
range partition, with `.mig` bucket output in M2. Until then checkout warns past 1 GB rather than
letting you find out from the OOM killer.

Fixed while rewriting: a child merged early could point at a parent merged later in the walk, so
`root` is now flattened in a final pass and the documented invariant actually holds.

### Speed and memory are reported on every run

Wall clock, reads/s, thread count, peak RSS and the UMI-counter share of it, in the report and in
`checkout.json`. `tests/benchmark/` (behind `RUN_BENCHMARK=1`) guards throughput, thread scaling,
bytes per distinct UMI, and output determinism across thread counts.

### Grouping accuracy against Calib

`scripts/compare_calib.py` scores read partitions against a known truth by adjusted Rand index,
reporting **splitting and merging separately** — splitting inflates the molecule count and is
recoverable, merging mixes molecules and destroys real variants. Calib clusters on barcode *and*
sequence; migec today groups on the barcode alone, and the measured gap is exactly the collision
rate: ARI 1.0000 on a clean 12 nt barcode, 0.8877 at 6 nt with 40% of reads merged. The migec
column is asserted on every test run; the Calib column needs Calib built locally.

The simulator gained an `adapter` field, without which its reads carry no constant region for a
barcode pattern to anchor on and cannot be checked out at all.

### Checkout: degenerate barcode patterns, trimming, header transfer

The barcode side of the pipeline, against the published MIGEC barcode tables, which are read
verbatim. Acceptance is a quality-aware log-likelihood ratio in bits rather than a mismatch count:
a match is worth +2.00 bits, a mismatch −9.55 bits at Q30 and −0.60 bits at Q2, so a mismatch on a
bad base is nearly free and one on a good base is fatal. Lowercase is a half weight rather than a
wildcard — v1 treated it as matching anything and discarded the evidence. Both of v1's real
defects are reproduced as tests: quality indexed from the read start rather than the match offset,
and a dangling `else` that meant low-quality mismatches were never counted.

`--trim pattern` leaves exactly the biological payload; the barcode travels on in SAM-style
RX/QX/BC tags, TAB-separated because `bwa mem -C` and `minimap2 -y` copy the comment verbatim into
the SAM record.

Note: A read equidistant from two sample tags is reported as **ambiguous**, not assigned to one of
them. That is a different diagnosis from "unmatched" — barcodes too close together versus a wrong
pattern — and one counter cannot say both.

### UMI statistics, and which entropy is allowed near a decision

Coverage histogram in MIGEC's power-of-two bins, per-position base composition with Shannon
entropy and information content for a logo, and the collision statistics.

Never: A logo draws Shannon entropy. The probability two molecules collide is the Rényi entropy of
order 2, `Π_j Σ_a p_j(a)²`. Since H₂ ≤ H₁, using Shannon overstates the usable barcode space and
understates collisions — the direction that silently merges distinct molecules. Both are reported;
only the collision form feeds the effective length, the correction, and the molecule count.

Count correction weighs three hypotheses per neighbour pair: a sequencing miscall (zero-truncated
Poisson at ε/3, since a miscall lands on one specific base), a polymerase error (a Luria–Delbrück
1/f² tail — an early-cycle error reaches a large share of the family and carries high quality in
every read, so a Poisson on the sequencing rate rejects it and it survives as a spurious second
molecule), and the barcode belonging to another real molecule whose size is drawn from the
library's own MIG size distribution. An isolated low-coverage UMI has no parent and keeps its
reads.

Note: The collision correction is declined above 90% occupancy rather than reported: the barcode space
is estimated from the observed barcodes, so at saturation the estimate collapses onto the observed
count and would report "no collisions" for the most collided library possible.

### Spike-in validation

`scripts/spikein_ratio.py` computes the published MIGEC metric — a real spike-in variant against
the worst *error* at the same substitution distance. Raw reads give V1/Err1 ≈ 1.4 and V2/Err2 ≈
0.3; UMI consensus is expected to reach 26–76 and 4.6–6.2. V2 is *less* abundant than the worst
2-substitution PCR error, so no abundance threshold can separate them — which is the entire
justification for molecular barcoding.

Note: The junction is anchored on its 3′ end only. Requiring the 5′ end too looks more robust and is
catastrophically wrong: V1 differs at position 4 and V2 at 7–8, so both variants count as zero and
the metric looks perfect.

### The rewrite

MIGEC and MAGERI are replaced by a single C++20 core with a pybind11 module and a typer CLI. The
old Groovy implementation is archived on `legacy-v1` (tag `v1-final`); MAGERI's alignment and
variant calling are out of scope, and the pipeline now ends at consensus FASTQ that arda,
minimap2 and bwa-meme consume directly.

### The `.mig` intermediate format

One format between all stages, frozen before the stages were written and pinned by a round-trip
test. Three choices that look wrong until measured, all documented in `docs/formats.rst`: raw
ASCII sequence rather than 2-bit packing (packing measured *worse* — 227 vs 197 B/pair — because
it destroys the cross-read redundancy a compressor finds in amplicon data), column-major block
layout (interleaving sequence and quality costs the compressor 10–20%), and a u64 `src_index`
(u32 overflows on a NovaSeq X run, and it is the tiebreak that makes output byte-identical across
thread counts).

Files are **range** partitions of the sort key, not hash partitions. A hash sends a barcode and
its one-mismatch neighbours to uncorrelated buckets, which makes barcode correction impossible to
apply locally and splits a molecule permanently — with each half looking like a well-formed MIG,
so nothing downstream detects it. Range partitioning costs the same and makes bucket order equal
key order, so the on-disk sort by sample/cell/UMI is a property of the layout rather than a pass.

### FASTQ IO

Plain or gzipped, decided by content rather than file name. The reader materialises a whole
four-line record before handing out any view, because the buffer compacts on refill and a record
straddling that boundary would otherwise dangle the header view — a bug that shows up on roughly
one record in ten million and looks like data corruption.

### Read simulator

`tests/synthetic/_sim.py` generates reads with full ground truth: templates, UMIs from a settable
base composition, RT errors (present in every read, therefore un-removable by consensus), early
PCR errors that reach consensus only when they exceed half the descendants, and lognormal MIG
sizes so that the 1–5 read regime the retention rule is about is actually populated.
