# Changelog

Hand-written and prose-heavy: each entry says what changed and, where it matters, which failure it
prevents. Releases before 2.0.0 are the Groovy MIGEC and are described by their git tags on the
`legacy-v1` branch.

## Unreleased — 2.0.0.dev0

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

⚠ **Correction is not bucketable by a plain range partition.** The top *b* bits of the key decide
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

⚠ At ~1 read/UMI **80% of barcode errors are unfixable in principle**. Of the rest migec fixes 11%
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
compared. ⚠ The writer buffer budget is split *across* buckets rather than being per-bucket: a
fixed per-writer block made cutting the input finer cost *more* memory (238 MB at 16 buckets
against 213 at one), which is backwards.

**The quality floor is added, not compared.** `Q(j) = −10 log10(p_cons(j) + p_floor)`. An RT or
first-cycle-PCR error is in every read of the molecule and no consensus removes it, so the two
failure modes are independent and the emitted quality carries both. Default 1e-4 from X2, so
nothing above ~Q38 is emitted.

**Splitting a group uses X3's measured threshold**, 8.68, two-sided and Bonferroni'd within the
group. ⚠ It implies a minimum group size: the strongest evidence a pair of columns can carry is
`log10 C(n, n/2)`, so a 50/50 split needs about 34 reads before it can clear 8.68 at all. Below
that the data cannot separate a subclone from two bad reads at a 1% false-positive rate.

**`--contig` for random-primed libraries.** Reads sharing a barcode tile the molecule rather than
starting at the same base (X1: true for 92% of 10x groups). They are placed against each other by
exact seed matching, cut into overlap components by a union-find that carries each read's offset,
and one consensus is emitted per component — never bridged across a gap, because 27.3% of groups
hold more than one and a single consensus over those asserts sequence no read covers. This is one
molecule's fragments and nothing more: full-length receptor assembly, doublet calling and
contaminating-chain filtering are arda's job. ⚠ It also needs a barcode that is not saturated, so
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

⚠ The first version of this null reported 1.04× and **all of it was artefact**. Two defects, both
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

⚠ **A tail quantile is not a constant until its Monte Carlo error is smaller than the digits
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

⚠ The pattern stops at the last *constant* run. Composition alone cannot tell a UMI from diverse
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
over the barcode bases plus `ε_pol × cycles`. ⚠ The Phred term is the mean of the *probabilities*,
not `10^(−mean Q/10)` — the function is convex, so half at Q40 and half at Q10 is 5%, not the 0.3%
"mean Q25" suggests.

Written to `checkout.barcode_space.tsv` and `checkout.umi_quality.tsv`, warned on in the report,
derived in `docs/barcode_space.rst`, drawn in `notebooks/barcode_space.py`, tested in
`tests/synthetic/test_barcode_space.py`.

### Two errors the checks found

⛔ **The UMI error estimator was 3× low, everywhere.** Its expected-children term used ε where it
had to use ε/3: a sequencing miscall has to land on one specific alternative base out of three, and
using ε makes the expectation 3× too large and the solved rate 3× too small. Against an injected
3·10⁻³ it returned 0.31× at every occupancy from 0.3% up. It now returns 0.92×.

⛔ **...and it fails downward as the barcode space fills.** The estimator subtracts the coincidence
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

⛔ **Two barcode rows declaring the same sample destroyed the output.** A MIGEC barcode table
writes a sample sequenced with more than one tag as several rows sharing the id — which
`sheet.py` documents — and checkout opened one file per *row*, so both rows `fopen`ed the same
path and interleaved two `FILE*` into it. The result was not a truncated FASTQ but a file that is
not a gzip stream at all, while the summary reported a clean run. Rows are now grouped by id: one
output file, one UMI counter, one summary row, and rows disagreeing about the UMI length are an
error rather than a counter holding two lengths at once.

⛔ **An exception on a worker thread aborted the process.** It propagated out of the thread
function and hit `std::terminate` — SIGABRT with no message and no output flushed. Reachable
today: `BarcodePattern::compile` accepted a pattern capturing more than 32 UMI bases, which threw
from `pack_barcode` on a worker. Workers now capture and the driver rethrows on the caller's
thread, and the length is checked where the pattern is compiled, so the error names the row that
caused it.

⚠ **The offset prune defeated the placement margin.** An offset was abandoned once it could no
longer reach the incumbent best — but a runner-up only has to reach `best − min_margin` to make
the placement ambiguous. Those offsets were dropped silently and the margin came back as the full
score, so a read with two placements 2.6 bits apart was reported as an unambiguous match at
whichever came first. The bar now leaves the margin's worth of room.

⚠ **The UMI error-rate estimator assumed a uniform base composition.** The distance-1 shell is
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

⚠ This is not a convenience. Amplicon libraries are sequenced in both orientations, and a MIG
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

⚠ 8.8 GB still does not fit a laptop, and the counters are **not yet partitioned**. The fix is the
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

⚠ A read equidistant from two sample tags is reported as **ambiguous**, not assigned to one of
them. That is a different diagnosis from "unmatched" — barcodes too close together versus a wrong
pattern — and one counter cannot say both.

### UMI statistics, and which entropy is allowed near a decision

Coverage histogram in MIGEC's power-of-two bins, per-position base composition with Shannon
entropy and information content for a logo, and the collision statistics.

⛔ A logo draws Shannon entropy. The probability two molecules collide is the Rényi entropy of
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

⚠ The collision correction is declined above 90% occupancy rather than reported: the barcode space
is estimated from the observed barcodes, so at saturation the estimate collapses onto the observed
count and would report "no collisions" for the most collided library possible.

### Spike-in validation

`scripts/spikein_ratio.py` computes the published MIGEC metric — a real spike-in variant against
the worst *error* at the same substitution distance. Raw reads give V1/Err1 ≈ 1.4 and V2/Err2 ≈
0.3; UMI consensus is expected to reach 26–76 and 4.6–6.2. V2 is *less* abundant than the worst
2-substitution PCR error, so no abundance threshold can separate them — which is the entire
justification for molecular barcoding.

⚠ The junction is anchored on its 3′ end only. Requiring the 5′ end too looks more robust and is
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
