# Changelog

Hand-written and prose-heavy: each entry says what changed and, where it matters, which failure it
prevents. Releases before 2.0.0 are the Groovy MIGEC and are described by their git tags on the
`legacy-v1` branch.

## Unreleased — 2.0.0.dev0

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
