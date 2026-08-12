# Changelog

Hand-written and prose-heavy: each entry says what changed and, where it matters, which failure it
prevents. Releases before 2.0.0 are the Groovy MIGEC and are described by their git tags on the
`legacy-v1` branch.

## Unreleased — 2.0.0.dev0

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
