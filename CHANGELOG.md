# Changelog

Hand-written and prose-heavy: each entry says what changed and, where it matters, which failure it
prevents. Releases before 2.0.0 are the Groovy MIGEC and are described by their git tags on the
`legacy-v1` branch.

## Unreleased — 2.0.0.dev0

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
