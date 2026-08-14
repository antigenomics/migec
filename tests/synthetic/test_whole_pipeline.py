"""The four stages as one chain, over a library whose truth is known.

Every stage is tested on its own, and `test_thread_invariance.py` asserts that each stage's own
output file does not move with `-t`. What neither covers is the chain: the invariants here are
about what one stage hands the next, and each of them can be violated while every stage still
reports a plausible summary of its own.

  * a read that a stage drops is a read no later stage can miss, because nothing downstream knows
    how many there should have been;
  * correction may only ever fold barcodes together, so a molecule count that RISES over `refine`
    is a merge that split something instead;
  * `cD` is the only place the depth of a molecule survives to, so if it stops summing back to the
    reads that entered, abundance is wrong and the consensus still looks right;
  * `subsample` keeps whole barcodes, which is worth nothing unless the rest of the chain then
    answers the same way about them;
  * and the bytes of EVERY file the chain writes -- not only the three the per-stage test digests
    -- must be a function of the input alone.

The corpus is ~9,000 reads over 2,000 molecules: the assertions are exact counts and file
comparisons rather than rates, so a larger library adds seconds and no evidence. What it does need
is a barcode error high enough that `refine` merges (otherwise the count invariant holds by having
nothing to correct) and a collision or two, so the truth is the number of distinct barcodes rather
than the number of molecules drawn.
"""

from __future__ import annotations

import collections
import gzip
import json

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"
THREADS = (1, 8)
KEEP_PERCENT = 25.0

# What `-t` is allowed to change, and the whole of it: two clocks, the memory high-water mark, and
# the thread count itself. Anything else that differs between the two runs is a result that moved.
CLOCKS = frozenset(
    {
        "threads", "wall_seconds", "match_seconds", "reads_per_second", "peak_rss_bytes",
        "table_seconds", "correct_seconds", "rewrite_seconds", "partition_seconds",
    }
)


@pytest.fixture(scope="module")
def chains(tmp_path_factory):
    """The same four-stage chain run at one thread and at eight, into two separate trees."""
    from migec.assemble import run as assemble
    from migec.checkout import run as checkout
    from migec.refine import run as refine
    from migec.subsample import run as subsample

    d = tmp_path_factory.mktemp("chain")
    cfg = SimConfig(adapter=ADAPTER, n_molecules=2000, n_clones=10, coverage=5.0, coverage_cv=0.8,
                    umi_len=12, umi_error=1.5e-3, seq_error=2e-3, seed=44)
    sim = simulate(cfg, d / "sim")
    (d / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")

    runs = {}
    for t in THREADS:
        root = d / f"t{t}"
        co = checkout(sim["reads"], d / "bc.txt", root / "co", threads=t)
        ss = subsample(root / "co" / "S1.fq.gz", root / "sub.fq.gz", keep_percent=KEEP_PERCENT)
        rf = refine(root / "co" / "S1.fq.gz", root / "rf", threads=t)
        asm = assemble(root / "rf" / "S1.fq.gz", root / "as", threads=t)
        # The same barcodes without correction, and the kept quarter of them, so the subsample can
        # be compared against a run that made the same merges -- which is none.
        whole = assemble(root / "co" / "S1.fq.gz", root / "as_all", threads=t)
        part = assemble(root / "sub.fq.gz", root / "as_sub", sample_id="S1", threads=t)
        runs[t] = (root, co, ss, rf, asm, whole, part)
    return d, sim, runs


def records(path):
    """FASTQ record count, read from the file rather than taken from a summary."""
    with gzip.open(path, "rt") as fh:
        return sum(1 for i, _ in enumerate(fh) if i % 4 == 0)


def reads_per_umi(path):
    counts = collections.Counter()
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                counts[next(f[5:] for f in line.split() if f.startswith("RX:Z:"))] += 1
    return counts


def consensus_records(path):
    """UMI -> (whole four-line record, cD depth)."""
    out, block = {}, []
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            block.append(line)
            if i % 4 == 3:
                umi = next(f[5:] for f in block[0].split() if f.startswith("RX:Z:"))
                depth = int(next(f[5:] for f in block[0].split() if f.startswith("cD:i:")))
                out[umi] = ("".join(block), depth)
                block = []
    return out


def test_no_read_is_created_or_destroyed(chains):
    """Counted from the files at every hand-off, because a summary counting its own output cannot
    catch a stage that wrote fewer records than it thinks it did."""
    _, sim, runs = chains
    root, co, ss, rf, asm, _, _ = runs[1]

    assert co["total"] == sim["n_reads"]
    assert co["assigned"] == co["total"], "the simulator writes only tagged reads"
    assert records(root / "co" / "S1.fq.gz") == co["assigned"]
    assert rf["reads"] == records(root / "rf" / "S1.fq.gz") == co["assigned"]
    assert asm["reads"] == rf["reads"]
    # Every read reached a group: one that did not is a read this pipeline has no other trace of.
    assert asm["reads_without_umi"] == 0
    assert asm["reads_dropped"] == 0

    # subsample is the one stage that is meant to drop reads, and it drops whole barcodes only.
    assert ss["reads"] == co["assigned"]
    assert ss["reads_kept"] == records(root / "sub.fq.gz") < ss["reads"]


def test_the_molecule_count_never_rises_after_correction(chains):
    """Correction folds barcodes together and can do nothing else.

    A count that went up would mean a merge had split a barcode's reads in two, which is the
    failure with no downstream symptom: both halves are well-formed MIGs.
    """
    _, sim, runs = chains
    root, co, _, rf, asm, _, _ = runs[1]

    observed = len(reads_per_umi(root / "co" / "S1.fq.gz"))
    corrected = len(reads_per_umi(root / "rf" / "S1.fq.gz"))
    assert observed == co["samples"][0]["umis"]
    assert corrected == rf["molecules"] <= observed
    assert rf["merged"] == observed - corrected > 0, "nothing was corrected; the test is vacuous"
    assert asm["groups"] == rf["molecules"]

    # ...and it moved towards the truth rather than merely downwards. The truth is the number of
    # DISTINCT barcodes drawn, not the number of molecules: two molecules on one barcode are one
    # molecule to every method there is.
    truth = sim["n_distinct_umis"]
    assert abs(corrected - truth) < abs(observed - truth)


def test_every_consensus_depth_sums_back_to_the_reads_that_entered(chains):
    """`cD` is where a molecule's abundance survives to, and it is derived, not carried.

    Checked per barcode against `refine`'s own output, not only in total: a read counted towards
    the wrong group leaves the sum right and the abundance wrong.
    """
    _, _, runs = chains
    root, _, _, _, asm, _, _ = runs[1]

    entered = reads_per_umi(root / "rf" / "S1.fq.gz")
    emitted = {umi: depth for umi, (_, depth) in consensus_records(
        root / "as" / "S1.consensus.fq.gz"
    ).items()}
    assert emitted == dict(entered)
    assert sum(emitted.values()) == asm["reads"]


def test_a_subsampled_library_gets_the_same_answer_for_the_barcodes_it_kept(chains):
    """What `subsample` is for, asserted where it matters rather than on its own output.

    A barcode is kept whole or not at all, so the reads under a kept barcode are exactly the reads
    that were under it in the full library -- and the consensus over them, its depth and its
    quality must therefore be the same record. Sampling READS would leave every one of these
    records intact in shape and wrong in depth, which is invisible in the subsample's own summary.
    """
    _, _, runs = chains
    root, _, ss, _, _, _, _ = runs[1]

    whole = consensus_records(root / "as_all" / "S1.consensus.fq.gz")
    part = consensus_records(root / "as_sub" / "S1.consensus.fq.gz")
    assert len(part) == ss["barcodes"] < len(whole)
    assert set(part) <= set(whole)
    assert {u: part[u] for u in part} == {u: whole[u] for u in part}


def test_the_whole_chain_is_byte_identical_at_one_and_eight_threads(chains):
    """Every file the chain writes, not the three a per-stage test digests.

    Never: a table that moves with `-t` is worse than a FASTQ that does. It is what a report is
    drawn from, nobody diffs it, and two runs of the same data then publish different numbers.
    """
    _, _, runs = chains
    (a, *_), (b, *_) = runs[1], runs[8]

    listing = {root: sorted(p.relative_to(root) for p in root.rglob("*") if p.is_file())
               for root in (a, b)}
    assert listing[a] == listing[b]
    assert len(listing[a]) > 20, "the chain wrote almost nothing; the comparison proves nothing"

    for rel in listing[a]:
        if rel.suffix == ".json":
            # The summaries carry the clocks by design. Everything else in them is a result, and
            # the path fields name the run's own directory -- so normalise that and compare the
            # rest key for key, which also catches a newly added field that is not deterministic.
            left = json.loads((a / rel).read_text().replace(str(a), "<root>"))
            right = json.loads((b / rel).read_text().replace(str(b), "<root>"))
            assert set(left) == set(right), rel
            for key in left:
                if key not in CLOCKS:
                    assert left[key] == right[key], f"{rel}:{key}"
        else:
            assert (a / rel).read_bytes() == (b / rel).read_bytes(), rel
