"""Paired input all the way down: checkout, refine, assemble, and the same chain over `.mig`.

Paired reads are tested inside checkout and nowhere after it -- `test_checkout_end_to_end.py`
asserts that a pair whose tag sits on R2 is swapped, and the story stops there. What that leaves
untested is the thing the swap exists for: a MIG holding both orientations of one molecule. Half
its reads are on the mate, so a normalisation that silently dropped them leaves a molecule that is
still a molecule, still consensed, still in the table -- at half its depth, with nothing in any
summary saying so. That is only visible once the reads reach `assemble`, which is where this file
looks.

`_sim.simulate` writes single-end only (its `paired` flag is unused), so the pair is built here the
way `test_dual_end.py` builds one.
"""

from __future__ import annotations

import gzip
import random

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import BASES

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"
MOLECULES = 800
DEPTH = 5
BARCODE_ERROR = 0.02  # per read, not per base: an outright miscall in the barcode


def paired_corpus(d, seed=17):
    """`MOLECULES` molecules, `DEPTH` reads each, every other read of each molecule flipped.

    The flip is WITHIN a molecule, not between molecules: that is what puts both orientations under
    one barcode, which is the case checkout's fallback exists for and the case a per-read test
    cannot see. Barcodes are drawn distinct, so the number of molecules is known exactly and the
    assertions below are equalities rather than tolerances; `BARCODE_ERROR` then gives `refine`
    real work, without which the molecule count would come out right by having nothing to correct.
    """
    rng = random.Random(seed)
    truth, seen, flipped, i = {}, set(), 0, 0
    with gzip.open(d / "R1.fq.gz", "wt") as f1, gzip.open(d / "R2.fq.gz", "wt") as f2:
        for _ in range(MOLECULES):
            umi = "".join(rng.choice(BASES) for _ in range(12))
            while umi in seen:
                umi = "".join(rng.choice(BASES) for _ in range(12))
            seen.add(umi)
            payload = "".join(rng.choice(BASES) for _ in range(60))
            truth[umi] = payload
            for k in range(DEPTH):
                observed = list(umi)
                if rng.random() < BARCODE_ERROR:
                    j = rng.randrange(len(observed))
                    observed[j] = rng.choice([b for b in BASES if b != observed[j]])
                tagged = "".join(observed) + ADAPTER + payload
                mate = "".join(rng.choice(BASES) for _ in range(60))
                a, b = (mate, tagged) if k % 2 == 0 else (tagged, mate)
                flipped += k % 2 == 0
                f1.write(f"@p{i}\n{a}\n+\n{'I' * len(a)}\n")
                f2.write(f"@p{i}\n{b}\n+\n{'I' * len(b)}\n")
                i += 1
    (d / "bc.txt").write_text(f"S1\t{'N' * 12}{ADAPTER.lower()}\n")
    return truth, flipped


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    """One paired library down both routes: R1 FASTQ, and the buckets `checkout --mig` writes.

    4,000 pairs is enough because every assertion here is exact -- every molecule, every read --
    rather than a rate with a confidence interval. It buys ~50 error children at `BARCODE_ERROR`,
    which is what makes the correction non-vacuous, and the whole fixture runs in under a second.
    """
    from migec.assemble import run as assemble
    from migec.checkout import run as checkout
    from migec.refine import run as refine

    d = tmp_path_factory.mktemp("paired")
    truth, flipped = paired_corpus(d)

    fastq = checkout(d / "R1.fq.gz", d / "bc.txt", d / "co", reads2=d / "R2.fq.gz")
    mig = checkout(d / "R1.fq.gz", d / "bc.txt", d / "co_mig", reads2=d / "R2.fq.gz", mig=True)
    out = {
        "fastq": (
            refine(d / "co" / "S1_R1.fq.gz", d / "rf"),
            assemble(d / "rf" / "S1.fq.gz", d / "as"),
        ),
    }
    mig_refined = refine(mig["mig_paths"][0], d / "rf_mig")
    out["mig"] = (mig_refined, assemble(mig_refined["mig_paths"][0], d / "as_mig"))
    return d, truth, flipped, fastq, mig, out


def consensus(path):
    """UMI -> (sequence, depth) from a consensus FASTQ."""
    records = {}
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                umi = next(f[5:] for f in line.split() if f.startswith("RX:Z:"))
                depth = int(next(f[5:] for f in line.split() if f.startswith("cD:i:")))
            elif i % 4 == 1:
                records[umi] = (line.rstrip("\n"), depth)
    return records


def headers(path):
    with gzip.open(path, "rt") as fh:
        return [line.rstrip("\n") for i, line in enumerate(fh) if i % 4 == 0]


def test_the_mates_stay_in_step(pipeline):
    """Two files written by two code paths, read back by a downstream tool that assumes rank order.

    Never: an aligner given R1 and R2 pairs them by POSITION. One mate filtered without the other
    shifts every pair after it, and the run does not fail -- it produces a BAM of mismatched pairs.
    """
    d, _, _, fastq, _, _ = pipeline
    r1, r2 = headers(d / "co" / "S1_R1.fq.gz"), headers(d / "co" / "S1_R2.fq.gz")
    assert len(r1) == len(r2) == fastq["assigned"] == MOLECULES * DEPTH
    assert [h.split()[0] for h in r1] == [h.split()[0] for h in r2]
    # Both mates carry the barcode, or nothing downstream can group the pair.
    assert all("RX:Z:" in h for h in r2)


def test_orientation_normalisation_loses_no_read_by_the_time_it_reaches_a_consensus(pipeline):
    """The end of the chain is where a lost orientation shows as a number rather than as an absence.

    Every molecule here has three reads tagged on R2 and two on R1. Without the fallback the R2
    ones never match, and what comes out is the same 800 molecules at depth 2 -- a library that
    looks well-formed and is missing 60% of its evidence.
    """
    d, truth, flipped, fastq, _, out = pipeline
    assert fastq["normalised"] == flipped > 0
    assert fastq["unmatched"] == 0
    _, asm = out["fastq"]
    assert asm["reads"] == MOLECULES * DEPTH
    assert asm["reads_without_umi"] == 0 and asm["reads_dropped"] == 0

    records = consensus(d / "as" / "S1.consensus.fq.gz")
    assert {umi: depth for umi, (_, depth) in records.items()} == dict.fromkeys(truth, DEPTH)
    # ...and the reads that were flipped carried the payload, not the mate: a consensus built from
    # a mate that was never reverse-complemented is a random sequence at full depth.
    assert {umi: seq for umi, (seq, _) in records.items()} == truth


def test_the_molecule_count_is_the_number_simulated(pipeline):
    _, truth, _, fastq, _, out = pipeline
    rf, asm = out["fastq"]
    # Non-vacuous: the barcodes seen are more than the molecules that exist, and correction is what
    # closes the gap. Without merges this equality would hold by there being nothing to get wrong.
    assert fastq["samples"][0]["umis"] > MOLECULES
    assert rf["merged"] == fastq["samples"][0]["umis"] - MOLECULES > 0
    assert rf["molecules"] == MOLECULES == len(truth)
    assert asm["groups"] == asm["molecules"] == MOLECULES


def test_a_bucket_holds_both_mates(pipeline):
    """One bucket file, not two: the `.mig` record carries the pair together, so there is no second
    partition to keep in step. Asserting the mechanism, because the consequence -- identical
    molecules below -- would also pass if the mate were silently dropped at checkout."""
    from migec import _core

    _, _, _, _, mig, _ = pipeline
    bucket = _core.MigFile(mig["mig_paths"][0])
    assert bucket.header["paired"]
    records = bucket.read_all()
    assert records, "the bucket is empty"
    assert all(len(r.seq2) == 60 and len(r.qual2) == len(r.seq2) for r in records)


def test_the_bucket_route_gives_the_same_molecules(pipeline):
    d, _, _, _, _, out = pipeline
    (rf_fq, as_fq), (rf_mig, as_mig) = out["fastq"], out["mig"]
    for key in ("reads", "barcodes", "merged", "merged_reads", "molecules"):
        assert rf_mig[key] == rf_fq[key], key
    for key in ("reads", "groups", "molecules", "groups_split"):
        assert as_mig[key] == as_fq[key], key
    # Byte for byte after decompression. The gzip framing may differ -- the bucket count does --
    # but not a base, a quality or a tag.
    assert gzip.open(d / "as_mig" / "S1.consensus.fq.gz", "rt").read() == gzip.open(
        d / "as" / "S1.consensus.fq.gz", "rt"
    ).read()
