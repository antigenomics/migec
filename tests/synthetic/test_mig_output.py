"""`checkout --mig` writes the partition `assemble` would have built, and assemble reads it.

The claim is not that it is faster -- that is `tests/benchmark/` -- but that it is the SAME
pipeline. A partition written by the wrong stage, or read back with the wrong bucket bits, groups
different reads together and produces a different set of molecules, which is the one failure a
consensus pipeline cannot survive. So what is asserted here is the molecule table, record for
record, against the FASTQ route.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


@pytest.fixture(scope="module")
def both_routes(tmp_path_factory):
    """One library, checked out twice -- to FASTQ and to buckets -- and assembled from each."""
    from migec.assemble import run as assemble
    from migec.checkout import run as checkout

    d = tmp_path_factory.mktemp("migout")
    cfg = SimConfig(n_molecules=4000, n_clones=20, umi_len=12, coverage=6.0, umi_error=1e-3,
                    adapter=ADAPTER, seed=23)
    sim = simulate(cfg, d / "sim")
    (d / "bc.txt").write_text(f"S1\t{'N' * cfg.umi_len}{ADAPTER.lower()}\n")

    fastq = checkout(sim["reads"], d / "bc.txt", d / "co_fastq")
    mig = checkout(sim["reads"], d / "bc.txt", d / "co_mig", mig=True)
    a_fastq = assemble(d / "co_fastq" / "S1.fq.gz", d / "as_fastq")
    a_mig = assemble(d / "co_mig" / "S1.000.mig", d / "as_mig")
    return d, fastq, mig, a_fastq, a_mig


def test_the_buckets_are_a_partition_of_the_sample(both_routes):
    d, fastq, mig, _, _ = both_routes
    assert not fastq["mig"]
    assert mig["mig"]
    assert mig["mig_paths"], "no bucket was written"
    # A bucket that received no read is absent rather than empty: a bucket is addressed by the
    # header inside it, so an absent one is an absent range of the key space.
    assert all(Path(p).exists() for p in mig["mig_paths"])
    assert not list((d / "co_mig").glob("S1.fq.gz")), "FASTQ was written as well"
    # Same reads, same barcodes: the route does not change what was matched.
    assert mig["assigned"] == fastq["assigned"]
    assert mig["samples"][0]["umis"] == fastq["samples"][0]["umis"]


def test_the_molecules_are_the_same_molecules(both_routes):
    _, _, _, a_fastq, a_mig = both_routes
    assert a_mig["reads"] == a_fastq["reads"]
    assert a_mig["groups"] == a_fastq["groups"]
    assert a_mig["molecules"] == a_fastq["molecules"]
    assert a_mig["groups_split"] == a_fastq["groups_split"]


def test_the_consensus_records_are_identical(both_routes):
    """Byte for byte after decompression, which is the assertion that leaves no room.

    The gzip framing does differ: the bucket count is not the same, and bucket boundaries are gzip
    member boundaries. What must not differ is a single base, quality or tag.
    """
    d, _, _, _, _ = both_routes
    a = gzip.open(d / "as_fastq" / "S1.consensus.fq.gz", "rt").read()
    b = gzip.open(d / "as_mig" / "S1.consensus.fq.gz", "rt").read()
    assert a == b
    assert (d / "as_fastq" / "S1.mig.tsv").read_text() == (d / "as_mig" / "S1.mig.tsv").read_text()


def test_assemble_does_not_eat_its_input(both_routes):
    """Never: only a bucket assemble WROTE may be removed.

    Pass 2 deletes each bucket as it finishes with it, which is right for the temporaries it made
    itself and catastrophic for checkout's output: the second run over the same sample would find
    an empty partition and report a smaller library, with nothing in the summary to say why.
    """
    d, _, mig, _, _ = both_routes
    from migec.assemble import run as assemble

    assert all(Path(p).exists() for p in mig["mig_paths"])
    again = assemble(d / "co_mig" / "S1.000.mig", d / "as_mig_again")
    assert again["molecules"] > 0
    assert all(Path(p).exists() for p in mig["mig_paths"])


def test_a_directory_of_two_samples_is_refused_by_name(both_routes, tmp_path):
    """A checkout output directory holds every sample's buckets, and assemble takes one sample.

    Never: assembling them together groups two samples' reads as one molecule wherever a UMI
    repeats -- which it does by design -- and nothing downstream can tell.
    """
    import shutil

    from migec.assemble import _mig_buckets

    d, _, mig, _, _ = both_routes
    for path in mig["mig_paths"][:2]:
        shutil.copy(path, tmp_path / Path(path).name)
        shutil.copy(path, tmp_path / Path(path).name.replace("S1.", "S2."))
    assert len(_mig_buckets(tmp_path / Path(mig["mig_paths"][0]).name)) == 2
    with pytest.raises(ValueError, match="per-sample stage"):
        _mig_buckets(tmp_path)


def test_a_limit_is_refused_on_a_partition(both_routes):
    """A limit is a prefix of the input, and a partition has no prefix left."""
    d, _, _, _, _ = both_routes
    from migec.assemble import run as assemble

    with pytest.raises(Exception, match="Limit at checkout"):
        assemble(d / "co_mig" / "S1.000.mig", d / "as_limited", limit_reads=100)
