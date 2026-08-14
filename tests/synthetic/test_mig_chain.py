"""The whole pipeline over `.mig` buckets agrees with the whole pipeline over FASTQ.

`checkout --mig` writes the partition, `refine` corrects on it and writes it back re-partitioned on
the corrected barcode, `assemble` consenses from it. Every stage in that chain has a FASTQ twin,
and the only thing that may differ between the two routes is the container: not a merge, not a
molecule, not a base.

This is the test that would catch the failure the chain is built to avoid -- a partition that stops
being a partition once the key changes, which loses exactly the reads whose barcode was corrected
across a bucket boundary, and loses them silently.
"""

from __future__ import annotations

import gzip

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


@pytest.fixture(scope="module")
def chain(tmp_path_factory):
    """One library with real barcode errors, run down both routes end to end."""
    from migec.assemble import run as assemble
    from migec.checkout import run as checkout
    from migec.refine import run as refine

    d = tmp_path_factory.mktemp("migchain")
    # A barcode error rate high enough that correction has real work: at 1e-3 per base over 12
    # bases, ~1.2% of barcodes carry one, and a quarter of those land in the partitioned prefix.
    cfg = SimConfig(n_molecules=6000, n_clones=25, umi_len=12, coverage=6.0, umi_error=2e-3,
                    adapter=ADAPTER, seed=31)
    sim = simulate(cfg, d / "sim")
    (d / "bc.txt").write_text(f"S1\t{'N' * cfg.umi_len}{ADAPTER.lower()}\n")

    checkout(sim["reads"], d / "bc.txt", d / "co_fq")
    checkout(sim["reads"], d / "bc.txt", d / "co_mig", mig=True)
    fq = {
        "refine": refine(d / "co_fq" / "S1.fq.gz", d / "re_fq"),
        "assemble": assemble(d / "re_fq" / "S1.fq.gz", d / "as_fq"),
    }
    mig = {
        "refine": refine(d / "co_mig" / "S1.000.mig", d / "re_mig"),
        "assemble": assemble(d / "re_mig" / "S1.000.mig", d / "as_mig"),
    }
    return d, fq, mig


def test_the_library_actually_needed_correcting(chain):
    """Without merges this whole comparison would pass vacuously."""
    _, fq, _ = chain
    assert fq["refine"]["merged"] > 50, "the simulated library carries no barcode errors to correct"


def test_refine_reports_the_same_numbers_from_buckets(chain):
    _, fq, mig = chain
    a, b = fq["refine"], mig["refine"]
    for key in ("reads", "barcodes", "merged", "merged_reads", "merged_by_payload", "molecules"):
        assert b[key] == a[key], key
    assert b["estimated_error"] == pytest.approx(a["estimated_error"], rel=1e-9)
    assert b["error_from_children"] == pytest.approx(a["error_from_children"], rel=1e-9)
    assert b["mig_paths"], "refine wrote no buckets"


def test_the_barcode_table_is_the_same_table(chain):
    """The audit trail. A `.mig` record has no OX:Z:, so this file IS the record of the merges."""
    d, _, _ = chain
    assert (d / "re_fq" / "S1.barcodes.tsv").read_text() == (
        d / "re_mig" / "S1.barcodes.tsv"
    ).read_text()


def test_no_read_is_lost_when_a_correction_crosses_a_bucket(chain):
    """The failure mode the re-partition exists for.

    A corrected barcode is a different key, and a key decides its bucket. Copying a bucket through
    unchanged would leave those reads addressed by their old bucket, where the next stage looks for
    a different range of the key space -- so they would be grouped with strangers or dropped.
    """
    _, fq, mig = chain
    assert mig["assemble"]["reads"] == fq["assemble"]["reads"]
    assert mig["assemble"]["groups"] == fq["assemble"]["groups"]
    assert mig["assemble"]["molecules"] == fq["assemble"]["molecules"]


def test_the_consensus_records_are_identical(chain):
    d, _, _ = chain
    a = gzip.open(d / "as_fq" / "S1.consensus.fq.gz", "rt").read()
    b = gzip.open(d / "as_mig" / "S1.consensus.fq.gz", "rt").read()
    assert a == b
    assert (d / "as_fq" / "S1.mig.tsv").read_text() == (d / "as_mig" / "S1.mig.tsv").read_text()


def test_the_barcode_quality_survives_the_format(chain):
    """`.mig` v2 carries the barcode's own quality, and refine's posterior weighs it.

    v1 stored only the minimum over the barcode, which is not the same evidence: the posterior
    reads the quality AT THE POSITION THAT DIFFERS. If the format dropped it, the two routes would
    disagree on the borderline merges -- which is exactly what the equality above would catch, and
    this asserts the mechanism rather than the consequence.
    """
    from migec import _core

    d, _, _ = chain
    f = _core.MigFile(str(d / "co_mig" / "S1.000.mig"))
    header = f.header
    assert header["format_version"] == _core.MIG_FORMAT_VERSION
    records = f.read_all()
    assert records, "the bucket is empty"
    assert all(len(r.qual_umi) == header["umi_len"] for r in records[:100])
