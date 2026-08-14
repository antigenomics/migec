"""A BAM carrying `RX` is the same input as the FASTQ it came from.

The stage under test is `refine`, deliberately: it opens its input THREE times in one call, so it
is exactly what a named-pipe implementation would fail. `assemble` follows, because its output is
barcode-ordered and therefore comparable byte for byte whatever the record order was.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess

import pytest

# Never: BEFORE the `migec` imports below, and not `pytestmark` alone. A module-scope
# import of a package whose extension is missing raises at COLLECTION, which pytest
# reports as an error rather than a skip.
pytest.importorskip("migec._core", reason="the C++ extension is not built: run `bash setup.sh`")

from migec.assemble import run as assemble_run
from migec.bam import _has_references, is_alignment
from migec.checkout import run as checkout_run
from migec.refine import run as refine_run

from ._sim import SimConfig, simulate

from tests.conftest import requires_core

pytestmark = [
    requires_core,
    pytest.mark.skipif(shutil.which("samtools") is None, reason="samtools is not on PATH"),
]

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def _corpus(tmp_path):
    """A checkout-tagged FASTQ, and the same records as an unaligned BAM."""
    sim = simulate(
        SimConfig(adapter=ADAPTER, n_molecules=2_000, n_clones=50, coverage=5.0, umi_len=12,
                  umi_error=2e-3, seed=11),
        tmp_path / "sim",
    )
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")
    fastq = tmp_path / "co" / "S1.fq.gz"
    ubam = tmp_path / "S1.bam"
    # `-T '*'` is samtools' own inverse of what `samtools fastq -T` writes: every SAM-style tag in
    # the comment becomes a real tag. This is the round trip, so if it needed a flag we do not set
    # on the way out, this is where it would show.
    subprocess.run(
        ["samtools", "import", "-T", "*", "-s", str(fastq), "-o", str(ubam)], check=True,
    )
    return fastq, ubam


def _records(path):
    with gzip.open(path, "rt") as fh:
        return fh.read()


def test_a_bam_is_recognised_and_a_fastq_is_not(tmp_path):
    fastq, ubam = _corpus(tmp_path)
    assert is_alignment(ubam)
    # The gzip magic is the same for both, so the check has to reach the BAM block's own magic.
    assert not is_alignment(fastq)
    assert not is_alignment(tmp_path / "bc.txt")


def test_refine_and_assemble_agree_field_for_field_and_byte_for_byte(tmp_path):
    fastq, ubam = _corpus(tmp_path)

    from_fastq = refine_run(fastq, tmp_path / "rf_fq", sample_id="S1")
    from_bam = refine_run(ubam, tmp_path / "rf_bam", sample_id="S1")

    ignored = {"input", "wall_seconds", "table_seconds", "correct_seconds", "rewrite_seconds",
               "peak_rss_bytes", "threads"}
    for key, value in from_fastq.items():
        if key in ignored:
            continue
        assert from_bam[key] == value, key
    assert from_bam["input"].endswith("S1.bam")
    # The audit trail is one row per barcode in barcode order, so it compares directly.
    assert ((tmp_path / "rf_bam" / "S1.barcodes.tsv").read_text()
            == (tmp_path / "rf_fq" / "S1.barcodes.tsv").read_text())

    a_fq = assemble_run(tmp_path / "rf_fq" / "S1.fq.gz", tmp_path / "as_fq", sample_id="S1")
    a_bam = assemble_run(tmp_path / "rf_bam" / "S1.fq.gz", tmp_path / "as_bam", sample_id="S1")
    assert a_bam["molecules"] == a_fq["molecules"]
    assert (_records(tmp_path / "as_bam" / "S1.consensus.fq.gz")
            == _records(tmp_path / "as_fq" / "S1.consensus.fq.gz"))


def test_assemble_takes_the_bam_directly(tmp_path):
    fastq, ubam = _corpus(tmp_path)
    direct = assemble_run(ubam, tmp_path / "as_bam", sample_id="S1")
    plain = assemble_run(fastq, tmp_path / "as_fq", sample_id="S1")
    assert direct["molecules"] == plain["molecules"]
    assert (_records(tmp_path / "as_bam" / "S1.consensus.fq.gz")
            == _records(tmp_path / "as_fq" / "S1.consensus.fq.gz"))


def _consensus_through(reads, tmp_path, name, reads2=None):
    """checkout -> refine -> assemble, returning the consensus FASTQ's text and the summary."""
    d = tmp_path / name
    d.mkdir()
    (d / "bc.txt").write_text(f"S1\t^NNNNNNNNNNNN{ADAPTER}\n")
    checkout_run(reads, d / "bc.txt", d / "co", reads2=reads2)
    refine_run(d / "co" / "S1_R1.fq.gz", d / "rf", sample_id="S1")
    summary = assemble_run(d / "rf" / "S1.fq.gz", d / "as", sample_id="S1")
    return _records(d / "as" / "S1.consensus.fq.gz"), summary


def test_a_paired_bam_carries_both_mates_and_collate_keeps_them_together(tmp_path):
    """The motivating case: one file, two mates, and the pairing has to survive the conversion."""
    from .test_paired_pipeline import paired_corpus

    src = tmp_path / "src"
    src.mkdir()
    paired_corpus(src)
    ubam = tmp_path / "paired.bam"
    subprocess.run(
        ["samtools", "import", "-1", str(src / "R1.fq.gz"), "-2", str(src / "R2.fq.gz"),
         "-o", str(ubam)], check=True,
    )
    # The same records with a reference in the header, which is what sends `as_fastq` down the
    # collate branch. That branch is the one a real aligned BAM takes, and mispairing there is the
    # failure nothing downstream could detect -- so it is compared against the unpaired-header
    # conversion of the identical reads.
    header = subprocess.run(
        ["samtools", "view", "-H", str(ubam)], capture_output=True, text=True, check=True,
    ).stdout
    (tmp_path / "hdr.sam").write_text(header + "@SQ\tSN:chr1\tLN:1000\n")
    aligned = tmp_path / "aligned.bam"
    with open(aligned, "wb") as fh:
        subprocess.run(
            ["samtools", "reheader", str(tmp_path / "hdr.sam"), str(ubam)], stdout=fh, check=True,
        )

    # Pin the branch, or this test would pass while only ever exercising one of the two.
    assert not _has_references(ubam)
    assert _has_references(aligned)

    from_fastq, s_fastq = _consensus_through(
        src / "R1.fq.gz", tmp_path, "fq", reads2=src / "R2.fq.gz"
    )
    from_ubam, s_ubam = _consensus_through(ubam, tmp_path, "ubam")
    from_aligned, s_aligned = _consensus_through(aligned, tmp_path, "aligned")

    assert s_ubam["molecules"] == s_fastq["molecules"]
    assert s_aligned["molecules"] == s_fastq["molecules"]
    assert from_ubam == from_fastq
    assert from_aligned == from_fastq


def test_a_umi_in_a_separate_index_read_reaches_the_same_molecules(tmp_path):
    """The motivating case, end to end: R1 plus an I1 holding the UMI, and no `checkout` at all.

    A capture, exome or ctDNA kit reads the UMI on the index, so it is never inside R1. One
    `samtools import` puts it in `RX` and the answer has to be the one `checkout` would have
    reached from the un-split reads -- otherwise the entry point is a different pipeline wearing
    the same name.
    """
    sim = simulate(
        SimConfig(adapter=ADAPTER, n_molecules=2_000, n_clones=50, coverage=5.0, umi_len=12,
                  umi_error=2e-3, seed=11),
        tmp_path / "sim",
    )
    r1, i1 = tmp_path / "R1.fq", tmp_path / "I1.fq"
    with gzip.open(sim["reads"], "rt") as fh, open(r1, "w") as fr, open(i1, "w") as fi:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].split()[0]
            elif i % 4 == 1:
                seq = line.strip()
            elif i % 4 == 3:
                qual = line.strip()
                fr.write(f"@{name}\n{seq[12:]}\n+\n{qual[12:]}\n")
                fi.write(f"@{name}\n{seq[:12]}\n+\n{qual[:12]}\n")

    tagged = tmp_path / "tagged.bam"
    # Never: `--barcode-tag RX`. Left at its default, `samtools import` writes an index read into
    # `BC`, which is the SAMPLE barcode -- refine would then find no RX and refuse the file.
    subprocess.run(
        ["samtools", "import", "-0", str(r1), "--i1", str(i1),
         "--barcode-tag", "RX", "--quality-tag", "QX", "-o", str(tagged)], check=True,
    )
    from_index = refine_run(tagged, tmp_path / "rf_idx", sample_id="S1")

    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")
    from_checkout = refine_run(tmp_path / "co" / "S1.fq.gz", tmp_path / "rf_co", sample_id="S1")

    for key in ("reads", "barcodes", "merged", "molecules"):
        assert from_index[key] == from_checkout[key], key
    assert from_index["molecules"] == 2_000


def test_a_bam_without_rx_is_refused_rather_than_reporting_no_molecules(tmp_path):
    sim = simulate(
        SimConfig(adapter=ADAPTER, n_molecules=200, n_clones=10, coverage=3.0, umi_len=12, seed=2),
        tmp_path / "sim",
    )
    bare = tmp_path / "bare.bam"
    subprocess.run(
        ["samtools", "import", "-s", str(sim["reads"]), "-o", str(bare)], check=True,
    )
    with pytest.raises(ValueError, match="no RX tag"):
        refine_run(bare, tmp_path / "rf")
    # ...but checkout is exactly the stage that extracts one, so it takes the same file.
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    summary = checkout_run(bare, tmp_path / "bc.txt", tmp_path / "co")
    assert summary["assigned"] > 0
    assert summary["input"].endswith("bare.bam")
