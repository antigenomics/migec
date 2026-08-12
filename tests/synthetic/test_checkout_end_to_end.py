"""checkout end to end, on reads built with the real MIGEC barcode table."""

from __future__ import annotations

import gzip
import random

import pytest

from tests.conftest import requires_core

pytestmark = requires_core

# misc/barcodes.txt from the MIGEC repository, verbatim.
BARCODES = """\
S1\taaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S2\taaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S3\taaGCCcagtggtatcaacgcagagtNNNNtNNNNtNNNN
S4\taaGGTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
"""
TAGS = {"S1": "ACT", "S2": "AGA", "S3": "GCC", "S4": "GGT"}
ADAPTER = "cagtggtatcaacgcagagt".upper()


def make_reads(path, n_per_sample=200, reads_per_umi=8, seed=0):
    """Reads with a known sample tag, a known 12 nt UMI, and a known payload."""
    rng = random.Random(seed)
    truth = []
    with gzip.open(path, "wt") as fh:
        i = 0
        for sample, tag in TAGS.items():
            for _ in range(n_per_sample // reads_per_umi):
                umi = "".join(rng.choice("ACGT") for _ in range(12))
                payload = "".join(rng.choice("ACGT") for _ in range(60))
                for _ in range(reads_per_umi):
                    seq = ("AA" + tag + ADAPTER + umi[:4] + "T" + umi[4:8] + "T" + umi[8:]
                           + payload)
                    fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
                    truth.append((f"r{i}", sample, umi, payload))
                    i += 1
    return truth


def test_checkout_assigns_extracts_and_trims(tmp_path):
    from migec.checkout import run

    reads = tmp_path / "in.fq.gz"
    truth = make_reads(reads)
    (tmp_path / "barcodes.txt").write_text(BARCODES)

    summary = run(reads, tmp_path / "barcodes.txt", tmp_path / "out")

    assert summary["total"] == len(truth)
    assert summary["assigned"] == len(truth), "every read carries a valid tag"
    assert summary["unmatched"] == 0
    assert summary["ambiguous"] == 0

    # Reads land in the right sample file, the UMI reaches the header, and the payload is exactly
    # what was planted -- no adapter, no tag, no UMI left in the sequence.
    by_read = {t[0]: t for t in truth}
    seen = 0
    for s in summary["samples"]:
        with gzip.open(tmp_path / "out" / f"{s['sample_id']}.fq.gz", "rt") as fh:
            lines = fh.read().splitlines()
        for j in range(0, len(lines), 4):
            name, _, comment = lines[j][1:].partition(" ")
            _, sample, umi, payload = by_read[name]
            assert sample == s["sample_id"]
            assert lines[j + 1] == payload
            assert f"RX:Z:{umi}" in comment
            assert f"BC:Z:{sample}" in comment
            assert "\t" in comment  # SAM-conformant separator, required by bwa -C
            seen += 1
    assert seen == len(truth)


def test_umi_statistics_are_reported_per_sample(tmp_path):
    from migec.checkout import run

    reads = tmp_path / "in.fq.gz"
    make_reads(reads, n_per_sample=800, reads_per_umi=8, seed=3)
    (tmp_path / "barcodes.txt").write_text(BARCODES)
    summary = run(reads, tmp_path / "barcodes.txt", tmp_path / "out")

    for s in summary["samples"]:
        assert s["umi_length"] == 12
        assert s["umis"] == 100
        assert s["mean_reads_per_umi"] == pytest.approx(8.0)
        assert s["over_sequenced"]
        # Random UMIs: close to 2 bits per position, so information content near zero.
        assert s["total_entropy"] == pytest.approx(24.0, abs=3.0)
        assert s["effective_length"] == pytest.approx(12.0, abs=1.5)
        assert len(s["composition"]) == 12
        # All reads sit in MIGs of 8.
        assert s["reads_in_migs_ge5"] == pytest.approx(1.0)


def test_qc_tables_are_written_and_parseable(tmp_path):
    from migec.checkout import run

    reads = tmp_path / "in.fq.gz"
    make_reads(reads, seed=5)
    (tmp_path / "barcodes.txt").write_text(BARCODES)
    run(reads, tmp_path / "barcodes.txt", tmp_path / "out")
    out = tmp_path / "out"

    summary_rows = (out / "checkout.summary.tsv").read_text().splitlines()
    assert len(summary_rows) == 5  # header + 4 samples

    cov = (out / "checkout.coverage.tsv").read_text().splitlines()
    assert cov[0].split("\t") == ["sample_id", "mig_size", "reads", "units"]
    assert len(cov) > 1

    comp = (out / "checkout.umi_composition.tsv").read_text().splitlines()
    assert comp[0].split("\t")[:2] == ["sample_id", "position"]
    assert len(comp) == 1 + 4 * 12  # header + 4 samples x 12 UMI positions

    assert (out / "checkout.json").exists()


def test_a_skewed_umi_reports_reduced_effective_length(tmp_path):
    """The number that matters: a UMI can be 12 nt and worth far less."""
    from migec.checkout import run

    rng = random.Random(11)
    reads = tmp_path / "in.fq.gz"
    with gzip.open(reads, "wt") as fh:
        for i in range(2000):
            # Only the last 4 positions vary; the first 8 are fixed.
            umi = "AAAACCCC" + "".join(rng.choice("ACGT") for _ in range(4))
            payload = "".join(rng.choice("ACGT") for _ in range(60))
            seq = "AAACT" + ADAPTER + umi[:4] + "T" + umi[4:8] + "T" + umi[8:] + payload
            fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    (tmp_path / "barcodes.txt").write_text(BARCODES)

    summary = run(reads, tmp_path / "barcodes.txt", tmp_path / "out")
    s = next(x for x in summary["samples"] if x["sample_id"] == "S1")
    assert s["umi_length"] == 12
    assert s["effective_length"] == pytest.approx(4.0, abs=0.5)
    assert s["effective_space"] == pytest.approx(256.0, rel=0.3)
    # Information content is where the skew shows up on a logo: 8 fixed positions x 2 bits.
    assert s["total_information"] == pytest.approx(16.0, abs=1.0)


def test_unmatched_reads_are_counted_and_optionally_written(tmp_path):
    from migec.checkout import run

    reads = tmp_path / "in.fq.gz"
    with gzip.open(reads, "wt") as fh:
        for i in range(50):
            seq = "GATTACA" * 12
            fh.write(f"@junk{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    (tmp_path / "barcodes.txt").write_text(BARCODES)

    summary = run(reads, tmp_path / "barcodes.txt", tmp_path / "out", write_unmatched=True)
    assert summary["assigned"] == 0
    assert summary["unmatched"] == 50
    with gzip.open(tmp_path / "out" / "unmatched.fq.gz", "rt") as fh:
        assert len(fh.read().splitlines()) == 50 * 4


def test_trim_none_keeps_the_read_whole(tmp_path):
    from migec.checkout import run

    reads = tmp_path / "in.fq.gz"
    truth = make_reads(reads, n_per_sample=8, reads_per_umi=8, seed=9)
    (tmp_path / "barcodes.txt").write_text(BARCODES)
    run(reads, tmp_path / "barcodes.txt", tmp_path / "out", trim="none")

    with gzip.open(tmp_path / "out" / "S1.fq.gz", "rt") as fh:
        lines = fh.read().splitlines()
    # Untrimmed reads still start with the adapter region.
    assert lines[1].startswith("AAACT" + ADAPTER)
    assert len(lines[1]) > len(truth[0][3])


def make_pairs(r1_path, r2_path, n=400, seed=11, flip_every=3):
    """Paired reads where every `flip_every`-th pair carries the tag on R2 instead of R1."""
    rng = random.Random(seed)
    truth = []
    with gzip.open(r1_path, "wt") as f1, gzip.open(r2_path, "wt") as f2:
        for i in range(n):
            sample = list(TAGS)[i % 4]
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(60))
            mate = "".join(rng.choice("ACGT") for _ in range(60))
            tagged = ("AA" + TAGS[sample] + ADAPTER + umi[:4] + "T" + umi[4:8] + "T" + umi[8:]
                      + payload)
            flipped = i % flip_every == 0
            a, b = (mate, tagged) if flipped else (tagged, mate)
            f1.write(f"@p{i}\n{a}\n+\n{'I' * len(a)}\n")
            f2.write(f"@p{i}\n{b}\n+\n{'I' * len(b)}\n")
            truth.append((f"p{i}", sample, umi, payload, mate, flipped))
    return truth


def test_paired_checkout_normalises_orientation(tmp_path):
    from migec.checkout import run

    r1, r2 = tmp_path / "R1.fq.gz", tmp_path / "R2.fq.gz"
    truth = make_pairs(r1, r2)
    (tmp_path / "barcodes.txt").write_text(BARCODES)
    summary = run(r1, tmp_path / "barcodes.txt", tmp_path / "out", reads2=r2)

    assert summary["assigned"] == len(truth)
    assert summary["paired"]
    # A pair whose tag sits on R2 must be swapped, not discarded: an amplicon library sequenced in
    # both orientations otherwise loses those reads at consensus without anything reporting it.
    assert summary["normalised"] == sum(1 for t in truth if t[5])

    by_read = {t[0]: t for t in truth}
    seen = 0
    for s in summary["samples"]:
        with gzip.open(tmp_path / "out" / f"{s['sample_id']}_R1.fq.gz", "rt") as fh:
            a = fh.read().splitlines()
        with gzip.open(tmp_path / "out" / f"{s['sample_id']}_R2.fq.gz", "rt") as fh:
            b = fh.read().splitlines()
        assert len(a) == len(b)
        for j in range(0, len(a), 4):
            name = a[j][1:].split(" ")[0]
            _, sample, umi, payload, mate, _flipped = by_read[name]
            assert sample == s["sample_id"]
            assert a[j + 1] == payload, "R1 always carries the tag after normalisation"
            assert b[j + 1] == mate, "the mate is passed through whole"
            # Both mates carry the barcode, or nothing downstream can group the pair.
            assert f"RX:Z:{umi}" in a[j] and f"RX:Z:{umi}" in b[j]
            seen += 1
    assert seen == len(truth)


def test_output_does_not_depend_on_thread_count(tmp_path):
    from migec.checkout import run

    reads = tmp_path / "in.fq.gz"
    make_reads(reads, n_per_sample=600, seed=13)
    (tmp_path / "barcodes.txt").write_text(BARCODES)

    digests = []
    for threads in (1, 4):
        out = tmp_path / f"t{threads}"
        summary = run(reads, tmp_path / "barcodes.txt", out, threads=threads)
        assert summary["threads"] == threads
        digests.append((tmp_path / f"t{threads}" / "S1.fq.gz").read_bytes())
    assert digests[0] == digests[1], "checkout output must not depend on -t"


def test_resource_use_is_reported(tmp_path):
    from migec.checkout import format_report, run

    reads = tmp_path / "in.fq.gz"
    make_reads(reads, seed=17)
    (tmp_path / "barcodes.txt").write_text(BARCODES)
    summary = run(reads, tmp_path / "barcodes.txt", tmp_path / "out")

    assert summary["wall_seconds"] > 0
    assert summary["reads_per_second"] > 0
    assert summary["peak_rss_bytes"] > 0
    assert 0 < summary["umi_memory_bytes"] < summary["peak_rss_bytes"]
    report = format_report(summary)
    assert "reads/s" in report and "peak RSS" in report
