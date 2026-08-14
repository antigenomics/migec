"""assemble, end to end against the simulator's truth.

The C++ unit tests cover the consensus itself (tests/cpp/test_consensus.cpp). What is checked here
is the part only an end-to-end run can check: that the reads of a molecule find each other through
the range partition, that the emitted consensus is the molecule's true sequence, and that the
quality cap survives the whole pipeline.
"""

from __future__ import annotations

import gzip

import pytest

from migec.assemble import format_report, run
from migec.checkout import run as checkout_run

from ._sim import SimConfig, simulate

from tests.conftest import requires_core

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def build(tmp_path, **kwargs):
    """Simulate, check out, and assemble. Returns (summary, consensus by UMI, truth by UMI)."""
    cfg = SimConfig(adapter=ADAPTER, **kwargs)
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")
    summary = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "asm")

    consensus = {}
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                umi = next(f[5:] for f in line.split() if f.startswith("RX:Z:"))
                depth = int(next(f[5:] for f in line.split() if f.startswith("cD:i:")))
            elif i % 4 == 1:
                consensus.setdefault(umi, []).append((line.rstrip("\n"), depth))

    truth = {}
    with open(sim["truth_consensus"]) as fh:
        umi = None
        for line in fh:
            if line.startswith(">"):
                umi = next(f[4:] for f in line.split() if f.startswith("umi="))
            else:
                truth.setdefault(umi, []).append(line.rstrip("\n"))
    return summary, consensus, truth


def test_every_molecule_comes_back(tmp_path):
    summary, consensus, truth = build(
        tmp_path, n_molecules=400, n_clones=6, coverage=8.0, umi_error=0.0
    )
    # A collision puts two molecules on one barcode, so groups <= molecules simulated.
    assert summary["groups"] == len(consensus)
    assert summary["molecules"] >= summary["groups"]
    assert set(consensus) <= set(truth)
    assert len(consensus) > 350


def test_the_consensus_is_the_molecule(tmp_path):
    _, consensus, truth = build(
        tmp_path, n_molecules=300, n_clones=4, coverage=8.0, seq_error=5e-3, umi_error=0.0
    )
    # Stratified by depth, because that is what the claim is about: a consensus over one read is
    # that read, and averaging it in hides the thing being measured.
    err = {}
    for umi, seqs in consensus.items():
        if len(truth[umi]) != 1 or len(seqs) != 1:
            continue  # a collision: the group is two molecules, and there is no single truth
        (a, depth), b = seqs[0], truth[umi][0]
        n = min(len(a), len(b))
        bucket = "deep" if depth >= 5 else "shallow"
        mm, bases = err.get(bucket, (0, 0))
        err[bucket] = (mm + sum(a[i] != b[i] for i in range(n)), bases + n)

    deep_mm, deep_bases = err["deep"]
    assert deep_bases > 10_000
    # The M1 gate: per-base error at or below 1e-5 once a molecule has five reads.
    assert deep_mm / deep_bases <= 1e-5
    # And the shallow tail is where the residual actually lives, which is the point of reporting
    # the MIG size histogram rather than one number.
    shallow_mm, shallow_bases = err["shallow"]
    assert shallow_mm / shallow_bases > deep_mm / max(deep_bases, 1)


def test_no_quality_exceeds_the_rt_floor(tmp_path):
    cfg = SimConfig(adapter=ADAPTER, n_molecules=200, coverage=20.0, mean_qual=40)
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")
    run(tmp_path / "co" / "S1.fq.gz", tmp_path / "asm", rt_floor=1e-3)  # -> Q30 cap

    worst = 0
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 3:
                worst = max(worst, max(ord(c) - 33 for c in line.rstrip("\n")))
    assert worst <= 30


def test_deep_groups_reach_the_cap(tmp_path):
    summary, _, _ = build(tmp_path, n_molecules=150, coverage=30.0, coverage_cv=0.3)
    # Deep, clean molecules should sit at the cap rather than short of it.
    assert summary["mean_quality"] > 0.9 * summary["quality_cap"]
    assert summary["mean_consensus_error"] < 1e-4


def test_tags_and_the_per_molecule_table(tmp_path):
    build(tmp_path, n_molecules=120, coverage=6.0)
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        header = fh.readline().rstrip("\n")
    name, comment = header[1:].split(" ", 1)
    assert name.startswith("S1.")
    # TAB between tags: bwa -C and minimap2 -y copy the comment verbatim into the SAM record.
    tags = dict(t.split(":Z:") if ":Z:" in t else t.split(":i:") for t in comment.split("\t"))
    assert set(tags) == {"RX", "BC", "MI", "cD"}
    assert tags["BC"] == "S1"
    assert int(tags["cD"]) >= 1

    rows = (tmp_path / "asm" / "S1.mig.tsv").read_text().splitlines()
    assert rows[0].split("\t")[:2] == ["cell", "umi"]
    assert len(rows) - 1 >= 100


def test_a_file_without_rx_tags_is_reported_not_assumed(tmp_path):
    plain = tmp_path / "plain.fq"
    plain.write_text("@r0\nACGTACGT\n+\nIIIIIIII\n")
    summary = run(plain, tmp_path / "asm")
    assert summary["reads"] == 1
    assert summary["reads_without_umi"] == 1
    assert summary["groups"] == 0
    assert "no RX tag" in format_report(summary)


def test_two_runs_with_different_umi_lengths_are_refused(tmp_path):
    mixed = tmp_path / "mixed.fq"
    mixed.write_text(
        "@r0 RX:Z:ACGTACGTACGT\nACGT\n+\nIIII\n@r1 RX:Z:ACGTAC\nACGT\n+\nIIII\n"
    )
    with pytest.raises(RuntimeError, match="two runs concatenated"):
        run(mixed, tmp_path / "asm")


def test_the_report_says_when_everything_is_a_singleton(tmp_path):
    summary, _, _ = build(tmp_path, n_molecules=300, coverage=1.0, coverage_cv=0.1)
    assert "seen once" in format_report(summary)


def test_the_same_umi_in_two_cells_is_two_molecules(tmp_path):
    """A UMI is only ever unique inside the compartment it was added to. Two cells reusing one
    barcode is the normal case, not an error, and grouping on the UMI alone would merge them."""
    reads = tmp_path / "cells.fq"
    payload_a, payload_b = "A" * 60, "G" * 60
    lines = []
    for cell, payload in (("AAAACCCCGGGGTTTT", payload_a), ("TTTTGGGGCCCCAAAA", payload_b)):
        for i in range(6):
            lines.append(
                f"@r{cell}{i} RX:Z:ACGTACGTACGT\tCB:Z:{cell}\tBC:Z:S1\n"
                f"{payload}\n+\n{'I' * len(payload)}\n"
            )
    reads.write_text("".join(lines))

    summary = run(reads, tmp_path / "asm", sample_id="S1")
    assert summary["cell_length"] == 16
    assert summary["groups"] == 2
    assert summary["molecules"] == 2

    seqs = {}
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                cb = next(f[5:] for f in line.split() if f.startswith("CB:Z:"))
            elif i % 4 == 1:
                seqs[cb] = line.rstrip("\n")
    assert seqs["AAAACCCCGGGGTTTT"] == payload_a
    assert seqs["TTTTGGGGCCCCAAAA"] == payload_b


def _tiled(molecule, starts, length, umi="ACGTACGTACGT"):
    out = []
    for i, s in enumerate(starts):
        frag = molecule[s : s + length]
        out.append(f"@t{i} RX:Z:{umi}\tBC:Z:S1\n{frag}\n+\n{'I' * len(frag)}\n")
    return "".join(out)


def test_contig_mode_assembles_tiled_reads_into_one_contig(tmp_path):
    import random

    rng = random.Random(3)
    molecule = "".join(rng.choice("ACGT") for _ in range(300))
    reads = tmp_path / "tiled.fq"
    # Overlapping by 40 nt at each step: one molecule, one contig, 300 nt long.
    reads.write_text(_tiled(molecule, [0, 30, 60, 90, 120, 150, 180, 210], 90))

    summary = run(reads, tmp_path / "asm", sample_id="S1", contig=True)
    assert summary["groups"] == 1
    assert summary["molecules"] == 1
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        fh.readline()
        assert fh.readline().rstrip("\n") == molecule


def test_contig_mode_never_bridges_a_gap(tmp_path):
    """Two islands of reads sharing a barcode but no sequence are two contigs. A single consensus
    over them would assert 100 nt that no read covers."""
    import random

    rng = random.Random(4)
    molecule = "".join(rng.choice("ACGT") for _ in range(400))
    reads = tmp_path / "gapped.fq"
    reads.write_text(_tiled(molecule, [0, 30, 60, 250, 280, 310], 90))

    summary = run(reads, tmp_path / "asm", sample_id="S1", contig=True)
    assert summary["groups"] == 1
    assert summary["groups_fragmented"] == 1
    assert summary["contigs"] == 2
    assert summary["molecules"] == 2

    lengths = []
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                lengths.append(len(line.rstrip("\n")))
    assert sorted(lengths) == [150, 150]


def test_amplicon_mode_is_the_default_and_does_not_place_reads(tmp_path):
    import random

    rng = random.Random(5)
    molecule = "".join(rng.choice("ACGT") for _ in range(300))
    reads = tmp_path / "tiled.fq"
    reads.write_text(_tiled(molecule, [0, 30, 60, 90], 90))

    summary = run(reads, tmp_path / "asm", sample_id="S1")
    assert summary["molecules"] == 1
    assert summary["contigs"] == 0
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        fh.readline()
        # Left-anchored, so the consensus is 90 nt of disagreeing reads -- which is why --contig
        # exists and why it is not the default.
        assert len(fh.readline().rstrip("\n")) == 90


def test_a_saturated_barcode_is_reported_not_assumed(tmp_path):
    """A 4 nt UMI over 200 molecules cannot tag them distinctly, and the report has to say so
    rather than presenting 256 groups as 256 molecules."""
    import random

    rng = random.Random(6)
    reads = tmp_path / "short.fq"
    lines = []
    for i in range(600):
        umi = "".join(rng.choice("ACGT") for _ in range(4))
        payload = "".join(rng.choice("ACGT") for _ in range(60))
        lines.append(f"@r{i} RX:Z:{umi}\tBC:Z:S1\n{payload}\n+\n{'I' * 60}\n")
    reads.write_text("".join(lines))

    summary = run(reads, tmp_path / "asm", sample_id="S1")
    assert summary["expected_molecules_per_group"] > 1.5
    assert "longer UMI" in format_report(summary)


def test_a_shallow_library_is_reported_not_thresholded(tmp_path):
    """1-3 reads per UMI, which is the common case rather than the exotic one: bulk repertoire
    profiling and shallow 3' GEX both look like this. Nothing may be silently dropped, and the
    report has to say that the UMI is buying counting rather than error correction."""
    summary, consensus, _ = build(
        tmp_path, n_molecules=8_000, n_clones=50, coverage=1.6, coverage_cv=0.35, umi_len=12
    )
    hist = {b["min_reads"]: b["groups"] for b in summary["coverage"]}
    assert hist[1] > 0.5 * summary["groups"], "this fixture is meant to be singleton-dominated"

    # Every read is accounted for: --min-reads defaults to 1, so nothing is thresholded away.
    assert summary["reads_dropped"] == 0
    assert summary["molecules"] == summary["groups"]
    assert sum(len(v) for v in consensus.values()) == summary["molecules"]

    report = format_report(summary)
    assert "seen once" in report
    assert "counting here, not error correction" in report


def test_the_split_threshold_is_inert_when_no_group_is_deep_enough(tmp_path):
    """A pair of columns can carry at most log10 C(n, n/2), so a 1-3 read library cannot reach
    8.68 however the reads disagree. Zero splits there is correct, and it is not evidence that the
    library has no doublets -- it is evidence there was nothing to test."""
    summary, _, _ = build(
        tmp_path, n_molecules=4_000, coverage=1.5, coverage_cv=0.3, seq_error=1e-2
    )
    assert summary["groups_split"] == 0


def test_shallow_consensus_error_is_the_read_error(tmp_path):
    """A consensus over one read IS that read, so its posterior is the read's own reported error
    and no better. The reported number must show that rather than one borrowed from the deep case,
    which is why the MIG size histogram is printed next to it."""
    shallow, _, _ = build(
        tmp_path / "shallow", n_molecules=4_000, coverage=1.2, coverage_cv=0.2, mean_qual=30
    )
    deep, _, _ = build(
        tmp_path / "deep", n_molecules=1_000, coverage=20.0, coverage_cv=0.3, mean_qual=30
    )
    # A single Q30 read cannot do better than 1e-3, and that is what comes out.
    assert shallow["mean_consensus_error"] == pytest.approx(1e-3, rel=0.5)
    assert deep["mean_consensus_error"] < 0.1 * shallow["mean_consensus_error"]
    # ...and the floor still caps what is claimed, however good or bad the input.
    assert shallow["mean_quality"] <= shallow["quality_cap"]
    assert deep["mean_quality"] <= deep["quality_cap"]


# ---------------------------------------------------------------------------------------------
# Counting mode (`--fast`): the modal exact sequence, and the best quality any read of it carried.


def test_fast_mode_takes_the_modal_sequence_and_the_best_quality(tmp_path):
    """Three reads of one molecule, two agreeing. The majority string wins whole, and each of its
    bases carries the best quality the reads that voted for it reported -- never a quality taken
    from the read that disagreed."""
    reads = tmp_path / "modal.fq"
    reads.write_text(
        "@a RX:Z:ACGTACGTACGT\tBC:Z:S1\nACGTACGT\n+\n5I5IIIII\n"
        "@b RX:Z:ACGTACGTACGT\tBC:Z:S1\nACGTACGT\n+\nI5IIIIII\n"
        "@c RX:Z:ACGTACGTACGT\tBC:Z:S1\nTTTTTTTT\n+\nIIIIIIII\n"
    )
    summary = run(reads, tmp_path / "asm", sample_id="S1", fast=True, rt_floor=1e-6)
    assert summary["groups"] == 1
    assert summary["molecules"] == 1
    # All three reads counted towards the molecule; two of them carried what was emitted.
    assert summary["mean_support"] == pytest.approx(2 / 3)

    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        _, seq, _, qual = (fh.readline().rstrip("\n") for _ in range(4))
    assert seq == "ACGTACGT"
    # Position 0 is Q40 in read b, position 1 is Q40 in read a: the per-base max over the two
    # reads that carried this sequence. The Q20s ('5') are gone, and the Q40s of the read that
    # voted for TTTTTTTT never entered.
    assert [ord(c) - 33 for c in qual] == [40] * 8
    assert summary["reads"] == 3


def test_fast_mode_still_obeys_the_rt_floor(tmp_path):
    reads = tmp_path / "one.fq"
    reads.write_text("@a RX:Z:ACGTACGTACGT\tBC:Z:S1\nACGTACGT\n+\n" + "I" * 8 + "\n")
    summary = run(reads, tmp_path / "asm", sample_id="S1", fast=True, rt_floor=1e-3)
    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        qual = fh.readlines()[3].rstrip("\n")
    # Q40 reported, Q30 floor: an error made before amplification is in every read, and the fast
    # path is no more entitled to claim past it than the full one.
    assert max(ord(c) - 33 for c in qual) <= 30
    assert summary["quality_cap"] == pytest.approx(30.0)


def test_fast_mode_counts_the_same_molecules_as_the_full_path(tmp_path):
    """What counting mode must preserve is the count. The sequence may differ from the full
    consensus on a group whose reads disagree -- that is the trade -- but a molecule must not
    appear or vanish."""
    cfg = SimConfig(adapter=ADAPTER, n_molecules=800, n_clones=5, coverage=6.0, seq_error=5e-3)
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")

    full = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "full")
    fast = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "fast", fast=True)

    assert fast["groups"] == full["groups"]
    assert fast["reads"] == full["reads"]
    # Never split, by construction: the linkage model is the column model, which fast mode skips.
    assert fast["molecules"] == fast["groups"]
    assert fast["molecules"] <= full["molecules"]
    assert "counting mode" in format_report(fast)


def test_fast_mode_does_not_correct_errors_and_says_so(tmp_path):
    """The trade, measured: at 5e-3 per base and 8 reads the full path removes essentially every
    sequencing error, and the fast path keeps whatever the majority string carried.

    `umi_error=0` so that the two paths are compared on the consensus alone -- a barcode error
    puts reads of one molecule under two keys and would charge both paths for it equally, which
    is a different measurement (docs/refine.rst)."""
    cfg = SimConfig(
        adapter=ADAPTER, n_molecules=600, n_clones=4, coverage=8.0, seq_error=5e-3, umi_error=0.0
    )
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")

    truth = {}
    with open(sim["truth_consensus"]) as fh:
        umi = None
        for line in fh:
            if line.startswith(">"):
                umi = next(f[4:] for f in line.split() if f.startswith("umi="))
            else:
                truth.setdefault(umi, []).append(line.rstrip("\n"))

    def error_rate(out_dir, **kw):
        """Per-base error against truth, over molecules seen five times or more.

        Stratified, because a consensus over one read IS that read in both modes: averaging the
        singletons in would measure the MIG size distribution rather than the two algorithms.
        """
        run(tmp_path / "co" / "S1.fq.gz", out_dir, **kw)
        mm = bases = 0
        with gzip.open(out_dir / "S1.consensus.fq.gz", "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    umi = next(f[5:] for f in line.split() if f.startswith("RX:Z:"))
                    depth = int(next(f[5:] for f in line.split() if f.startswith("cD:i:")))
                elif i % 4 == 1 and depth >= 5 and len(truth.get(umi, [])) == 1:
                    a, b = line.rstrip("\n"), truth[umi][0]
                    n = min(len(a), len(b))
                    mm += sum(a[j] != b[j] for j in range(n))
                    bases += n
        return mm / bases

    full = error_rate(tmp_path / "full")
    fast = error_rate(tmp_path / "fast", fast=True)
    # The M1 gate for the column model, and what majority-voting whole strings costs against it.
    assert full <= 1e-5
    assert fast > 10 * full


def test_coverage_is_capped_for_the_consensus_but_never_for_the_count(tmp_path):
    """10x's rule -- past ~10,000 reads a barcode adds nothing to the consensus and costs time and
    memory -- applied to the reads that are consensed only. The molecule's depth is the other half
    of what this pipeline produces, and capping it would flatten the abundance of exactly the
    most-amplified molecules."""
    over = 50
    reads = tmp_path / "deep.fq"
    payload = "ACGT" * 15
    with open(reads, "w") as fh:
        for i in range(10_000 + over):
            fh.write(f"@r{i} RX:Z:ACGTACGTACGT\tBC:Z:S1\n{payload}\n+\n{'I' * len(payload)}\n")

    summary = run(reads, tmp_path / "asm", sample_id="S1")
    assert summary["groups"] == 1
    assert summary["groups_capped"] == 1
    assert summary["reads_over_cap"] == over
    assert summary["max_reads_per_group"] == 10_000
    assert "still counted" in format_report(summary)

    with gzip.open(tmp_path / "asm" / "S1.consensus.fq.gz", "rt") as fh:
        header = fh.readline()
    assert "cD:i:10050" in header
    row = (tmp_path / "asm" / "S1.mig.tsv").read_text().splitlines()[1].split("\t")
    assert row[5] == "10050"
