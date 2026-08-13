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
