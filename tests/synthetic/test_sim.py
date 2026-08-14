"""The simulator is a test asset, so it needs its own tests: if its truth is wrong, every accuracy
number measured against it is wrong in the same direction and nothing downstream will notice."""

from __future__ import annotations

import gzip
from collections import Counter

import pytest

from tests.synthetic._sim import SimConfig, simulate


def test_reads_match_truth_table(tmp_path):
    cfg = SimConfig(n_molecules=200, n_clones=5, coverage=6.0, seed=1)
    info = simulate(cfg, tmp_path)

    with gzip.open(info["reads"], "rt") as fh:
        names = [line[1:].split()[0] for i, line in enumerate(fh) if i % 4 == 0]
    with open(info["truth_reads"]) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh][1:]

    assert len(names) == len(rows) == info["n_reads"]
    assert names == [r[0] for r in rows]


def test_molecule_sizes_sum_to_read_count(tmp_path):
    cfg = SimConfig(n_molecules=300, coverage=4.0, seed=2)
    info = simulate(cfg, tmp_path)
    with open(info["truth_molecules"]) as fh:
        sizes = [int(line.split("\t")[3]) for line in fh.readlines()[1:]]
    assert sum(sizes) == info["n_reads"]
    assert len(sizes) == cfg.n_molecules
    # The low-coverage regime the retention rule is about must actually be populated: a Poisson
    # at mean 4 would leave almost nothing at 1-2 reads, which would make the rule untestable.
    assert sum(1 for s in sizes if s <= 2) > 0


def test_umi_errors_are_injected_and_recorded(tmp_path):
    cfg = SimConfig(n_molecules=400, coverage=8.0, umi_error=0.01, seed=3)
    info = simulate(cfg, tmp_path)
    with open(info["truth_reads"]) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh.readlines()[1:]]
    mismatched = sum(1 for r in rows if r[3] != r[4])
    assert mismatched > 0, "umi_error=0.01 over 12 bases must corrupt some UMIs"
    # Roughly 1 - (1-0.01)^12 = 11.4% of reads; allow a wide band, this is a smoke check.
    assert 0.03 < mismatched / len(rows) < 0.25


def test_zero_error_is_exactly_reproducible(tmp_path):
    cfg = SimConfig(n_molecules=50, coverage=5.0, seq_error=0.0, umi_error=0.0, rt_error=0.0,
                    pcr_error=0.0, seed=7)
    info = simulate(cfg, tmp_path / "a")
    cons = {}
    with open(info["truth_consensus"]) as fh:
        name = None
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                cons[name] = line.strip()

    # With no error of any kind, every read of a molecule equals that molecule's consensus.
    with gzip.open(info["reads"], "rt") as fh:
        lines = fh.read().splitlines()
    with open(info["truth_reads"]) as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh.readlines()[1:]]
    umi_len = cfg.umi_len
    for i, row in enumerate(rows):
        seq = lines[4 * i + 1]
        assert seq[:umi_len] == row[4]                      # observed UMI is the prefix
        assert seq[umi_len:] == cons[f"mol{row[1]}"]        # payload is the true consensus


def test_same_seed_same_output(tmp_path):
    a = simulate(SimConfig(n_molecules=100, seed=11), tmp_path / "a")
    b = simulate(SimConfig(n_molecules=100, seed=11), tmp_path / "b")
    assert open(a["truth_reads"]).read() == open(b["truth_reads"]).read()
    c = simulate(SimConfig(n_molecules=100, seed=12), tmp_path / "c")
    assert open(a["truth_reads"]).read() != open(c["truth_reads"]).read()


def test_collisions_appear_when_the_umi_space_is_small(tmp_path):
    # 4^4 = 256 UMIs for 400 molecules: collisions are guaranteed by the birthday bound, and the
    # correction step must not treat them as errors.
    cfg = SimConfig(n_molecules=400, umi_len=4, coverage=3.0, umi_error=0.0, seed=5)
    info = simulate(cfg, tmp_path)
    assert info["n_umi_collisions"] > 0
    assert info["n_distinct_umis"] < cfg.n_molecules

    with open(info["truth_molecules"]) as fh:
        umis = Counter(line.split("\t")[2] for line in fh.readlines()[1:])
    assert max(umis.values()) > 1


def test_variant_clones_are_one_reference_plus_point_variants(tmp_path):
    """`variant_af` is what makes a variant CALLER comparable, so its truth has to be exact."""
    cfg = SimConfig(n_molecules=4000, n_clones=4, coverage=5.0, seq_len=120,
                    variant_af=0.05, seed=3)
    info = simulate(cfg, tmp_path)

    reference = "".join(
        line.strip() for line in open(info["reference"]).read().splitlines()
        if not line.startswith(">")
    )
    rows = [line.split("\t") for line in open(info["truth_variants"]).read().splitlines()[1:]]
    assert len(rows) == info["n_variants"] == cfg.n_clones - 1

    clones = {}
    name = None
    for line in open(info["clones"]).read().splitlines():
        if line.startswith(">"):
            name = line[1:]
            clones[name] = ""
        else:
            clones[name] += line.strip()
    assert clones["clone0"] == reference

    for i, (chrom, pos, ref_base, alt, af) in enumerate(rows, start=1):
        p = int(pos) - 1
        assert chrom == "clone0"
        assert reference[p] == ref_base            # the ref base is the reference's base
        assert clones[f"clone{i}"][p] == alt       # and the variant clone carries the alt
        # ...and differs from the reference nowhere else, or "the variant" would be plural.
        assert sum(a != b for a, b in zip(clones[f"clone{i}"], reference)) == 1
        assert float(af) == cfg.variant_af

    # The molecules actually drawn follow the requested fraction, which is what a called allele
    # fraction is compared against.
    per_clone = Counter(
        int(line.split("\t")[1])
        for line in open(info["truth_molecules"]).read().splitlines()[1:]
    )
    for i in range(1, cfg.n_clones):
        assert per_clone[i] / cfg.n_molecules == pytest.approx(cfg.variant_af, rel=0.3)
