"""refine, against the simulator's truth.

The stage that decides how many molecules there were, so what is checked is the molecule count and
what it cost: how many real molecules were destroyed to get it, and whether the correction can be
audited afterwards.
"""

from __future__ import annotations

import gzip

import pytest

from migec.checkout import run as checkout_run
from migec.refine import format_report, run

from ._sim import SimConfig, simulate

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def build(tmp_path, **kwargs):
    cfg = SimConfig(adapter=ADAPTER, **kwargs)
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")
    return sim, cfg


def test_the_molecule_count_comes_back(tmp_path):
    sim, cfg = build(
        tmp_path, n_molecules=20_000, n_clones=100, coverage=6.0, coverage_cv=0.4,
        umi_len=12, umi_error=3e-3, seed=5,
    )
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    # Before correction the barcode errors inflate the count badly; after it, within a few percent.
    assert s["barcodes"] > 1.1 * cfg.n_molecules
    assert s["molecules"] == pytest.approx(cfg.n_molecules, rel=0.05)
    # The M3 gate: the estimated per-base UMI error within 20% of what was injected.
    assert s["estimated_error"] == pytest.approx(cfg.umi_error, rel=0.2)


def test_clonality_is_measured_not_assumed(tmp_path):
    """Payload agreement is worth log(1/clonality), so the clonality has to be right."""
    build(tmp_path / "diverse", n_molecules=8_000, n_clones=200, coverage=5.0, umi_len=12, seed=3)
    diverse = run(tmp_path / "diverse" / "co" / "S1.fq.gz", tmp_path / "diverse" / "ref")
    build(tmp_path / "clonal", n_molecules=8_000, n_clones=1, coverage=5.0, umi_len=12, seed=3)
    clonal = run(tmp_path / "clonal" / "co" / "S1.fq.gz", tmp_path / "clonal" / "ref")

    assert diverse["payload_clonality"] == pytest.approx(1 / 200, abs=0.01)
    assert clonal["payload_clonality"] > 0.9
    # ...and the report says what that buys, rather than leaving a bare number.
    assert "clonal" in format_report(clonal)


def test_a_corrected_read_records_what_it_was(tmp_path):
    """A correction nobody can audit is a correction nobody can check."""
    build(tmp_path, n_molecules=5_000, coverage=6.0, umi_len=12, umi_error=5e-3, seed=4)
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    assert s["merged"] > 0

    corrected, distances = 0, []
    with gzip.open(tmp_path / "ref" / "S1.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            tags = dict(f.split(":Z:") for f in line.split() if ":Z:" in f)
            if "OX" in tags:
                corrected += 1
                assert tags["OX"] != tags["RX"]
                distances.append(sum(a != b for a, b in zip(tags["OX"], tags["RX"])))
    assert corrected == s["merged_reads"]
    # Merges chain: x folds into y, and y is later folded into z, so a read can end up two
    # substitutions from where it started. Each STEP is distance 1; the chain is not, and a long
    # one would mean correction is walking away from the molecule rather than towards it.
    assert min(distances) == 1
    assert max(distances) <= 3
    assert distances.count(1) > 0.9 * len(distances)


def test_the_barcode_table_is_complete(tmp_path):
    build(tmp_path, n_molecules=3_000, coverage=5.0, umi_len=12, seed=6)
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    rows = (tmp_path / "ref" / "S1.barcodes.tsv").read_text().splitlines()
    assert rows[0] == "umi\treads\tcorrected_reads\tparent"
    assert len(rows) - 1 == s["barcodes"]
    merged = sum(1 for r in rows[1:] if r.split("\t")[3] != ".")
    assert merged == s["merged"]


def test_the_evidence_earns_its_place(tmp_path):
    """Turning the two new terms off must measurably cost something, or they are decoration."""
    build(tmp_path, n_molecules=20_000, n_clones=100, coverage=2.0, coverage_cv=0.4,
          umi_len=12, umi_error=3e-3, seed=5)
    reads = tmp_path / "co" / "S1.fq.gz"
    full = run(reads, tmp_path / "full")
    counts_only = run(reads, tmp_path / "counts", use_quality=False, use_payload=False)
    assert full["merged"] > counts_only["merged"]
    assert full["merged_by_payload"] > 0
    # ...and payload evidence is what lifts the count gates, so it is the term that matters here.
    assert counts_only["merged_by_payload"] == 0


def test_a_shallow_library_corrects_conservatively(tmp_path):
    """1-3 reads per UMI. Most barcode errors have no observable parent, so the honest outcome is
    a molecule count that is still inflated -- not one bought by deleting real molecules."""
    sim, cfg = build(
        tmp_path, n_molecules=20_000, n_clones=100, coverage=1.3, coverage_cv=0.35,
        umi_len=12, umi_error=3e-3, seed=5,
    )
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    # Never below the truth: that would mean real molecules were merged away.
    assert s["molecules"] >= 0.98 * cfg.n_molecules
    assert "reads per barcode" in format_report(s)


def test_a_file_without_rx_tags_is_refused(tmp_path):
    plain = tmp_path / "plain.fq"
    plain.write_text("@r0\nACGTACGT\n+\nIIIIIIII\n")
    with pytest.raises(RuntimeError, match="RX:Z:"):
        run(plain, tmp_path / "ref")


def test_refine_output_feeds_assemble(tmp_path):
    """The pipeline contract: checkout -> refine -> assemble, with the tags intact throughout."""
    from migec.assemble import run as assemble_run

    build(tmp_path, n_molecules=5_000, n_clones=50, coverage=6.0, umi_len=12,
          umi_error=3e-3, seed=7)
    refined = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    assembled = assemble_run(tmp_path / "ref" / "S1.fq.gz", tmp_path / "asm")
    assert assembled["groups"] == refined["molecules"]
    assert assembled["molecules"] >= assembled["groups"]
