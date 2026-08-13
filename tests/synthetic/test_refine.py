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
    assert rows[0] == "cell\tumi\treads\tcorrected_reads\tparent"
    assert len(rows) - 1 == s["barcodes"]
    merged = sum(1 for r in rows[1:] if r.split("\t")[4] != ".")
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


def _cell_reads(path, n_cells, per_cell, umi_error=0.0, seed=0):
    """A 10x-shaped corpus: every cell reuses the SAME UMIs, which is the case that breaks a
    UMI-only barcode table."""
    import gzip
    import random

    rng = random.Random(seed)
    umis = ["".join(rng.choice("ACGT") for _ in range(12)) for _ in range(per_cell)]
    cells = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(n_cells)]
    with gzip.open(path, "wt") as fh:
        i = 0
        for cell in cells:
            for umi in umis:
                obs = list(umi)
                if rng.random() < umi_error:
                    j = rng.randrange(12)
                    obs[j] = rng.choice([b for b in "ACGT" if b != obs[j]])
                payload = "".join(rng.choice("ACGT") for _ in range(60))
                for _ in range(4):
                    fh.write(
                        f"@r{i} RX:Z:{''.join(obs)}\tQX:Z:{'I' * 12}\tCB:Z:{cell}\t"
                        f"CY:Z:{'I' * 16}\tBC:Z:S1\n{payload}\n+\n{'I' * 60}\n"
                    )
                    i += 1
    return len(cells) * len(umis)


def test_the_same_umi_in_two_cells_is_not_corrected_away(tmp_path):
    """Every cell reuses the same 200 UMIs. Keyed on the UMI alone the table would hold 200
    barcodes for 10 cells' worth of molecules and correct them against each other."""
    reads = tmp_path / "cells.fq.gz"
    molecules = _cell_reads(reads, n_cells=10, per_cell=200, seed=1)
    s = run(reads, tmp_path / "ref")
    assert s["cell_length"] == 16
    assert s["barcodes"] == molecules == 2_000
    assert s["merged"] == 0
    assert s["molecules"] == molecules


def test_a_cell_barcode_error_is_corrected_and_audited(tmp_path):
    """The concatenated key means a substitution in the CELL barcode is a distance-1 neighbour too,
    so it is found -- and the read records what it was in OC:Z:."""
    import gzip

    reads = tmp_path / "cb.fq.gz"
    good = "ACGTACGTACGTACGT"
    bad = "ACGTACGTACGTACGA"  # distance 1 in the cell barcode
    payload = "ACGT" * 15
    with gzip.open(reads, "wt") as fh:
        for i in range(60):
            cb = bad if i == 0 else good
            fh.write(
                f"@r{i} RX:Z:AAACCCGGGTTT\tQX:Z:{'I' * 12}\tCB:Z:{cb}\tCY:Z:{'I' * 16}\t"
                f"BC:Z:S1\n{payload}\n+\n{'I' * len(payload)}\n"
            )
        # Background molecules, so the library is not one clone in one cell.
        import random

        rng = random.Random(2)
        for i in range(4_000):
            cb = "".join(rng.choice("ACGT") for _ in range(16))
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            pay = "".join(rng.choice("ACGT") for _ in range(60))
            fh.write(
                f"@b{i} RX:Z:{umi}\tQX:Z:{'I' * 12}\tCB:Z:{cb}\tCY:Z:{'I' * 16}\t"
                f"BC:Z:S1\n{pay}\n+\n{'I' * 60}\n"
            )
    s = run(reads, tmp_path / "ref")
    assert s["merged"] >= 1

    fixed = 0
    with gzip.open(tmp_path / "ref" / "S1.fq.gz", "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 != 0:
                continue
            tags = dict(f.split(":Z:") for f in line.split() if ":Z:" in f)
            if tags.get("OC") == bad:
                fixed += 1
                assert tags["CB"] == good
    assert fixed == 1


def test_the_diagnostics_say_where_the_errors_are(tmp_path):
    """The shape that matters: error children pile up at one read. A flat curve, or one rising at
    high counts, means correction is eating real molecules rather than errors."""
    build(tmp_path, n_molecules=20_000, n_clones=100, coverage=6.0, coverage_cv=0.5,
          umi_len=12, umi_error=3e-3, seed=5)
    run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")

    rows = [
        dict(zip(h, r.split("\t")))
        for h in [(tmp_path / "ref" / "S1.bins.tsv").read_text().splitlines()[0].split("\t")]
        for r in (tmp_path / "ref" / "S1.bins.tsv").read_text().splitlines()[1:]
    ]
    by_size = {int(r["min_reads"]): r for r in rows}
    assert float(by_size[1]["fraction_erroneous"]) > 0.5
    assert float(by_size[4]["fraction_erroneous"]) < 0.05
    # Simulated payloads are random, so every bin sits near two bits. The point of the column is
    # that a bin holding one sequence repeated drops well below the rest.
    for r in rows:
        if int(r["barcodes"]) > 50:
            assert 1.5 < float(r["payload_entropy_bits"]) <= 2.0


def test_the_rank_curve_is_log_spaced_and_monotone(tmp_path):
    """Cell Ranger's plot. One row per barcode would be hundreds of millions for a figure read on
    a log axis, so the ranks are log-spaced -- and the curve still has to be a curve."""
    build(tmp_path, n_molecules=20_000, coverage=6.0, coverage_cv=0.6, umi_len=12, seed=5)
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")

    lines = (tmp_path / "ref" / "S1.rank.tsv").read_text().splitlines()
    assert lines[0] == "rank\treads\tcumulative_reads\tcumulative_fraction"
    rows = [line.split("\t") for line in lines[1:]]
    assert 20 < len(rows) < 0.05 * s["molecules"], "the rank table is not log-spaced"

    ranks = [int(r[0]) for r in rows]
    reads = [int(r[1]) for r in rows]
    fractions = [float(r[3]) for r in rows]
    assert ranks == sorted(ranks)
    assert reads == sorted(reads, reverse=True), "barcodes must be ranked by depth"
    assert fractions == sorted(fractions)
    assert fractions[-1] == pytest.approx(1.0, abs=1e-6)


def _droplets(path, n_cells, cell_molecules, n_ambient, ambient_molecules, seed=0):
    """A 10x-shaped library: a few hundred real cells over a swamp of ambient barcodes."""
    import gzip
    import random

    rng = random.Random(seed)
    out, i = [], 0
    def emit(cb, n):
        nonlocal i
        for _ in range(n):
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            pay = "".join(rng.choice("ACGT") for _ in range(60))
            for _ in range(2):
                out.append(
                    f"@r{i} RX:Z:{umi}\tQX:Z:{'I' * 12}\tCB:Z:{cb}\tCY:Z:{'I' * 16}\t"
                    f"BC:Z:S1\n{pay}\n+\n{'I' * 60}\n"
                )
                i += 1
    real = []
    for _ in range(n_cells):
        cb = "".join(rng.choice("ACGT") for _ in range(16))
        real.append(cb)
        emit(cb, rng.randint(*cell_molecules))
    for _ in range(n_ambient):
        emit("".join(rng.choice("ACGT") for _ in range(16)), rng.randint(*ambient_molecules))
    with gzip.open(path, "wt") as fh:
        fh.write("".join(out))
    return set(real)


def test_cells_are_called_off_the_molecule_curve(tmp_path):
    """500 real cells at 200-600 molecules, 20,000 ambient barcodes at 1-3. OrdMag has to find the
    500 and nothing else."""
    reads = tmp_path / "d.fq.gz"
    real = _droplets(reads, 500, (200, 600), 20_000, (1, 3), seed=0)
    s = run(reads, tmp_path / "ref", expect_cells=500)

    assert s["cells_observed"] == 20_500
    assert s["cells_called"] == 500
    assert s["cell_threshold"] > 3, "the threshold must sit above the ambient barcodes"

    rows = [
        line.split("\t")
        for line in (tmp_path / "ref" / "S1.cells.tsv").read_text().splitlines()[1:]
    ]
    assert len(rows) == s["cells_observed"]
    called = {r[0] for r in rows if r[2] == "1"}
    # Exactly the cells that were put in: no ambient barcode called, no real cell missed.
    assert called == real


def test_the_knee_is_reported_next_to_the_call_not_instead_of_it(tmp_path):
    """OrdMag makes the call; the knee describes the curve. Both are reported so a disagreement is
    visible rather than resolved silently."""
    reads = tmp_path / "d.fq.gz"
    _droplets(reads, 400, (300, 900), 10_000, (1, 2), seed=1)
    s = run(reads, tmp_path / "ref", expect_cells=400)
    assert 0 < s["knee_rank"] <= s["cells_observed"]
    assert s["knee_molecules"] > 0
    # On a library with a clean separation the two agree to well within the factor of three the
    # report warns at.
    assert max(s["knee_rank"], s["cells_called"]) < 3 * min(s["knee_rank"], s["cells_called"])


def test_expect_cells_being_wrong_does_not_move_the_call_much(tmp_path):
    """OrdMag takes the 99th percentile of the top N, so over-stating N by 400x still lands the
    index inside the real cells. That robustness is the reason the rule is worth having, and it
    is worth a test rather than an assumption."""
    reads = tmp_path / "d.fq.gz"
    real = _droplets(reads, 500, (200, 600), 20_000, (1, 3), seed=0)
    tight = run(reads, tmp_path / "tight", expect_cells=500)
    loose = run(reads, tmp_path / "loose", expect_cells=200_000)
    assert tight["cells_called"] == len(real)
    assert loose["cells_called"] == pytest.approx(len(real), rel=0.02)


def test_a_disagreement_between_ordmag_and_the_knee_is_reported(tmp_path):
    """The call is a rule and the knee is what the data says on its own. When they part company
    the report has to say so rather than presenting the rule's answer alone."""
    summary = {
        "reads": 1000, "barcodes": 100, "merged": 0, "merged_reads": 0, "merged_by_payload": 0,
        "molecules": 100, "molecules_corrected": 100.0, "saturated": False,
        "estimated_error": 1e-3, "payload_clonality": 0.01, "wall_seconds": 1.0,
        "peak_rss_bytes": 1 << 20, "table_bytes": 1 << 10, "coverage": [],
        "cell_length": 16, "cells_observed": 10_000, "cells_called": 5_000,
        "molecules_in_called": 90, "cell_threshold": 2, "knee_rank": 300,
        "knee_molecules": 40, "min_posterior": 0.95,
        "whitelist": {"barcodes": 0, "exact": 0, "corrected": 0, "off_list": 0,
                      "reads_corrected": 0, "far": 0, "background_prior": 0.0},
        "suspected_residual": 0, "residual_fdr_at_one": 0.0, "mig_size_threshold": 1,
        "target_fdr": 0.05,
    }
    report = format_report(summary)
    assert "OrdMag calls" in report
    assert "16.7" in report  # 5000 / 300
    summary["cells_called"] = 320  # now they agree
    assert "OrdMag calls" not in format_report(summary)


def test_a_bulk_library_has_no_cells_to_call(tmp_path):
    build(tmp_path, n_molecules=2_000, coverage=5.0, umi_len=12, seed=8)
    s = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "ref")
    assert s["cell_length"] == 0
    assert s["cells_observed"] == 0
    assert s["cells_called"] == 0
    assert not (tmp_path / "ref" / "S1.cells.tsv").exists()


def test_a_whitelist_snaps_errors_and_refuses_strangers(tmp_path):
    """The two halves of the same decision. A barcode one substitution off a heavily-used list
    entry is a miscall; a barcode that is nobody's neighbour is a real off-list barcode and must
    be left exactly where it is."""
    import gzip
    import random

    rng = random.Random(4)
    listed = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(5_000)]
    (tmp_path / "wl.txt").write_text("\n".join(listed) + "\n")

    reads, i = [], 0
    def emit(cb, n, qual_at_end=40):
        nonlocal i
        for _ in range(n):
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            pay = "".join(rng.choice("ACGT") for _ in range(60))
            cy = "I" * 15 + chr(33 + qual_at_end)
            reads.append(
                f"@r{i} RX:Z:{umi}\tQX:Z:{'I' * 12}\tCB:Z:{cb}\tCY:Z:{cy}\tBC:Z:S1\n"
                f"{pay}\n+\n{'I' * 60}\n"
            )
            i += 1

    # 200 real cells, deeply used, from the list.
    for cb in listed[:200]:
        emit(cb, 400)
    # One of them mis-sequenced at the last base, called at Q3 -- a base the instrument itself
    # says it is unsure of. The posterior scales as (parent reads) x e/3, so a snap needs both a
    # well-used parent and a genuinely poor base; a confident base is never overridden.
    off_by_one = listed[0][:15] + ("A" if listed[0][15] != "A" else "C")
    emit(off_by_one, 20, qual_at_end=3)
    # ...and a swarm of barcodes that are nobody's neighbour: ambient, hopped, undeclared.
    strangers = []
    while len(strangers) < 3_000:
        cb = "".join(rng.choice("ACGT") for _ in range(16))
        if cb not in set(listed):
            strangers.append(cb)
    for cb in strangers:
        emit(cb, 1)

    with gzip.open(tmp_path / "r.fq.gz", "wt") as fh:
        fh.write("".join(reads))

    s = run(tmp_path / "r.fq.gz", tmp_path / "ref", cell_whitelist=tmp_path / "wl.txt")
    w = s["whitelist"]
    assert w["exact"] == 200
    assert w["barcodes"] == 200 + 1 + len(strangers)
    # The mis-sequenced barcode is snapped back...
    assert w["corrected"] == 1
    assert w["reads_corrected"] == 20
    # ...and not one of the 3,000 strangers is absorbed into a neighbour.
    assert w["off_list"] >= len(strangers)
    assert w["background_prior"] > 0

    # The snap is applied to the reads and audited, whether or not the posterior also merged it.
    fixed = 0
    with gzip.open(tmp_path / "ref" / "S1.fq.gz", "rt") as fh:
        for j, line in enumerate(fh):
            if j % 4 != 0:
                continue
            tags = dict(f.split(":Z:") for f in line.split() if ":Z:" in f)
            if tags.get("OC") == off_by_one:
                fixed += 1
                assert tags["CB"] == listed[0]
    assert fixed == 20


def test_a_whitelist_of_the_wrong_length_is_refused(tmp_path):
    import gzip

    (tmp_path / "wl.txt").write_text("ACGTACGT\nTTTTGGGG\n")
    with gzip.open(tmp_path / "r.fq.gz", "wt") as fh:
        fh.write(
            f"@r0 RX:Z:ACGTACGTACGT\tQX:Z:{'I' * 12}\tCB:Z:ACGTACGTACGTACGT\t"
            f"CY:Z:{'I' * 16}\tBC:Z:S1\nACGT\n+\nIIII\n"
        )
    with pytest.raises(RuntimeError, match="whitelist holds"):
        run(tmp_path / "r.fq.gz", tmp_path / "ref", cell_whitelist=tmp_path / "wl.txt")


def test_the_residual_fdr_finds_what_correction_left_behind(tmp_path):
    """The number that says whether the molecule count can be trusted. It must rise where
    correction is known to fail -- 1-3 reads per UMI -- and fall where it does not."""
    build(tmp_path / "shallow", n_molecules=30_000, n_clones=100, coverage=1.6,
          coverage_cv=0.4, umi_len=12, umi_error=5e-3, seed=5)
    shallow = run(tmp_path / "shallow" / "co" / "S1.fq.gz", tmp_path / "shallow" / "ref")
    build(tmp_path / "deep", n_molecules=30_000, n_clones=100, coverage=6.0,
          coverage_cv=0.4, umi_len=12, umi_error=5e-3, seed=5)
    deep = run(tmp_path / "deep" / "co" / "S1.fq.gz", tmp_path / "deep" / "ref")

    assert shallow["residual_fdr_at_one"] > 10 * deep["residual_fdr_at_one"]
    assert shallow["suspected_residual"] > 0
    # Never: Reported, never applied: every molecule survives whatever the threshold says.
    assert shallow["mig_size_threshold"] >= 1
    hist = {b["min_reads"]: b["molecules"] for b in shallow["coverage"]}
    assert hist[1] > 0, "1-read molecules must still be in the output"
    assert "REPORTED, not applied" in format_report(shallow)


def test_a_count_ratio_alone_would_report_no_residual_at_all(tmp_path):
    """The trap this estimator had to avoid. At 1-3 reads per UMI nothing is 20x anything, so a
    pure count-ratio criterion reports zero residual in exactly the regime where it is worst."""
    build(tmp_path, n_molecules=30_000, n_clones=100, coverage=1.6, coverage_cv=0.4,
          umi_len=12, umi_error=5e-3, seed=5)
    with_payload = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "with")
    counts_only = run(tmp_path / "co" / "S1.fq.gz", tmp_path / "without", use_payload=False)
    assert with_payload["suspected_residual"] > 0
    assert counts_only["suspected_residual"] == 0
