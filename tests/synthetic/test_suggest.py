"""`migec suggest`: recovering the barcode layout from the per-cycle base composition.

The signal is a 1/4 PWM trace. A cycle the synthesiser mixed shows all four bases near 25%; a
cycle that is a fixed adapter base shows one near 100%; biological payload is neither. Segmenting
on that gives back the pattern the library was built with, which beats trusting a protocol
description written for the bench rather than for the file.
"""

from __future__ import annotations

import pytest

from tests.conftest import requires_core
from tests.synthetic._sim import SimConfig, simulate

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def _suggest(tmp_path, cfg, **kw):
    from migec import _core

    sim = simulate(cfg, tmp_path / "sim")
    return _core.suggest(str(sim["reads"]), **kw), sim


def test_the_pattern_is_recovered_from_the_composition(tmp_path):
    cfg = SimConfig(n_molecules=3000, umi_len=12, coverage=6.0, adapter=ADAPTER, seed=1)
    s, sim = _suggest(tmp_path, cfg, cycles=40)
    assert s["umi_length"] == 12
    assert s["pattern"].startswith("N" * 12)
    assert s["pattern"][12:].upper().startswith(ADAPTER)
    # ...and it is the pattern the simulator says it wrote.
    assert s["pattern"][: len(sim["pattern"])].upper() == sim["pattern"].upper()


def test_umi_and_constant_cycles_are_told_apart_by_their_traces(tmp_path):
    cfg = SimConfig(n_molecules=3000, umi_len=12, coverage=6.0, adapter=ADAPTER, seed=2)
    s, _ = _suggest(tmp_path, cfg, cycles=40)
    umi = [c for c in s["cycles"] if c["cycle"] < 12]
    const = [c for c in s["cycles"] if 12 <= c["cycle"] < 32]
    # A UMI cycle is four flat lines at 25%: near-zero deviation, ~2 bits, collision ~1/4.
    assert all(c["deviation"] < 0.06 for c in umi), [c["deviation"] for c in umi]
    assert all(c["entropy"] > 1.95 for c in umi)
    assert all(abs(c["collision"] - 0.25) < 0.01 for c in umi)
    # A constant cycle is one base at ~100%: maximal deviation, ~0 bits, collision ~1.
    assert all(c["deviation"] > 0.6 for c in const)
    assert all(c["entropy"] < 0.2 for c in const)
    assert all(c["collision"] > 0.9 for c in const)


def test_a_skewed_synthesiser_mix_is_still_read_as_umi(tmp_path):
    # An oligo synthesiser does not deliver 25% of each base. 20/30/30/20 is ordinary and must not
    # be mistaken for payload -- it is exactly the case `effective_length` exists to quantify.
    cfg = SimConfig(
        n_molecules=3000, umi_len=12, coverage=6.0, adapter=ADAPTER, seed=3,
        umi_base_freqs=(0.20, 0.30, 0.30, 0.20),
    )
    s, _ = _suggest(tmp_path, cfg, cycles=40)
    assert s["umi_length"] == 12, s["pattern"]


def test_the_segments_cover_every_profiled_cycle(tmp_path):
    cfg = SimConfig(n_molecules=2000, umi_len=12, coverage=5.0, adapter=ADAPTER, seed=4)
    s, _ = _suggest(tmp_path, cfg, cycles=50)
    assert s["segments"][0]["begin"] == 0
    assert s["segments"][-1]["end"] == 50
    for a, b in zip(s["segments"], s["segments"][1:]):
        assert a["end"] == b["begin"], "segments must tile the profile with no gap or overlap"
    kinds = [x["kind"] for x in s["segments"]]
    assert kinds[0] == "umi" and kinds[1] == "constant"


def test_a_read_with_no_umi_says_so_rather_than_inventing_one(tmp_path):
    # This is the paired-end case: the tag is on the other mate, and the honest answer is to say
    # nothing was found rather than to label payload as a barcode.
    path = tmp_path / "payload.fq"
    body = "ACGTTTGCAAGGCCTTAACCGGTTACGTACGTACGT"
    path.write_text("".join(f"@r{i}\n{body}\n+\n{'I' * len(body)}\n" for i in range(2000)))
    from migec import _core

    s = _core.suggest(str(path), cycles=30)
    assert s["umi_length"] == 0
    assert "other mate" in s["note"]


def test_the_suggested_pattern_actually_checks_out(tmp_path):
    """The point of the pattern is to be used, so run it."""
    from migec.checkout import run
    from migec.suggest import run as suggest_run

    cfg = SimConfig(n_molecules=3000, umi_len=12, coverage=6.0, adapter=ADAPTER, seed=5)
    sim = simulate(cfg, tmp_path / "sim")
    s = suggest_run(sim["reads"], tmp_path / "sug", cycles=40)
    (tmp_path / "bc.txt").write_text(f"S1\t{s['pattern']}\n")
    summary = run(sim["reads"], tmp_path / "bc.txt", tmp_path / "out")
    assert summary["assigned"] / summary["total"] > 0.98, summary
    assert summary["samples"][0]["umi_length"] == 12
    assert (tmp_path / "sug" / "suggest.cycles.tsv").exists()
    assert (tmp_path / "sug" / "suggest.segments.tsv").exists()


def test_the_report_renders(tmp_path):
    from migec.suggest import format_report, run

    cfg = SimConfig(n_molecules=1500, umi_len=12, coverage=5.0, adapter=ADAPTER, seed=6)
    sim = simulate(cfg, tmp_path / "sim")
    text = format_report(run(sim["reads"], cycles=40))
    assert "UMI" in text and "pattern" in text
    assert "4^12" in text and f"{4**12:,}" in text


def test_an_empty_file_is_an_error(tmp_path):
    from migec import _core

    path = tmp_path / "empty.fq"
    path.write_text("")
    with pytest.raises(RuntimeError, match="no reads"):
        _core.suggest(str(path))
