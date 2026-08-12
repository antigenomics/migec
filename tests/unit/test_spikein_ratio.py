"""The spike-in metric is the project's headline claim, so the code computing it needs its own
test -- on synthetic reads where the answer is known by construction."""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
spec = importlib.util.spec_from_file_location("spikein_ratio", ROOT / "scripts" / "spikein_ratio.py")
sr = importlib.util.module_from_spec(spec)
sys.modules["spikein_ratio"] = sr
spec.loader.exec_module(sr)

TRUTH = sr.load_truth(ROOT / "tests" / "data" / "migec_spikein_truth.tsv")
EHEB = TRUTH["EHEB"]
V1 = TRUTH["EHEB-V1"]
V2 = TRUTH["EHEB-V2"]


def sub(seq: str, pos: int, base: str) -> str:
    return seq[:pos] + base + seq[pos + 1 :]


def write_reads(path, junction_counts, flank_left="CCCCCCCCCC", flank_right="GGGGGGGGGG"):
    with gzip.open(path, "wt") as fh:
        i = 0
        for junction, n in junction_counts.items():
            for _ in range(n):
                seq = flank_left + junction + flank_right
                fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")
                i += 1


def test_truth_table_is_self_consistent():
    # The declared substitution counts must reproduce from the sequences themselves.
    assert sr.hamming(EHEB, V1) == 1
    assert sr.hamming(EHEB, V2) == 2
    assert len(EHEB) == len(V1) == len(V2) == 48
    assert EHEB.startswith("TGT") and EHEB.endswith("TGG")  # junction, anchors included


def test_ratios_on_a_constructed_library(tmp_path):
    err1 = sub(EHEB, 20, "A" if EHEB[20] != "A" else "C")   # 1 sub, not V1
    err2 = sub(err1, 30, "A" if err1[30] != "A" else "C")   # 2 subs, not V2
    reads = tmp_path / "r.fq.gz"
    write_reads(reads, {EHEB: 10000, V1: 500, V2: 40, err1: 100, err2: 10})

    counts = sr.collect_junctions(str(reads), EHEB)
    r = sr.ratios(counts, TRUTH)

    assert r["EHEB"] == 10000
    assert r["EHEB-V1"] == 500
    assert r["EHEB-V2"] == 40
    assert r["Err1"] == 100
    assert r["Err2"] == 10
    assert r["V1/Err1"] == 5.0
    assert r["V2/Err2"] == 4.0
    assert r["EHEB/Err1"] == 100.0


def test_known_variants_are_excluded_from_the_error_set(tmp_path):
    """The metric is meaningless if V1 is allowed to be its own competitor."""
    reads = tmp_path / "r.fq.gz"
    write_reads(reads, {EHEB: 1000, V1: 500, V2: 40})
    r = sr.ratios(sr.collect_junctions(str(reads), EHEB), TRUTH)
    assert r["Err1"] == 0          # V1 is real, so there is no competing 1-sub error
    assert r["V1/Err1"] is None    # reported as undefined rather than as infinity


def test_reverse_complemented_libraries_are_found(tmp_path):
    """Experiment 1 carries this clone reverse complemented; missing it reports a false zero."""
    reads = tmp_path / "r.fq.gz"
    rc = sr.revcomp(EHEB)
    with gzip.open(reads, "wt") as fh:
        for i in range(300):
            seq = "TTTTTTTTTT" + rc + "AAAAAAAAAA"
            fh.write(f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n")

    counts = sr.collect_junctions(str(reads), EHEB)
    assert counts[EHEB] == 300  # counted in the forward orientation


def test_consensus_improves_the_ratio(tmp_path):
    """The shape of the claim: collapsing an error cloud that consensus would remove raises
    V1/Err1, while the real variant is untouched."""
    err1 = sub(EHEB, 20, "A" if EHEB[20] != "A" else "C")

    raw = tmp_path / "raw.fq.gz"
    write_reads(raw, {EHEB: 10000, V1: 500, err1: 400})
    r_raw = sr.ratios(sr.collect_junctions(str(raw), EHEB), TRUTH)

    # After consensus the PCR error, spread thinly across many molecules, mostly disappears while
    # the real variant keeps its molecules.
    cons = tmp_path / "cons.fq.gz"
    write_reads(cons, {EHEB: 1200, V1: 60, err1: 2})
    r_cons = sr.ratios(sr.collect_junctions(str(cons), EHEB), TRUTH)

    assert r_raw["V1/Err1"] == 1.25
    assert r_cons["V1/Err1"] == 30.0
    assert r_cons["V1/Err1"] > r_raw["V1/Err1"]
    lo, hi = sr.PUBLISHED["migec"]["V1/Err1"]
    assert lo <= r_cons["V1/Err1"] <= hi
