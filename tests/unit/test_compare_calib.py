"""The grouping-accuracy metric used to compare against Calib."""

from __future__ import annotations

import gzip
import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "compare_calib.py"
spec = importlib.util.spec_from_file_location("compare_calib", SCRIPT)
cc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cc)


def write_truth(path, pairs):
    with open(path, "w") as fh:
        fh.write("read_id\tmolecule_id\tclone_id\n")
        for r, m in pairs:
            fh.write(f"{r}\t{m}\tC1\n")


def test_a_perfect_partition_scores_one(tmp_path):
    truth_path = tmp_path / "t.tsv"
    write_truth(truth_path, [("r0", "m0"), ("r1", "m0"), ("r2", "m1"), ("r3", "m1")])
    truth = cc.read_truth(truth_path)
    pred = {"r0": "a", "r1": "a", "r2": "b", "r3": "b"}
    m = cc.adjusted_rand(truth, pred)
    assert m["ari"] == pytest.approx(1.0)
    assert m["reads_in_split_molecules"] == 0.0
    assert m["reads_in_merged_clusters"] == 0.0


def test_split_and_merge_are_reported_separately(tmp_path):
    # ARI alone cannot say whether a tool over-split or over-merged, and the two failures have
    # opposite costs: splitting inflates the molecule count, merging destroys real variants.
    truth_path = tmp_path / "t.tsv"
    write_truth(truth_path, [("r0", "m0"), ("r1", "m0"), ("r2", "m1"), ("r3", "m1")])
    truth = cc.read_truth(truth_path)

    over_split = cc.adjusted_rand(truth, {"r0": "a", "r1": "b", "r2": "c", "r3": "c"})
    assert over_split["reads_in_split_molecules"] == pytest.approx(0.5)
    assert over_split["reads_in_merged_clusters"] == 0.0
    assert over_split["predicted_clusters"] == 3

    over_merged = cc.adjusted_rand(truth, {"r0": "a", "r1": "a", "r2": "a", "r3": "a"})
    assert over_merged["reads_in_merged_clusters"] == pytest.approx(1.0)
    assert over_merged["reads_in_split_molecules"] == 0.0
    assert over_merged["predicted_clusters"] == 1


def test_calib_cluster_file_is_parsed(tmp_path):
    # Nine columns: cluster_id, node_id, read_id, then name/seq/qual per mate. The mate suffix has
    # to come off or nothing lines up with the truth file.
    path = tmp_path / "x.cluster"
    path.write_text(
        "0\t0\t0\tr0/1\tACGT\tIIII\tr0/2\tACGT\tIIII\n"
        "0\t1\t1\tr1/1\tACGT\tIIII\tr1/2\tACGT\tIIII\n"
        "1\t2\t2\tr2/1\tACGT\tIIII\tr2/2\tACGT\tIIII\n"
    )
    assert cc.read_calib(path) == {"r0": "0", "r1": "0", "r2": "1"}


def test_migec_groups_come_from_the_header_tags(tmp_path):
    path = tmp_path / "S1.fq.gz"
    with gzip.open(path, "wt") as fh:
        for name, umi in (("r0", "AAAACCCCGGGG"), ("r1", "AAAACCCCGGGG"), ("r2", "TTTTGGGGCCCC")):
            fh.write(f"@{name} RX:Z:{umi}\tQX:Z:{'I' * 12}\tBC:Z:S1\nACGT\n+\nIIII\n")
    groups = cc.read_migec([path])
    assert groups["r0"] == groups["r1"] != groups["r2"]
    assert groups["r0"] == "S1:AAAACCCCGGGG"


def test_untagged_input_is_rejected_rather_than_scored(tmp_path):
    # Silently scoring every read into one group would report a spectacular merge rate rather than
    # the actual mistake, which is that the file never went through checkout.
    path = tmp_path / "raw.fq"
    path.write_text("@r0\nACGT\n+\nIIII\n")
    with pytest.raises(SystemExit):
        cc.read_migec([path])
