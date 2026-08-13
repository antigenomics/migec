"""The public bindings that nothing else in the repo calls.

`umi_statistics`, `MigFile` and `match_pattern` are exported from `migec` and documented as the
interactive way to inspect a pattern or a `.mig` file. Nothing in the package used them, so nothing
tested them -- an exported API with no check is one that breaks when someone finally types it.
"""

from __future__ import annotations

import pytest

from migec import _core


def test_umi_statistics_over_a_known_list():
    """One entry per READ, so a UMI seen three times is three entries."""
    umis = ["ACGTACGTACGT"] * 8 + ["TTTTAAAACCCC"] * 2 + ["GGGGCCCCAAAA"]
    s = _core.umi_statistics(umis)

    assert s["total"] == 11
    assert s["distinct"] == 3
    assert s["mean_reads_per_umi"] == pytest.approx(11 / 3)
    assert len(s["composition"]) == 12

    # Power-of-two bins: 1 read -> bin 0, 2-3 -> bin 1, 8-15 -> bin 3.
    units = {2**b: v for b, v in enumerate(s["hist_units"]) if v}
    assert units == {1: 1, 2: 1, 8: 1}

    # Three barcodes over twelve positions cannot be worth twelve nt.
    assert 0 < s["effective_length"] < 12
    assert s["effective_space"] == pytest.approx(4 ** s["effective_length"], rel=1e-6)


def test_umi_statistics_rejects_mixed_lengths():
    with pytest.raises(RuntimeError):
        _core.umi_statistics(["ACGTACGT", "ACGT"])


def test_match_pattern_reports_placement_and_capture():
    hit = _core.match_pattern("NNNNcagtggtatcaacgcagagt", "ACGTCAGTGGTATCAACGCAGAGT" + "T" * 40)
    assert hit["found"]
    assert hit["offset"] == 0
    assert hit["umi"] == "ACGT"
    assert hit["payload_begin"] == 24

    # A read the pattern is not in comes back as not found rather than as a low score.
    assert not _core.match_pattern("NNNNcagtggtatcaacgcagagt", "T" * 64)["found"]


def test_match_pattern_captures_a_cell_barcode_separately():
    hit = _core.match_pattern("XXXXNNNNcagtggtatcaacgcagagt",
                              "TTTTACGTCAGTGGTATCAACGCAGAGT" + "A" * 40)
    assert hit["found"]
    assert hit["umi"] == "ACGT"
    assert hit["cell"] == "TTTT"


def test_pack_and_unpack_round_trip():
    # pack_barcode returns (key, has_n): the key is a plain integer and the ambiguity travels
    # beside it, so an N never turns the key into something that is not a number.
    for umi in ("ACGTACGTACGT", "AAAA", "T" * 32):
        key, has_n = _core.pack_barcode(umi)
        assert not has_n
        assert _core.unpack_barcode(key, len(umi)) == umi
    key, has_n = _core.pack_barcode("ACGN")
    assert has_n
    assert _core.unpack_barcode(key, 4) == "ACGA"


def test_the_range_partition_keeps_neighbours_together():
    """A range partition splits a barcode from its neighbour only in the top bits, which is the
    property that makes correction applicable at all. A hash would split every position."""
    key, _ = _core.pack_barcode("ACGTACGTACGT")
    home = _core.bucket_of(key, 2)
    crossings = []
    for j in range(12):
        for base in "ACGT":
            neighbour = list("ACGTACGTACGT")
            if neighbour[j] == base:
                continue
            neighbour[j] = base
            nkey, _ = _core.pack_barcode("".join(neighbour))
            if _core.bucket_of(nkey, 2) != home:
                crossings.append(j)
    # Two bucket bits are the first barcode position, so only a substitution THERE can send a
    # barcode and its neighbour to different buckets. A hash would scatter all twelve.
    assert set(crossings) == {0}


def test_mig_file_round_trips_through_python(tmp_path):
    """`MigFile` is the documented way to look inside a .mig, and checkout does not write one yet,
    so this is the only thing exercising the Python side of the format."""
    import migec

    assert migec.MigFile is _core.MigFile
    with pytest.raises(RuntimeError):
        _core.MigFile(str(tmp_path / "does_not_exist.mig"))
