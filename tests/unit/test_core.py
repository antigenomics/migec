"""Tests for the native extension as seen from Python."""

from __future__ import annotations

import gzip

import pytest

from tests.conftest import requires_core

pytestmark = requires_core


def test_versions_agree():
    import migec
    from migec import _core

    assert migec.__version__ == _core.__version__


def test_pack_unpack_round_trip():
    from migec import pack_barcode, unpack_barcode

    packed, has_n = pack_barcode("ACGTACGTACGT")
    assert not has_n
    assert unpack_barcode(packed, 12) == "ACGTACGTACGT"

    packed, has_n = pack_barcode("ACGN")
    assert has_n
    assert unpack_barcode(packed, 4) == "ACGA"


def test_packed_order_matches_lexicographic_order():
    # The range partition and the on-disk sort both depend on this.
    from migec import pack_barcode

    words = ["AAAA", "AACA", "ACGT", "CAAA", "GTTT", "TTTT"]
    packed = [pack_barcode(w)[0] for w in words]
    assert packed == sorted(packed)


def test_bucket_of_is_order_preserving():
    from migec import bucket_of, pack_barcode

    lo = bucket_of(pack_barcode("AAAAAAAAAAAA")[0], 4)
    hi = bucket_of(pack_barcode("TTTTTTTTTTTT")[0], 4)
    assert lo == 0
    assert hi == 15
    assert bucket_of(pack_barcode("ACGT")[0], 0) == 0


def test_reverse_complement_reverses_quality():
    from migec import reverse_complement

    seq, qual = reverse_complement("ACGTN", "12345")
    assert seq == "NACGT"
    assert qual == "54321"


def test_count_fastq(tmp_path):
    from migec import count_fastq

    p = tmp_path / "x.fq.gz"
    with gzip.open(p, "wt") as fh:
        for i in range(50):
            fh.write(f"@r{i}\nACGT\n+\nIIII\n")
    assert count_fastq(str(p)) == 50


def test_malformed_fastq_raises(tmp_path):
    from migec import MigecError, count_fastq

    p = tmp_path / "bad.fq"
    p.write_text("this is not a fastq\n")
    with pytest.raises(MigecError):
        count_fastq(str(p))
