"""The i7 x i5 contingency table, and the one contamination a per-sample count cannot see.

On a patterned flowcell a free index primer can prime a neighbouring cluster, so a molecule from
one sample is read carrying its own i7 and another sample's i5. The read then lands in that other
sample and looks exactly like one of its reads -- there is nothing in the sequence, the barcode or
the count to say otherwise. The only evidence is the index pair in the instrument's own header,
which is why this costs nothing to measure and is worth measuring on every run.

The library here is built with a known hopped fraction, so the estimate has a truth to be right
about rather than a plausible-looking number.
"""

from __future__ import annotations

import gzip
import random

import pytest

pytest.importorskip("migec._core", reason="the C++ extension is not built: run `bash setup.sh`")

from tests.conftest import requires_core  # noqa: E402

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"
# Two dual-indexed samples, as an ordinary sheet declares them.
DECLARED = [("ATCACG", "CGTGAT"), ("CGATGT", "ACATCG")]
HOPPED_FRACTION = 0.02


def _corpus(path, n_reads=20000, seed=7, hop=HOPPED_FRACTION):
    """A dual-indexed run with a known share of reads carrying a swapped i5."""
    rng = random.Random(seed)
    hopped = 0
    with gzip.open(path, "wt", compresslevel=1) as fh:
        for i in range(n_reads):
            sample = i % 2
            i7, i5 = DECLARED[sample]
            if rng.random() < hop:
                # The hop: this sample's i7 with the OTHER sample's i5.
                i5 = DECLARED[1 - sample][1]
                hopped += 1
            umi = "".join(rng.choice("ACGT") for _ in range(8))
            payload = "".join(rng.choice("ACGT") for _ in range(60))
            seq = umi + ADAPTER + payload
            fh.write(f"@r{i} 1:N:0:{i7}+{i5}\n{seq}\n+\n{'I' * len(seq)}\n")
    return hopped


def _sheet(path):
    path.write_text(f"S1\t{'N' * 8}{ADAPTER.lower()}\n")


@pytest.fixture(scope="module")
def run_summary(tmp_path_factory):
    from migec.checkout import run

    d = tmp_path_factory.mktemp("hopping")
    hopped = _corpus(d / "reads.fq.gz")
    _sheet(d / "bc.txt")
    return d, hopped, run(d / "reads.fq.gz", d / "bc.txt", d / "out")


def test_the_table_is_read_off_the_header(run_summary):
    _, _, s = run_summary
    ix = s["index_hopping"]
    assert ix["estimable"]
    assert ix["i7_indices"] == 2
    assert ix["i5_indices"] == 2
    # Two declared combinations and two that nobody ordered: the full 2x2.
    assert ix["declared_pairs"] == 2
    assert ix["hopped_pairs"] == 2


def test_the_rate_is_the_rate_that_was_injected(run_summary):
    """The number, against the truth. A plausible rate computed the wrong way is the failure mode."""
    _, hopped, s = run_summary
    ix = s["index_hopping"]
    assert ix["reads_hopped"] == hopped
    assert ix["rate"] == pytest.approx(HOPPED_FRACTION, abs=0.005)


def test_a_single_indexed_run_is_not_estimable(tmp_path):
    """Never: unmeasurable is not zero.

    With one index there are no combinations, so nothing can be off-diagonal. Reporting 0.0% there
    would say the run is clean when in fact the question cannot be asked of it.
    """
    from migec.checkout import run

    with gzip.open(tmp_path / "reads.fq.gz", "wt") as fh:
        for i in range(200):
            seq = "ACGTACGT" + ADAPTER + "A" * 60
            fh.write(f"@r{i} 1:N:0:ATCACG\n{seq}\n+\n{'I' * len(seq)}\n")
    _sheet(tmp_path / "bc.txt")
    s = run(tmp_path / "reads.fq.gz", tmp_path / "bc.txt", tmp_path / "out")
    assert not s["index_hopping"]["estimable"]
    assert s["index_hopping"]["rate"] == 0.0


def test_a_header_with_no_index_says_nothing(tmp_path):
    from migec.checkout import run

    with gzip.open(tmp_path / "reads.fq.gz", "wt") as fh:
        for i in range(200):
            seq = "ACGTACGT" + ADAPTER + "A" * 60
            fh.write(f"@read_{i}\n{seq}\n+\n{'I' * len(seq)}\n")
    _sheet(tmp_path / "bc.txt")
    s = run(tmp_path / "reads.fq.gz", tmp_path / "bc.txt", tmp_path / "out")
    assert s["index_hopping"]["pairs"] == []
    assert not s["index_hopping"]["estimable"]


def test_unmatched_reads_are_counted_too(tmp_path):
    """Never: restricting the table to ASSIGNED reads hides the population it exists to measure.

    A hopped read whose in-line barcode does not match any pattern is still evidence that indices
    hopped on this run.
    """
    from migec.checkout import run

    with gzip.open(tmp_path / "reads.fq.gz", "wt") as fh:
        for i in range(100):
            # No adapter: nothing will match the sheet's pattern.
            seq = "".join(random.Random(i).choice("ACGT") for _ in range(80))
            fh.write(f"@r{i} 1:N:0:ATCACG+CGTGAT\n{seq}\n+\n{'I' * len(seq)}\n")
    _sheet(tmp_path / "bc.txt")
    s = run(tmp_path / "reads.fq.gz", tmp_path / "bc.txt", tmp_path / "out")
    assert s["assigned"] == 0
    assert sum(row["reads"] for row in s["index_hopping"]["pairs"]) == 100


def test_the_table_is_written(run_summary):
    d, _, s = run_summary
    rows = (d / "out" / "checkout.index_pairs.tsv").read_text().splitlines()
    assert rows[0].split("\t") == ["i7", "i5", "reads", "share_of_i7", "share_of_i5", "declared"]
    assert len(rows) == 1 + len(s["index_hopping"]["pairs"])
    # Sorted by reads, so the declared combinations are the first rows a reader sees.
    counts = [int(r.split("\t")[2]) for r in rows[1:]]
    assert counts == sorted(counts, reverse=True)


# --- the other half of the same header: where on the flowcell ------------------------------------


def test_the_tile_map_is_read_off_the_header(tmp_path):
    """A tile is a physical patch of the flowcell, and a read count cannot show a dead one.

    Never: this is a MAP, not a total. Two lanes at the same yield and one tile of the second
    delivering a tenth of its neighbours is a run with a problem, and every other number migec
    reports is identical between that run and a healthy one.
    """
    from migec.checkout import run

    with gzip.open(tmp_path / "reads.fq.gz", "wt") as fh:
        i = 0
        for lane in (1, 2):
            for tile in (1101, 1102):
                # The second lane's second tile is starved, as a bubble would leave it.
                n = 10 if (lane, tile) == (2, 1102) else 100
                for _ in range(n):
                    umi = "".join(random.Random(i).choice("ACGT") for _ in range(8))
                    seq = umi + ADAPTER + "A" * 60
                    fh.write(
                        f"@M01:1:FC1:{lane}:{tile}:{1000 + i}:{2000 + i} 1:N:0:ATCACG+CGTGAT\n"
                        f"{seq}\n+\n{'I' * len(seq)}\n"
                    )
                    i += 1
    _sheet(tmp_path / "bc.txt")
    s = run(tmp_path / "reads.fq.gz", tmp_path / "bc.txt", tmp_path / "out")

    assert s["reads_with_coordinates"] == s["total"] == 310
    by_key = {(r["lane"], r["tile"]): r["reads"] for r in s["tiles"]}
    assert by_key == {(1, 1101): 100, (1, 1102): 100, (2, 1101): 100, (2, 1102): 10}

    rows = (tmp_path / "out" / "checkout.tiles.tsv").read_text().splitlines()
    assert rows[0].split("\t") == ["lane", "tile", "reads", "share_of_lane"]
    # The starved tile is a tenth of its lane's other tile, and the share column says so.
    starved = [r for r in rows[1:] if r.startswith("2\t1102")][0]
    assert float(starved.split("\t")[3]) == pytest.approx(10 / 110, abs=1e-6)


def test_a_header_the_instrument_did_not_write_gives_no_map(tmp_path):
    """SRA rewrites headers to `@SRR1763769.1 1/2` and the coordinates are gone for good.

    Never invent them: a run whose headers did not survive reports no map, rather than a map of
    one tile that reads as a single-tile flowcell.
    """
    from migec.checkout import run

    with gzip.open(tmp_path / "reads.fq.gz", "wt") as fh:
        for i in range(50):
            seq = "ACGTACGT" + ADAPTER + "A" * 60
            fh.write(f"@SRR1763769.{i} {i}/1\n{seq}\n+\n{'I' * len(seq)}\n")
    _sheet(tmp_path / "bc.txt")
    s = run(tmp_path / "reads.fq.gz", tmp_path / "bc.txt", tmp_path / "out")
    assert s["tiles"] == []
    assert s["reads_with_coordinates"] == 0
