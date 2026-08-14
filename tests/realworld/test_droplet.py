"""sc5p_v2_hs_PBMC_1k VDJ-T, a 10x droplet library, through checkout and refine.

The fixture is `ci/sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz`: 1% of the barcodes of 10x's public 5' v2
T-cell library, 18,036 reads over 1,014 cells and 2,636 (cell, UMI) barcodes (`SOURCES.md`).

**The pair is reconstructed, and here that is lossless.** The fixture is checkout's R2 output, so
the barcode lives in the `CB`/`CY` and `RX`/`QX` tags. On 5' v2 R1 is *exactly* 26 nt and holds
nothing but the barcode -- 16 nt cell + 10 nt UMI -- so concatenating the tags back reproduces the
original R1 byte for byte, and R2 was never trimmed (the pattern matches on R1). Unlike the
amplicon fixture, nothing here is re-synthesised: `checkout` sees the library's own bases and
qualities.

**What this fixture cannot show is cell calling with a knee.** Keeping 1% of the barcodes divides
molecules per cell by ~100: the full library is 305,702 molecules over 813 called cells (~376
each), the fixture is 2,628 over 1,014 (~2.6). OrdMag has nothing to find in that -- it calls
every observed cell at a threshold of 1 molecule, against a knee at 7 -- so the assertions below
are that cell calling RUNS, that its key is the whole cell barcode, and that what it counts is
molecules rather than reads. A cell/background separation belongs on the full library, not here.
"""

from __future__ import annotations

import gzip

import pytest

from tests.conftest import UMI_DATA, requires_core, requires_umi_data

from ._fixture import records

pytestmark = [requires_core, requires_umi_data]

READS = UMI_DATA / "ci" / "sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz"


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    """checkout on the reconstructed pair, then refine on R2. Under two seconds end to end."""
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run
    from migec.sheet import parse_layout, preset

    work = tmp_path_factory.mktemp("droplet")
    truth: dict[str, tuple[str, str]] = {}
    with gzip.open(work / "R1.fq.gz", "wt") as r1, gzip.open(work / "R2.fq.gz", "wt") as r2:
        for name, tags, sequence, quality in records(READS):
            truth[name] = (tags["CB"], tags["RX"])
            r1.write(f"@{name}\n{tags['CB']}{tags['RX']}\n+\n{tags['CY']}{tags['QX']}\n")
            r2.write(f"@{name}\n{sequence}\n+\n{quality}\n")

    pattern, anchored = parse_layout(preset("10x-v2")[0])
    (work / "bc.txt").write_text(f"PBMC\t{pattern}\n")
    checkout = checkout_run(work / "R1.fq.gz", work / "bc.txt", work / "co", work / "R2.fq.gz")
    # Never: refine takes R2. Trimming the pattern leaves R1 empty, so the payload-agreement term
    # would have nothing to read and the reported clonality would come out 1.0 saying so.
    refine = refine_run(work / "co" / "PBMC_R2.fq.gz", work / "ref")
    return {
        "truth": truth,
        "pattern": pattern,
        "anchored": anchored,
        "checkout": checkout,
        "sample": checkout["samples"][0],
        "refine": refine,
        "dir": work,
    }


def test_the_preset_describes_this_chemistry(pipeline):
    """`10x-v2` has to be 16 cell bases then 10 UMI bases, at the read start and nowhere else.

    The layout is purely positional: there is no constant sequence anywhere in R1, so a free scan
    has no evidence to choose an offset with and the pattern must anchor. That anchoring is what
    `parse_layout` returns, and getting it wrong is not a subtle failure -- a scored-nothing
    pattern scanned freely is refused on every read.
    """
    pattern = pipeline["pattern"]
    assert pipeline["anchored"]
    assert pattern == "X" * 16 + "N" * 10
    assert pipeline["sample"]["umi_length"] == 10
    assert pipeline["refine"]["cell_length"] == 16


def test_every_read_is_assigned_and_both_barcodes_round_trip(pipeline):
    """Cell barcode and UMI come back exactly as 10x wrote them, for all 18,036 reads.

    R1 is 26 nt of barcode and the pattern consumes all of it, so its payload is empty -- which is
    a state the reader has to carry rather than treat as a short read. The cDNA is on R2 and comes
    through untouched.
    """
    checkout, truth = pipeline["checkout"], pipeline["truth"]
    assert checkout["paired"]
    assert checkout["total"] == len(truth)
    assert checkout["assigned"] == checkout["total"]
    assert checkout["unmatched"] == 0
    assert checkout["short_payload"] == 0
    assert pipeline["sample"]["mean_payload_length"] == 0.0, "R1 is barcode and nothing else"

    seen = 0
    for name, tags, sequence, _ in records(pipeline["dir"] / "co" / "PBMC_R2.fq.gz"):
        cell, umi = truth[name]
        assert tags["CB"] == cell
        assert tags["RX"] == umi
        assert len(sequence) > 26, "R2 is 90 nt of cDNA -- it must not have been trimmed"
        seen += 1
    assert seen == len(truth)


def test_the_cell_key_is_the_whole_cell_barcode(pipeline):
    """Cells counted by migec must equal the distinct cell barcodes in the file.

    The cell key is the top `2 * cell_length` bits of the packed barcode. Taking "everything above
    the UMI" instead coincides with that only when cell + UMI fill all 32 bases -- here they fill
    26 -- and getting it wrong shatters each cell into fragments that still look like cells, so
    only an independent count catches it. Folding N to A is what `pack_barcode` does: 1,020
    distinct `CB` strings, 7 of them carrying an N, are 1,014 cells.
    """
    cells, barcodes = set(), set()
    for _, tags, _, _ in records(READS):
        cell, umi = tags["CB"].replace("N", "A"), tags["RX"].replace("N", "A")
        cells.add(cell)
        barcodes.add((cell, umi))

    refine = pipeline["refine"]
    assert refine["cells_observed"] == len(cells)
    assert refine["barcodes"] == len(barcodes)
    assert refine["barcodes"] > refine["cells_observed"], "a cell holds several molecules"


def test_cell_calling_counts_molecules_and_not_reads(pipeline):
    """Molecules per called cell, against the reads that made them.

    The library averages 6.9 reads per barcode, so a cell caller that counted reads would report
    ~17.8 per cell where the molecule count is ~2.6. That gap is the whole assertion: one
    over-amplified molecule must not lift an empty droplet up the rank curve, which is the
    artefact the knee plot exists to show.
    """
    refine, sample = pipeline["refine"], pipeline["sample"]
    assert 0 < refine["cells_called"] <= refine["cells_observed"]
    assert refine["molecules_in_called"] <= refine["molecules"]

    per_cell = refine["molecules_in_called"] / refine["cells_called"]
    reads_per_cell = sample["reads"] / refine["cells_observed"]
    assert 1.0 <= per_cell < reads_per_cell / 2
    assert per_cell < sample["mean_reads_per_umi"]

    ranks = (pipeline["dir"] / "ref" / "PBMC.cells.tsv").read_text().splitlines()
    assert len(ranks) - 1 == refine["cells_observed"]
    called = sum(1 for row in ranks[1:] if row.split("\t")[2] == "1")
    assert called == refine["cells_called"]


def test_the_t_cells_read_as_diverse(pipeline):
    """The contrast the amplicon test is the other half of.

    Payload agreement between two barcodes is worth `log(1 / clonality)`, so on a clonal library
    it is worth nothing and refine must measure the clonality rather than assume it. A thousand T
    cells sequenced over their rearranged receptors are the diverse end: measured 5.7e-4 here
    against 0.80 for the HIV amplicon, from the same field.
    """
    refine = pipeline["refine"]
    assert refine["payload_clonality"] < 0.01
    assert refine["molecules"] == refine["barcodes"] - refine["merged"]
    assert refine["molecules"] <= refine["barcodes"]
    # Payload agreement is what makes a merge possible at 1-3 reads per barcode, where the count
    # ratio is no evidence at all. On a library this diverse it has to fire at least once.
    assert refine["merged_by_payload"] > 0
