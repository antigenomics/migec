"""subsample: a smaller library that is still a library.

The failure this guards against is silent. Sampling reads, or taking the first N barcodes, both
produce a file that looks fine and has the MIG size distribution destroyed -- so every fixture
built from it tests a library nobody has.
"""

from __future__ import annotations

import collections
import gzip

import pytest

# Never: BEFORE the `migec` imports below, and not `pytestmark` alone. A module-scope
# import of a package whose extension is missing raises at COLLECTION, which pytest
# reports as an error rather than a skip -- so a machine without the built extension
# fails the suite instead of saying it cannot run it.
pytest.importorskip("migec._core", reason="the C++ extension is not built: run `bash setup.sh`")

from migec.checkout import run as checkout_run
from migec.subsample import format_report, run

from ._sim import SimConfig, simulate

from tests.conftest import requires_core

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def umis_of(path):
    counts = collections.Counter()
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                counts[next(f[5:] for f in line.split() if f.startswith("RX:Z:"))] += 1
    return counts


@pytest.fixture
def library(tmp_path):
    sim = simulate(
        SimConfig(adapter=ADAPTER, n_molecules=4_000, n_clones=20, coverage=6.0,
                  coverage_cv=0.6, umi_len=12, seed=11),
        tmp_path / "sim",
    )
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co")
    return tmp_path / "co" / "S1.fq.gz"


def test_every_read_of_a_kept_barcode_is_kept(library, tmp_path):
    """The definition. A barcode is in or out; its reads never split."""
    full = umis_of(library)
    s = run(library, tmp_path / "s.fq.gz", keep_percent=10.0)
    kept = umis_of(tmp_path / "s.fq.gz")
    assert set(kept) <= set(full)
    for umi, n in kept.items():
        assert n == full[umi], f"{umi} kept {n} of its {full[umi]} reads"
    assert s["barcodes"] == len(kept)
    assert s["reads_kept"] == sum(kept.values())


def test_the_mig_size_distribution_survives(library, tmp_path):
    """What sampling reads destroys, and the reason this command exists."""
    full = umis_of(library)
    run(library, tmp_path / "s.fq.gz", keep_percent=20.0)
    kept = umis_of(tmp_path / "s.fq.gz")

    full_mean = sum(full.values()) / len(full)
    kept_mean = sum(kept.values()) / len(kept)
    assert kept_mean == pytest.approx(full_mean, rel=0.15)

    # ...and the shape, not just the mean: the share of molecules seen more than once.
    full_multi = sum(1 for v in full.values() if v > 1) / len(full)
    kept_multi = sum(1 for v in kept.values() if v > 1) / len(kept)
    assert kept_multi == pytest.approx(full_multi, abs=0.05)


def test_sampling_reads_would_have_destroyed_it(library, tmp_path):
    """The comparison that justifies the rule, run rather than asserted."""
    import random

    full = umis_of(library)
    run(library, tmp_path / "s.fq.gz", keep_percent=20.0)
    by_barcode = umis_of(tmp_path / "s.fq.gz")

    rng = random.Random(0)
    by_read = collections.Counter()
    for umi, n in full.items():
        drawn = sum(1 for _ in range(n) if rng.random() < 0.2)
        if drawn:
            by_read[umi] = drawn

    full_mean = sum(full.values()) / len(full)
    assert sum(by_barcode.values()) / len(by_barcode) == pytest.approx(full_mean, rel=0.15)
    # Sampling reads at the same rate collapses the mean towards one read per molecule.
    assert sum(by_read.values()) / len(by_read) < 0.5 * full_mean


def test_selection_is_deterministic(library, tmp_path):
    run(library, tmp_path / "a.fq.gz", keep_percent=10.0)
    run(library, tmp_path / "b.fq.gz", keep_percent=10.0)
    assert (tmp_path / "a.fq.gz").read_bytes() == (tmp_path / "b.fq.gz").read_bytes()
    # ...and nested: a smaller fraction is a subset of a larger one, since the test is on the
    # same hash. That is what makes a fixture shrinkable without re-deriving it.
    run(library, tmp_path / "c.fq.gz", keep_percent=5.0)
    assert set(umis_of(tmp_path / "c.fq.gz")) <= set(umis_of(tmp_path / "a.fq.gz"))


def test_a_cell_is_kept_whole(tmp_path):
    """Sampling molecules independently gives thousands of cells holding one molecule each, which
    is the read-sampling mistake wearing a different hat."""
    reads = tmp_path / "cells.fq.gz"
    import random

    rng = random.Random(3)
    cells = ["".join(rng.choice("ACGT") for _ in range(16)) for _ in range(200)]
    with gzip.open(reads, "wt") as fh:
        i = 0
        for cell in cells:
            for _ in range(30):
                umi = "".join(rng.choice("ACGT") for _ in range(12))
                for _ in range(3):
                    fh.write(
                        f"@r{i} RX:Z:{umi}\tQX:Z:{'I' * 12}\tCB:Z:{cell}\tBC:Z:S1\n"
                        f"ACGT\n+\nIIII\n"
                    )
                    i += 1
    run(reads, tmp_path / "s.fq.gz", keep_percent=20.0)

    per_cell = collections.Counter()
    with gzip.open(tmp_path / "s.fq.gz", "rt") as fh:
        for j, line in enumerate(fh):
            if j % 4 == 0:
                cb = next(f[5:] for f in line.split() if f.startswith("CB:Z:"))
                per_cell[cb] += 1
    assert per_cell, "nothing was kept"
    # Every cell that survived kept all 30 of its molecules, 3 reads each.
    for cell, n in per_cell.items():
        assert n == 90, f"cell {cell} kept {n} of its 90 reads"


def test_the_keep_fraction_is_validated(library, tmp_path):
    with pytest.raises(ValueError, match="ten-thousandths"):
        run(library, tmp_path / "s.fq.gz", keep_percent=0.0)
    with pytest.raises(ValueError, match="ten-thousandths"):
        run(library, tmp_path / "s.fq.gz", keep_percent=200.0)


def test_the_report_says_what_was_kept(library, tmp_path):
    s = run(library, tmp_path / "s.fq.gz", keep_percent=10.0)
    report = format_report(s)
    assert "reads per barcode" in report
    assert f"{s['barcodes']:,}" in report
