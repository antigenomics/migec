"""Every stage, end to end, at several thread counts: the bytes must not move.

The C++ suite asserts this on the stages directly (tests/cpp/test_parallel_stages.cpp) and under
the thread sanitizer. This is the same property at the CLI level, over a simulated library that
goes through all three stages in order -- because the thing that actually breaks is not a race
inside one stage but a stage whose output shifts and a downstream stage that silently accepts it.

Never: `-t` may change the wall clock and nothing else. A pipeline whose results depend on how many
cores the machine had cannot be compared between runs, and no summary reports it.
"""

from __future__ import annotations

import hashlib

import pytest

# Never: BEFORE the `migec` imports below, and not `pytestmark` alone. A module-scope
# import of a package whose extension is missing raises at COLLECTION, which pytest
# reports as an error rather than a skip -- so a machine without the built extension
# fails the suite instead of saying it cannot run it.
pytest.importorskip("migec._core", reason="the C++ extension is not built: run `bash setup.sh`")

from migec.assemble import run as assemble_run
from migec.checkout import run as checkout_run
from migec.refine import run as refine_run

from ._sim import SimConfig, simulate

from tests.conftest import requires_core

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"
THREADS = [1, 2, 3, 8]


def digest(path) -> str:
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    """One simulated library with barcode errors, uneven depth, and enough molecules to bucket."""
    d = tmp_path_factory.mktemp("threads")
    cfg = SimConfig(
        adapter=ADAPTER, n_molecules=3000, n_clones=8, coverage=4.0, coverage_cv=0.9,
        seq_error=3e-3, umi_error=1e-3,
    )
    sim = simulate(cfg, d / "sim")
    (d / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    return d, sim


def test_checkout_is_byte_identical_at_every_thread_count(library):
    d, sim = library
    seen = {}
    for t in THREADS:
        summary = checkout_run(sim["reads"], d / "bc.txt", d / f"co{t}", threads=t)
        seen[t] = digest(d / f"co{t}" / "S1.fq.gz")
        assert summary["threads"] == t
    assert len(set(seen.values())) == 1, seen


def test_refine_is_byte_identical_at_every_thread_count(library):
    d, _ = library
    seen, decisions = {}, {}
    for t in THREADS:
        summary = refine_run(d / "co1" / "S1.fq.gz", d / f"rf{t}", threads=t)
        seen[t] = digest(d / f"rf{t}" / "S1.fq.gz")
        # The decisions, not only the bytes: a merge that moved would change the molecule count.
        decisions[t] = (summary["merged"], summary["merged_reads"], summary["molecules"],
                        round(summary["estimated_error"], 12))
    assert len(set(seen.values())) == 1, seen
    assert len(set(decisions.values())) == 1, decisions


def test_assemble_is_byte_identical_at_every_thread_count(library):
    d, _ = library
    fastq, table, stats = {}, {}, {}
    for t in THREADS:
        summary = assemble_run(d / "rf1" / "S1.fq.gz", d / f"as{t}", threads=t)
        fastq[t] = digest(d / f"as{t}" / "S1.consensus.fq.gz")
        table[t] = digest(d / f"as{t}" / "S1.mig.tsv")
        stats[t] = (summary["groups"], summary["molecules"], summary["groups_split"],
                    summary["buckets"])
    assert len(set(fastq.values())) == 1, fastq
    assert len(set(table.values())) == 1, table
    # The bucket count is deliberately NOT a function of -t: if it were, it would choose the gzip
    # member boundaries and the files would differ byte-wise while holding identical records.
    assert len(set(stats.values())) == 1, stats


def test_the_consensus_is_sorted_by_barcode_however_many_buckets(library):
    """Bucket order is key order, so concatenating the buckets in order yields a sorted file. That
    is what lets the merge be a concatenation rather than a k-way merge."""
    import gzip

    d, _ = library
    for t in (1, 8):
        umis = []
        with gzip.open(d / f"as{t}" / "S1.consensus.fq.gz", "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    umis.append(next(f[5:] for f in line.split() if f.startswith("RX:Z:")))
        assert umis == sorted(umis)
        assert len(umis) > 100


def test_a_limited_run_says_so_rather_than_looking_like_a_small_library(library):
    d, _ = library
    from migec.assemble import format_report

    summary = assemble_run(d / "rf1" / "S1.fq.gz", d / "limited", limit_reads=500)
    assert summary["reads"] == 500
    assert summary["limited"]
    assert "limited" in format_report(summary)

    by_umi = assemble_run(d / "rf1" / "S1.fq.gz", d / "limited_umi", limit_umis=100)
    assert by_umi["limited"]
    assert by_umi["groups"] <= 100
    # The barcodes brought their reads with them, which is what makes this a usable smoke test.
    assert by_umi["reads"] > by_umi["groups"]


def test_the_whole_pipeline_agrees_with_itself_at_sixteen_threads(library):
    """The end-to-end check: run all three stages at -t 16 and compare against the -t 1 chain."""
    d, sim = library
    checkout_run(sim["reads"], d / "bc.txt", d / "e2e_co", threads=16)
    refine_run(d / "e2e_co" / "S1.fq.gz", d / "e2e_rf", threads=16)
    assemble_run(d / "e2e_rf" / "S1.fq.gz", d / "e2e_as", threads=16)

    assert digest(d / "e2e_co" / "S1.fq.gz") == digest(d / "co1" / "S1.fq.gz")
    assert digest(d / "e2e_rf" / "S1.fq.gz") == digest(d / "rf1" / "S1.fq.gz")
    assert digest(d / "e2e_as" / "S1.consensus.fq.gz") == digest(d / "as1" / "S1.consensus.fq.gz")
