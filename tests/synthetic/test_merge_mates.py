"""`--mate2` / `--merge-mates`: the two ends of one molecule, placed rather than merged.

An overlapping pair comes back as one consensus spanning the insert. A pair whose mates do NOT
overlap comes back as two contigs, because bridging them would assert bases that no read covers --
the same rule `--contig` follows, applied to the pair.
"""

from __future__ import annotations

import gzip

import pytest

pytest.importorskip("migec._core", reason="the C++ extension is not built: run `bash setup.sh`")

from migec.assemble import format_report, run
from migec.checkout import run as checkout_run

from ._sim import SimConfig, simulate

from tests.conftest import requires_core

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"
SEQ_LEN = 120


def build(tmp_path, read_len, **kwargs):
    """Simulate a paired library, check it out paired, and return (sim, checkout dir)."""
    cfg = SimConfig(
        adapter=ADAPTER, seq_len=SEQ_LEN, paired=True, read_len=read_len,
        n_clones=8, n_molecules=300, coverage=6.0, coverage_cv=0.3, **kwargs,
    )
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co", reads2=sim["reads2"])
    return sim, tmp_path / "co"


def consensuses(path):
    """name -> sequence, from a consensus FASTQ."""
    out = {}
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0:
                name = line[1:].split()[0]
            elif i % 4 == 1:
                out[name] = line.rstrip("\n")
    return out


def truth(sim):
    """umi -> full insert."""
    out = {}
    with open(sim["truth_consensus"]) as fh:
        umi = None
        for line in fh:
            if line.startswith(">"):
                umi = next(f[4:] for f in line.split() if f.startswith("umi="))
            else:
                out[umi] = line.rstrip("\n")
    return out


def test_overlapping_mates_come_back_as_one_consensus_of_the_whole_insert(tmp_path):
    # 80 + 80 over a 120 nt insert: the mates overlap by 40.
    sim, co = build(tmp_path, read_len=80)
    summary = run(co / "S1_R1.fq.gz", tmp_path / "asm", mate2=co / "S1_R2.fq.gz")
    got = consensuses(tmp_path / "asm" / "S1.consensus.fq.gz")
    want = truth(sim)

    # One contig per molecule: a merged pair has no `.c1` suffix in its name.
    assert summary["groups"] > 200
    assert summary["groups_fragmented"] == 0
    assert all(".c" not in name for name in got)
    # ...and it spans the insert, not one mate.
    full = [s for s in got.values() if len(s) == SEQ_LEN]
    assert len(full) == len(got)
    exact = sum(1 for name, seq in got.items() if seq == want.get(name.split(".")[-1]))
    assert exact / len(got) > 0.95, f"only {exact}/{len(got)} consensuses match the true insert"
    assert "mates merged" in format_report(summary)


def test_mates_that_do_not_overlap_stay_two_contigs(tmp_path):
    # 45 + 45 over a 120 nt insert leaves 30 bases that neither mate covers.
    _, co = build(tmp_path, read_len=45)
    summary = run(co / "S1_R1.fq.gz", tmp_path / "asm", mate2=co / "S1_R2.fq.gz")
    got = consensuses(tmp_path / "asm" / "S1.consensus.fq.gz")

    assert summary["groups_fragmented"] == summary["groups"]
    assert all(".c" in name for name in got)
    # Never bridged: no emitted contig is longer than the mate it came from.
    assert max(len(s) for s in got.values()) <= 45


def test_the_mig_route_merges_the_pair_already_in_the_record(tmp_path):
    """`checkout --mig` stores both mates, so --merge-mates alone is enough there."""
    cfg = SimConfig(
        adapter=ADAPTER, seq_len=SEQ_LEN, paired=True, read_len=80,
        n_clones=8, n_molecules=300, coverage=6.0, coverage_cv=0.3,
    )
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    # The same reads through both of checkout's output routes: buckets, and the FASTQ pair.
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co_mig", reads2=sim["reads2"],
                 mig=True)
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co", reads2=sim["reads2"])
    buckets = sorted((tmp_path / "co_mig").glob("S1.*.mig"))
    assert buckets, "checkout --mig wrote no buckets"

    fastq = run(
        (tmp_path / "co") / "S1_R1.fq.gz", tmp_path / "asm_fq",
        mate2=(tmp_path / "co") / "S1_R2.fq.gz",
    )
    mig = run(buckets[0], tmp_path / "asm_mig", merge_mates=True)
    # The two routes are the same partition of the same reads, so they are the same molecules.
    assert mig["groups"] == fastq["groups"]
    assert mig["molecules"] == fastq["molecules"]
    assert (consensuses(tmp_path / "asm_mig" / "S1.consensus.fq.gz")
            == consensuses(tmp_path / "asm_fq" / "S1.consensus.fq.gz"))


def test_merging_is_byte_identical_at_any_thread_count(tmp_path):
    """`-t` changes the wall clock and nothing else, on this path as on every other."""
    _, co = build(tmp_path, read_len=80)
    one = run(co / "S1_R1.fq.gz", tmp_path / "t1", mate2=co / "S1_R2.fq.gz", threads=1)
    eight = run(co / "S1_R1.fq.gz", tmp_path / "t8", mate2=co / "S1_R2.fq.gz", threads=8)
    assert one["molecules"] == eight["molecules"]
    assert ((tmp_path / "t1" / "S1.consensus.fq.gz").read_bytes()
            == (tmp_path / "t8" / "S1.consensus.fq.gz").read_bytes())
    assert ((tmp_path / "t1" / "S1.mig.tsv").read_bytes()
            == (tmp_path / "t8" / "S1.mig.tsv").read_bytes())


def test_a_single_end_input_is_refused_rather_than_half_assembled(tmp_path):
    cfg = SimConfig(adapter=ADAPTER, seq_len=SEQ_LEN, n_clones=4, n_molecules=50, coverage=4.0)
    sim = simulate(cfg, tmp_path / "sim")
    (tmp_path / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
    checkout_run(sim["reads"], tmp_path / "bc.txt", tmp_path / "co", mig=True)
    bucket = sorted((tmp_path / "co").glob("S1.*.mig"))[0]
    with pytest.raises(RuntimeError, match="single-end"):
        run(bucket, tmp_path / "asm", merge_mates=True)


def test_mates_of_different_lengths_are_refused(tmp_path):
    """Matched by POSITION, so a file that runs out early is two runs rather than a pair."""
    sim, co = build(tmp_path, read_len=80)
    short = tmp_path / "short_R2.fq.gz"
    with gzip.open(co / "S1_R2.fq.gz", "rt") as src, gzip.open(short, "wt") as dst:
        for i, line in enumerate(src):
            if i >= 400:
                break
            dst.write(line)
    with pytest.raises(RuntimeError, match="not two mates of one run"):
        run(co / "S1_R1.fq.gz", tmp_path / "asm", mate2=short)
