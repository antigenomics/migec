"""Telling a real variant from a systematic artifact with no panel of normals.

The discriminator is how a variant's allele frequency moves as `--min-reads` rises. A real variant
sits in a fixed fraction of molecules whatever built them; a context artifact rides on molecules
made from few reads, because a singleton consensus IS one raw read and carries the raw error rate.

The numbers in the fixture are measured, from `SRR17220923` (certified 1%, 20 ng, 10x).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
spec = importlib.util.spec_from_file_location(
    "blind_artifact_filter", ROOT / "scripts" / "blind_artifact_filter.py")
bf = importlib.util.module_from_spec(spec)
sys.modules["blind_artifact_filter"] = bf
spec.loader.exec_module(bf)

VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def write_vcf(path: Path, calls: list[tuple[str, str, str, str, float]]) -> Path:
    body = "".join(f"{c}\t{p}\t.\t{r}\t{a}\t100\tPASS\tDP=8000;AF={af}\n" for c, p, r, a, af in calls)
    path.write_text(VCF_HEAD + body)
    return path


def test_a_stable_vaf_is_real_and_a_vanishing_one_is_artifact(tmp_path):
    """The measured case: H1047R holds its frequency, the dark-G calls disappear."""
    mr1 = write_vcf(tmp_path / "1.vcf", [
        ("3", "179234297", "A", "G", 0.011839),   # PIK3CA H1047R, real
        ("17", "7674220", "C", "T", 0.013857),    # real
        ("17", "7673768", "T", "G", 0.004094),    # dark-G artifact
        ("4", "54733163", "A", "G", 0.006812),    # dark-G artifact
    ])
    mr3 = write_vcf(tmp_path / "3.vcf", [
        ("3", "179234297", "A", "G", 0.011819),
        ("17", "7674220", "C", "T", 0.013744),
    ])
    mr5 = write_vcf(tmp_path / "5.vcf", [
        ("3", "179234297", "A", "G", 0.012533),
        ("17", "7674220", "C", "T", 0.013802),
    ])
    t = {1: bf.read_vcf(mr1), 3: bf.read_vcf(mr3), 5: bf.read_vcf(mr5)}

    def verdict(key):
        return bf.classify({mr: t[mr].get(key) for mr in (1, 3, 5)}, drop=0.5, min_seen=2)[0]

    assert verdict(("3", "179234297", "A", "G")) == "real"
    assert verdict(("17", "7674220", "C", "T")) == "real"
    assert verdict(("17", "7673768", "T", "G")) == "artifact"
    assert verdict(("4", "54733163", "A", "G")) == "artifact"


def test_a_call_absent_at_the_lowest_threshold_is_inconclusive_not_real():
    """Never: a variant appearing only at a HIGHER threshold is the caller gaining power, not a
    variant surviving a filter. Scoring it as real would invert the test's logic."""
    vafs = {1: None, 3: 0.0018, 5: 0.0020}
    assert bf.classify(vafs, drop=0.5, min_seen=2)[0] == "inconclusive"


def test_a_halved_vaf_sits_on_the_boundary():
    """`--drop` is the knob; at 0.5 an exactly halved frequency is kept, less is not."""
    assert bf.classify({1: 0.010, 3: 0.005}, drop=0.5, min_seen=2)[0] == "real"
    assert bf.classify({1: 0.010, 3: 0.004}, drop=0.5, min_seen=2)[0] == "artifact"


def test_seen_at_too_few_thresholds_is_artifact():
    assert bf.classify({1: 0.01, 3: None, 5: None}, drop=0.5, min_seen=2)[0] == "artifact"


def test_read_vcf_takes_af_from_info_not_from_a_substring(tmp_path):
    """Never: `AF=` also matches inside `MQ_AF=` or `SOMEAF=`. Anchor it."""
    p = tmp_path / "x.vcf"
    p.write_text(VCF_HEAD + "1\t100\t.\tA\tG\t100\tPASS\tOTHERAF=0.9;AF=0.0123;DP=10\n")
    assert bf.read_vcf(p)[("1", "100", "A", "G")] == pytest.approx(0.0123)
