"""The per-target analysis, whose job is to replace an average with a count.

`ctdna_titration.py` divides a library total by an amplicon count inferred from consensus
prefixes. This one takes molecules actually aligned to each target. The tests pin the parts where
getting it wrong would quietly reintroduce the average, or mix two studies' panels together.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent
spec = importlib.util.spec_from_file_location("ctdna_persite", ROOT / "scripts" / "ctdna_persite.py")
cp = importlib.util.module_from_spec(spec)
sys.modules["ctdna_persite"] = cp
spec.loader.exec_module(cp)


@pytest.mark.parametrize(
    "alias, ng, arm, vaf, depth, diluted",
    [
        ("20ng_WT_10x_rep_3", 20, "WT", 0.0, 10.0, True),
        ("5ng_0.125_3.3x_rep_1", 5, "0.125", 0.00125, 3.3, True),
        ("20ng_0.25_10x_rep_2", 20, "0.25", 0.0025, 10.0, True),
        ("80ng_1_30x_rep_1", 80, "1", 0.01, 30.0, True),
    ],
)
def test_the_certified_arms_parse_with_their_frequencies(alias, ng, arm, vaf, depth, diluted):
    d = cp.parse_alias(alias)
    assert (d["input_ng"], d["arm"], d["vaf"], d["depth"], d["diluted"]) == \
        (ng, arm, vaf, depth, diluted)


def test_the_undiluted_arm_has_no_certified_frequency():
    """Never: `cell_line` is raw material. It is neither a 0% control nor a certified frequency,
    and coding it as either puts a fabricated number into the comparison against truth."""
    d = cp.parse_alias("20ng_cell_line_30x_rep_2")
    assert d["arm"] == "cell_line"
    assert d["vaf"] is None
    assert d["diluted"] is False
    assert d["input_ng"] == 20


@pytest.mark.parametrize(
    "alias",
    [
        "Phusion1|Wildtype",                                  # PRJNA507366 polymerase arm
        "3plx Platinum superfi 70ng 0.0625% (3)|0.0625% VAF",  # its dilution arm
        "SeraCare_Reference_Material",                         # its uninformative sample_alias
    ],
)
def test_a_different_study_is_excluded_rather_than_misparsed(alias):
    """Never: the two studies run DIFFERENT panels -- 5 targets against a 3-plex. Letting the
    other study's runs through would count its molecules against this panel's intervals and
    silently average two assays together."""
    assert cp.parse_alias(alias)["input_ng"] is None


def test_the_certified_frequencies_are_the_vendors():
    assert cp.ARMS == {"WT": 0.0, "0.125": 0.00125, "0.25": 0.0025, "1": 0.01}


def test_poisson_tail_is_the_closed_form():
    import math
    for lam in (0.5, 2.0, 8.0):
        assert cp.poisson_at_least(lam, 1) == pytest.approx(1 - math.exp(-lam))
    assert cp.poisson_at_least(0.0, 3) == 0.0


def test_read_tsv_keys_on_the_header(tmp_path):
    """The columns come off a cluster job whose field order is not guaranteed; reading them
    positionally is how a molecule count ends up in the gene column."""
    p = tmp_path / "m.tsv"
    p.write_text("run\tchrom\tstart\tend\tgene\tmolecules\nS1\t17\t100\t200\tTP53\t42\n")
    rows = cp.read_tsv(p)
    assert rows == [{"run": "S1", "chrom": "17", "start": "100", "end": "200",
                     "gene": "TP53", "molecules": "42"}]
