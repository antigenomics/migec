"""The pure functions behind the ctDNA titration and the SRA bootstrap.

Both scripts talk to the network for their real work, so what is testable here is the part that
decides what the network answers MEAN: the three design grammars, the read-structure verdict, the
amplicon tally and the detectability arithmetic. Each of these got something wrong once and the
wrong answer looked plausible, which is why they are pinned.
"""

from __future__ import annotations

import gzip
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent.parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ct = _load("ctdna_titration")
sf = _load("sra_fetch")


# --- the design grammars ------------------------------------------------------------------------


@pytest.mark.parametrize(
    "label, ng, arm, vaf, rep",
    [
        # PRJNA788522 puts the design in sample_alias.
        ("20ng_0.125_10x_rep_1", 20, "0.125", 0.00125, 1),
        ("5ng_WT_3.3x_rep_2", 5, "WT", 0.0, 2),
        ("80ng_1_10x_rep_3", 80, "1", 0.01, 3),
    ],
)
def test_the_titration_alias_parses(label, ng, arm, vaf, rep):
    d = ct.parse_label(label)
    assert (d["input_ng"], d["arm"], d["vaf"], d["replicate"]) == (ng, arm, vaf, rep)


def test_an_undiluted_arm_has_no_frequency_rather_than_a_guessed_one():
    """`cell_line` is undiluted reference material and the paper does not state its VAF.

    Never: an unknown frequency must stay None. Filling it in with 1.0, or with the highest
    dilution, would put a fabricated number straight into the detectability column.
    """
    d = ct.parse_label("20ng_cell_line_30x_rep_1")
    assert d["vaf"] is None
    assert d["input_ng"] == 20


@pytest.mark.parametrize(
    "label, enzyme, vaf",
    [
        # PRJNA507366 puts it in library_name instead -- every sample_alias there is the same
        # string, so a parser that only reads the alias calls a designed study undesigned.
        ("Phusion1|Wildtype", "Phusion", 0.0),
        ("Accuprime_hifi2|Wildtype", "Accuprime_hifi", 0.0),
        ("3plx Platinum superfi 70ng 0.0625% (3)|0.0625% VAF", "Platinum superfi", 0.000625),
        ("3plx Platinum superfi 50 ng 0.125%|0.125% VAF", "Platinum superfi", 0.00125),
    ],
)
def test_the_polymerase_label_parses(label, enzyme, vaf):
    d = ct.parse_label(label)
    assert d["enzyme"] == enzyme
    assert d["vaf"] == pytest.approx(vaf)


@pytest.mark.parametrize(
    "label, enzyme, ng, rep",
    [
        # A two-enzyme comparison at fixed input. The enzyme that matters is the one AFTER the
        # colon: the prefix names the prep, the suffix names the polymerase actually used.
        ("Seracare WT (1) 80ng PlatSuperfi:PlatRegular|Wildtype", "PlatRegular", 80, 1),
        ("Seracare WT (2) 80ng PlatSuperfi:PlatSuperfi|Wildtype", "PlatSuperfi", 80, 2),
    ],
)
def test_the_paired_enzyme_label_parses(label, enzyme, ng, rep):
    d = ct.parse_label(label)
    assert (d["enzyme"], d["input_ng"], d["replicate"], d["vaf"]) == (enzyme, ng, rep, 0.0)


def test_an_unparseable_label_survives_as_itself():
    """A label nobody wrote a grammar for must still run, carrying its raw text and no design.

    Never: the fallback must not invent a frequency. A run whose label says nothing about VAF
    gets `None` and is reported without a detectability column, rather than being assigned 0 --
    which would silently enrol it in the true-negative arm.
    """
    d = ct.parse_label("some sample nobody documented")
    assert d["vaf"] is None
    assert d["enzyme"] == ""
    assert d["arm"] == "some sample nobody documented"


# --- the detectability arithmetic ---------------------------------------------------------------


def test_enough_molecules_makes_detection_certain_and_too_few_makes_it_a_coin_flip():
    """The number the whole page turns on: support is a Poisson draw, not a guarantee."""
    plenty = ct.detectable(10_000, 0.01, min_support=3)      # expect 100
    assert plenty["variant_molecules"] == 100
    assert plenty["p_enough"] == pytest.approx(1.0, abs=1e-6)

    marginal = ct.detectable(2_000, 0.00125, min_support=3)  # expect 2.5
    assert marginal["variant_molecules"] == 2.5
    assert 0.4 < marginal["p_enough"] < 0.6


def test_a_zero_frequency_arm_expects_zero_variant_molecules():
    """The WT arm is the true negative; its expectation must be exactly zero, not a small number."""
    assert ct.detectable(50_000, 0.0, min_support=3)["variant_molecules"] == 0.0


def test_an_unknown_frequency_reports_nothing_rather_than_zero():
    assert ct.detectable(50_000, None, min_support=3) == {"variant_molecules": "", "p_enough": ""}


def test_genome_equivalents_uses_the_haploid_mass():
    # 20 ng / 3.3 pg is ~6,060 copies -- the ceiling on molecules before any loss.
    assert ct.genome_equivalents(20) == pytest.approx(6060.6, rel=1e-3)


# --- the amplicon tally -------------------------------------------------------------------------


def _consensus(tmp_path: Path, seqs: list[str]) -> Path:
    p = tmp_path / "c.fq.gz"
    with gzip.open(p, "wt") as fh:
        for i, s in enumerate(seqs):
            fh.write(f"@m{i}\n{s}\n+\n{'I' * len(s)}\n")
    return p


def test_the_share_floor_separates_amplicons_from_payload_error(tmp_path):
    """Four real amplicons at 24% each, plus a long tail of one-off error prefixes.

    Never: the floor cannot sit on the error slope. At 1% the singletons here get counted as
    amplicons, which inflates the divisor and deflates molecules-per-site on exactly the runs
    where the evidence is thinnest -- shallow ones, whose consensus carries the most payload error.
    """
    seqs = []
    for base in ("AAAA", "CCCC", "GGGG", "TTTT"):
        seqs += [base + "ACGTACGTACGTACGTACGTA"] * 60
    seqs += [f"ACGT{i:017d}" for i in range(40)]      # 40 distinct singletons, 14% together

    fq = _consensus(tmp_path, seqs)
    n_gap, share = ct.count_amplicons(fq, prefix_len=25, min_share=0.05)
    assert n_gap == 4
    assert share == pytest.approx(240 / 280, rel=1e-3)

    n_slope, _ = ct.count_amplicons(fq, prefix_len=25, min_share=0.001)
    assert n_slope == 44, "a floor on the error slope counts singletons as amplicons"


def test_an_empty_consensus_is_zero_amplicons_not_a_crash(tmp_path):
    assert ct.count_amplicons(_consensus(tmp_path, []), prefix_len=25, min_share=0.05) == (0, 0.0)


# --- the read-structure verdict -----------------------------------------------------------------


@pytest.mark.parametrize(
    "row, expect",
    [
        # Three reads per spot: the UMI may be on its own index read.
        ({"nreads": "3", "read_lengths": "8,101,101"}, "3+ reads per spot"),
        # Mates of visibly different length: one may still carry an inline barcode. This is
        # SRR15081472, which Maruzani used.
        ({"nreads": "2", "read_lengths": "101,110"}, "mates differ in length"),
        # Two equal mates: no separate UMI read. SRR10296599.
        ({"nreads": "2", "read_lengths": "96,96"}, "two equal mates"),
        ({"nreads": "1", "read_lengths": "151"}, "single read"),
    ],
)
def test_the_verdict_never_claims_more_than_metadata_can_show(row, expect):
    """Every verdict ends in "peek to confirm" on purpose.

    Never: read structure rules a run OUT, it cannot rule one IN -- an inline UMI is invisible in
    the metadata, which is exactly how SiMSen-Seq's 12 nt barcode went unnoticed in a literature
    that concluded no public ctDNA data had usable UMIs.
    """
    v = sf.verdict(row)
    assert expect in v
    assert v.endswith("peek to confirm")


# --- assay recommendations ------------------------------------------------------------------------


def test_every_layout_preset_carries_a_variant_calling_recommendation():
    """A preset says where the barcode is; an assay also needs to say what a consensus is worth.
    Never: leaving one out means `migec sheet --presets` silently recommends nothing for it."""
    from migec.sheet import ASSAYS, PRESETS

    assert set(PRESETS) <= set(ASSAYS), f"no recommendation for {set(PRESETS) - set(ASSAYS)}"
    for name, a in ASSAYS.items():
        assert a["min_reads_variant"] >= 3, (
            f"{name}: --min-reads 1 means a singleton consensus IS one raw read, with no error "
            f"correction -- measured to carry every 2-colour dark-G false positive")
        assert a["rt_error"], f"{name}: needs a pre-amplification floor"
        assert a["note"], f"{name}: a recommendation without a reason is folklore"
