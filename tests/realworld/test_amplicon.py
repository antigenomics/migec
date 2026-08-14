"""SRR1763769, an HIV-1 Primer ID amplicon, through checkout -> refine -> assemble.

The fixture is `ci/SRR1763769_umi0.5pct.fq.gz`: all the reads of 0.5% of the barcodes of a 2.12 M
read library, 9,824 reads over 623 barcodes at 15.8 reads each (`SOURCES.md`, and the dataset's
own README). Whole barcodes, never a read sample -- which is the only reason a consensus test on
it means anything.

Two things this fixture cannot do, stated here so nothing downstream over-reads it:

* **It is post-checkout.** The 9 nt Primer ID and the `CAGTTTAACTTTTGGGCCATCCA` primer were
  trimmed off when the fixture was made, so `checkout` is exercised by putting the primer back
  (`_restore_untrimmed`). The UMI bases and their qualities are the library's own -- they come
  back out of `RX`/`QX` -- and every number asserted below is computed from those. The primer's
  own bases are re-synthesised from the barcode table at a flat Q40, so "assigned" here tests
  placement, trimming and header transfer, NOT that the primer is present in the library. Never
  assert on `checkout.quality_calibration.tsv` from this file: it measures the error rate against
  the pattern's constant bases, and here those bases are ours.
* **It cannot pin the barcode error RATE.** The 0.5% subsample keeps a barcode when
  `hash(barcode)` falls in the kept range, and an error child hashes independently of its parent,
  so ~99.5% of the children the estimator counts were dropped with their parents kept. Measured
  below: 21 distance-1 pairs where 20.0 are expected by chance alone. `estimate_umi_error` finds
  no excess, returns its 1e-4 floor, and the ratio to the Phred prediction is 0.05 -- so the test
  asserts a DIRECTION (finite, positive, at or below what the chemistry predicts) and checks the
  missing excess itself, rather than quoting an order-of-magnitude bound the fixture cannot back.
"""

from __future__ import annotations

import gzip
import itertools
import math
from pathlib import Path

import pytest

from tests.conftest import UMI_DATA, requires_core, requires_umi_data

from ._fixture import records

pytestmark = [requires_core, requires_umi_data]

READS = UMI_DATA / "ci" / "SRR1763769_umi0.5pct.fq.gz"
BARCODES = UMI_DATA / "ci" / "SRR1763769_barcodes.txt"

# Column 2 of the barcode table, upper-cased: `NNNNNNNNNcagtttaacttttgggccatcca`. The 9 N are the
# Primer ID, the rest is the gene-specific primer that places it.
PRIMER = "CAGTTTAACTTTTGGGCCATCCA"


def _restore_untrimmed(source: Path, dest: Path) -> dict[str, tuple[str, str]]:
    """Re-prefix each read with the UMI and the primer, undoing the fixture's own checkout.

    Returns `name -> (umi, payload)`, which is the truth the round-trip is checked against.

    Note: the primer is given a flat Q40 rather than a quality borrowed from the read. Acceptance
    is a quality-weighted likelihood ratio, so a borrowed quality decides the answer: charging the
    23 primer bases the first payload base's Phred refused 5 of the 9,824 reads, every one of them
    a read whose payload starts at Q2. That refusal would have been an artefact of the
    reconstruction and nothing else.
    """
    truth: dict[str, tuple[str, str]] = {}
    with gzip.open(dest, "wt") as out:
        for name, tags, sequence, quality in records(source):
            truth[name] = (tags["RX"], sequence)
            out.write(
                f"@{name}\n{tags['RX']}{PRIMER}{sequence}\n+\n"
                f"{tags['QX']}{'I' * len(PRIMER)}{quality}\n"
            )
    return truth


@pytest.fixture(scope="module")
def pipeline(tmp_path_factory):
    """checkout -> refine -> assemble, run once for the module. Under two seconds end to end."""
    from migec.assemble import run as assemble_run
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run

    work = tmp_path_factory.mktemp("amplicon")
    truth = _restore_untrimmed(READS, work / "untrimmed.fq.gz")
    checkout = checkout_run(work / "untrimmed.fq.gz", BARCODES, work / "co")
    refine = refine_run(work / "co" / "CTRL.fq.gz", work / "ref")
    assemble = assemble_run(work / "ref" / "CTRL.fq.gz", work / "asm")
    return {
        "truth": truth,
        "checkout": checkout,
        "sample": checkout["samples"][0],
        "refine": refine,
        "assemble": assemble,
        "dir": work,
    }


def test_the_fixture_is_whole_barcodes_and_not_a_read_sample(pipeline):
    """Everything below assumes deep MIGs. This is the assertion that says they are still there.

    Subsampling reads instead of barcodes gives ~1 read per molecule, at which point the consensus
    is the read, the split threshold is inert and the quality cap is never reached -- the file
    would still parse and every other test here would still pass while testing nothing.
    """
    sample = pipeline["sample"]
    assert sample["mean_reads_per_umi"] > 10, "the published figure is 15.77 reads per barcode"
    assert sample["reads_in_migs_ge5"] > 0.5, "most reads must sit in a MIG worth assembling"


def test_every_read_is_assigned_and_the_umi_and_payload_round_trip(pipeline):
    """The pattern places on this library, and what it captures and trims is exactly right.

    The round-trip is the load-bearing half: for all 9,824 reads the UMI `checkout` extracts must
    equal the one the fixture carries, and the payload it leaves must equal the fixture's sequence
    byte for byte. An off-by-one in the match offset, in the trim point, or in the quality
    indexing (MIGEC v1 indexed quality from the read start rather than the match offset) moves one
    of the two and nothing else in the pipeline would notice.
    """
    checkout, truth = pipeline["checkout"], pipeline["truth"]
    assert checkout["total"] == len(truth)
    assert checkout["assigned"] == checkout["total"]
    assert checkout["unmatched"] == 0
    assert checkout["ambiguous"] == 0
    assert checkout["bad_umi"] == 0

    seen = 0
    for name, tags, sequence, _ in records(pipeline["dir"] / "co" / "CTRL.fq.gz"):
        umi, payload = truth[name]
        assert tags["RX"] == umi
        assert tags["BC"] == "CTRL"
        assert sequence == payload
        seen += 1
    assert seen == len(truth)


def test_the_umi_is_nine_bases_and_the_effective_length_stays_under_it(pipeline):
    """Recovered length and composition, both from the library's own barcodes.

    `effective_length = log4(prod m_j)` and `m_j >= 1/4` by construction, so it can never exceed
    the nominal 9 -- exceeding it is how the N-as-a-fifth-base bug announced itself (9.01 nt for a
    9 nt barcode). Under it by a little is what a real, slightly biased library looks like: 8.94
    here against 8.97 on the full run.

    The barcode count is checked against the fixture's own tags, folding N to A, because that is
    what `pack_barcode` does: 628 distinct `RX` strings, 5 of which carry an N, are 623 barcodes.
    """
    sample = pipeline["sample"]
    assert sample["umi_length"] == 9
    assert 8.5 < sample["effective_length"] <= sample["umi_length"]
    assert sample["barcode_space"]["effective_space"] <= sample["barcode_space"]["nominal_space"]

    distinct = {tags["RX"].replace("N", "A") for _, tags, _, _ in records(READS)}
    assert sample["umis"] == len(distinct)


def test_the_error_budget_is_finite_and_no_larger_than_the_chemistry_predicts(pipeline):
    """What the reported Phred predicts, against what the distance-1 excess can still find.

    Two properties hold on this fixture:

    * `from_phred` is the mean of 10^(-Q/10), not 10^(-mean Q/10). On real quality the two differ
      by 5x here (1.7e-3 against 3.3e-4) because the low-Q tail carries the errors; on flat
      simulated quality they coincide, so this is a check only real data can make.
    * the estimate is positive, finite and does not exceed the prediction. It cannot be pinned
      tighter: the subsample dropped the error children, which the second half of this test
      measures directly rather than asserting on faith.
    """
    budget = pipeline["sample"]["error_budget"]
    assert 1e-4 < budget["from_phred"] < 1e-2
    assert budget["from_phred"] > 3 * 10 ** (-budget["mean_phred"] / 10)
    assert budget["predicted"] == pytest.approx(budget["from_phred"] + budget["from_polymerase"])
    assert 0.0 < budget["estimated"] <= budget["predicted"]
    # The distance-1 estimator collapses once most of a barcode's 3L shell is occupied. At 0.26%
    # it has not collapsed, so a zero excess is the fixture's, not the estimator's.
    assert budget["neighbour_occupancy"] < 0.05
    assert not budget["estimate_unreliable"]

    barcodes = sorted({tags["RX"].replace("N", "A") for _, tags, _, _ in records(READS)})
    observed = sum(
        1
        for a, b in itertools.combinations(barcodes, 2)
        if sum(x != y for x, y in zip(a, b)) == 1
    )
    # Chance alone, over the nominal 4^9 space: n(n-1)/2 pairs, each with 3L of 4^L neighbours.
    chance = len(barcodes) * (len(barcodes) - 1) / 2 * (3 * 9) / 4**9
    assert observed < 2 * chance, (
        "the fixture has real error children again -- it was regenerated by something other than "
        "an independent hash over barcodes. Tighten the assertion above to a ratio bound."
    )


def test_the_hiv_amplicon_reads_as_clonal(pipeline):
    """Payload agreement is worth log(1/clonality), so the clonality has to come from the data.

    Every molecule here is the same 268 nt of pol, which is what makes 0.80 the right answer and
    why a residual-FDR estimator that assumed diverse payload once called 97.4% of this library's
    singletons error children. Asserted as a property, not a constant: an amplicon is clonal, the
    10x T cells in `test_droplet.py` are not, and the same field has to separate them.
    """
    assert pipeline["refine"]["payload_clonality"] > 0.5


def test_correction_only_ever_removes_barcodes(pipeline):
    """A merge deletes a molecule, so the count may only go down, and by exactly what was merged.

    Reads are conserved through the rewrite -- correction moves a read to another barcode, it
    never drops one -- and the audit trail has one row per surviving barcode plus the merged ones,
    which is what makes a merge reviewable after the fact.
    """
    refine = pipeline["refine"]
    assert refine["reads"] == pipeline["sample"]["reads"]
    assert refine["reads_without_umi"] == 0
    assert refine["barcodes"] == pipeline["sample"]["umis"]
    assert refine["molecules"] == refine["barcodes"] - refine["merged"]
    assert refine["molecules"] <= refine["barcodes"]
    assert refine["merged_reads"] >= refine["merged"], "a merged barcode carries at least one read"
    # Collisions hide molecules, they never invent them, so the corrected count is the floor.
    assert refine["molecules_corrected"] >= refine["molecules"]
    assert not refine["saturated"], "623 barcodes in 4^9 is 0.26% occupancy"

    rewritten = sum(1 for _ in records(pipeline["dir"] / "ref" / "CTRL.fq.gz"))
    assert rewritten == refine["reads"]

    table = (pipeline["dir"] / "ref" / "CTRL.barcodes.tsv").read_text().splitlines()
    assert len(table) - 1 == refine["barcodes"]


def test_the_consensus_is_a_valid_fastq_capped_at_the_rt_floor(pipeline):
    """One record per molecule, and no emitted quality above the pre-amplification floor.

    An error made in reverse transcription or the first PCR cycle is in every read of the
    molecule, so no depth of coverage removes it: with `--rt-error rt` (1e-4) nothing may be
    emitted above Q40. Both halves matter -- the cap must hold, and it must be REACHED, or a
    consensus that silently lost its depth would satisfy it just as well.

    `cD` is the true depth, including reads past the 10,000 consensus cap, so it sums to the
    input. Splitting a group at the 8.68 linkage threshold makes molecules from groups, never
    the other way, so molecules >= groups.
    """
    assemble = pipeline["assemble"]
    assert assemble["reads"] == pipeline["refine"]["reads"]
    assert assemble["molecules"] >= assemble["groups"]
    assert assemble["quality_cap"] == pytest.approx(-10 * math.log10(assemble["rt_floor"]))

    cap = assemble["quality_cap"]
    names: set[str] = set()
    depth, highest = 0, 0
    for name, tags, sequence, quality in records(
        pipeline["dir"] / "asm" / "CTRL.consensus.fq.gz"
    ):
        assert len(sequence) == len(quality)
        # The name IS the molecule id -- `dnaio` drops FASTQ comments, so it has to stand alone,
        # and a split group suffixes it rather than repeating it.
        assert name == f"CTRL.{tags['RX']}" or name.startswith(f"CTRL.{tags['RX']}.")
        assert tags["MI"] == name
        assert name not in names, "two molecules under one id is a silently merged pair"
        names.add(name)
        best = max(ord(c) - 33 for c in quality)
        assert best <= cap
        highest = max(highest, best)
        depth += int(tags["cD"])
    assert len(names) == assemble["molecules"]
    assert depth == assemble["reads"]
    assert highest == cap, "at 15.8 reads per molecule the floor has to be reached somewhere"
