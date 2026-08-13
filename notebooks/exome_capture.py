# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "migec", "polars"]
# ///
"""Duplicates or molecules? UMI-aware deduplication on a hybrid-capture library.

Run with:  uv run marimo edit notebooks/exome_capture.py

A capture panel or exome is the case where position-based deduplication and UMI-based grouping
disagree most visibly, because capture probes make reads pile up on the same coordinates whether
or not they came from the same molecule. This notebook builds a library where the truth is known
and shows what each answer costs.

Simulated on purpose: the point is a controlled comparison against a known molecule count, and no
public exome deposits its UMIs (the same gap `docs/variants.rst` documents for ctDNA).
"""

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # Capture panels: coordinates lie, barcodes do not

        Deduplicating a hybrid-capture library by alignment coordinate rests on an assumption:
        *two reads starting at the same base are copies of one molecule.* Capture is exactly where
        that assumption is weakest. Probes tile fixed positions, fragment ends cluster on a handful
        of preferred cut sites, and at any useful depth two independent molecules routinely share a
        start.

        Coordinate deduplication cannot tell those apart, so it **deletes real molecules** and
        undercounts. A UMI can tell them apart, which is the entire reason capture panels carry one.

        Three numbers say whether it matters for your library:

        1. **How many molecules share a start position** -- set by depth and fragment diversity.
        2. **How long is the barcode** -- whether two molecules can collide on the UMI as well.
        3. **How often is the barcode miscalled** -- an error splits one molecule into two.

        migec measures all three. Below, on a library whose true molecule count is known.
        """
    )
    return (mo,)


@app.cell
def _():
    import sys
    import tempfile
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from tests.synthetic._sim import SimConfig, simulate

    from migec.assemble import run as assemble_run
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run

    # A capture-like shape: a few thousand molecules over a small panel, ~12 reads each, and a
    # 12 nt UMI -- the length IDT, Twist and Roche all use for cfDNA and exome kits.
    ADAPTER = "CAGTGGTATCAACGCAGAGT"
    tmp = Path(tempfile.mkdtemp(prefix="migec-capture-"))

    cfg = SimConfig(
        n_molecules=30_000, n_clones=250, coverage=12.0, coverage_cv=0.8,
        umi_len=12, umi_error=1e-3, adapter=ADAPTER, seed=17,
    )
    sim = simulate(cfg, tmp / "sim")
    (tmp / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")

    co = checkout_run(sim["reads"], tmp / "bc.txt", tmp / "co")
    rf = refine_run(tmp / "co" / "S1.fq.gz", tmp / "rf")
    asm = assemble_run(tmp / "rf" / "S1.fq.gz", tmp / "as")

    TRUE_MOLECULES = cfg.n_molecules
    return (
        ADAPTER, Path, SimConfig, TRUE_MOLECULES, asm, assemble_run, cfg, co,
        checkout_run, refine_run, rf, sim, simulate, sys, tempfile, tmp,
    )


@app.cell
def _(mo):
    mo.md(
        """
        ## What the barcode recovered

        `refine` corrects barcode errors and counts molecules; `assemble` builds one consensus
        record per molecule. Against a simulator that knows the answer, the interesting column is
        the recovery, not the raw count.
        """
    )
    return


@app.cell
def _(TRUE_MOLECULES, asm, co, rf):
    import polars as pl

    pl.DataFrame([
        {"quantity": "reads", "value": rf["reads"]},
        {"quantity": "distinct barcodes seen", "value": rf["barcodes"]},
        {"quantity": "barcodes merged as errors", "value": rf["merged"]},
        {"quantity": "molecules called", "value": asm["molecules"]},
        {"quantity": "molecules simulated (truth)", "value": TRUE_MOLECULES},
        {"quantity": "recovery", "value": round(asm["molecules"] / TRUE_MOLECULES, 4)},
        {"quantity": "barcode error, Phred", "value": round(rf["error_phred"], 2)},
        {"quantity": "reads per molecule", "value": round(rf["reads"] / asm["molecules"], 2)},
    ])
    return (pl,)


@app.cell
def _(mo):
    mo.md(
        """
        ## The collision arithmetic, which is what coordinate deduplication ignores

        `checkout` reports how full the barcode space is. Occupancy is what decides whether two
        molecules could have drawn the same UMI -- the birthday problem, not intuition. Note that
        the usable space is the **collision (Renyi-2) entropy** of the observed base composition,
        not `4^L`: a real oligo mix is skewed, and using the nominal space overstates how many
        distinct barcodes there are and so understates collisions.
        """
    )
    return


@app.cell
def _(co, pl):
    b = co["samples"][0]["barcode_space"]
    pl.DataFrame([
        {"quantity": "nominal space (4^L)", "value": f"{b['nominal_space']:,.0f}"},
        {"quantity": "effective space (collision entropy)", "value": f"{b['effective_space']:,.0f}"},
        {"quantity": "occupancy", "value": f"{100 * b['occupancy']:.2f}%"},
        {"quantity": "MIGs holding >1 molecule", "value": f"{100 * b['p_multi']:.2f}%"},
        {"quantity": "molecules hidden by collision", "value": f"{b['hidden']:,.0f}"},
        {"quantity": "saturated?", "value": str(bool(b["saturated"]))},
    ])
    return (b,)


@app.cell
def _(mo):
    mo.md(
        """
        ## What this means for the pipeline

        **Do not run a coordinate deduplicator on this data, and do not run one after migec.**
        `MarkDuplicates` without UMI awareness collapses independent molecules that share a start;
        with `BARCODE_TAG=RX` it uses the UMI and becomes redundant with `assemble`. After
        `assemble` there is exactly one record per molecule, so any further deduplication removes
        real molecules.

        The order that works:

        ```bash
        migec checkout R1.fq.gz R2.fq.gz --bc-pattern '0:12' --sample S1 -o co/
        migec refine   co/S1.fq.gz -o rf/
        migec assemble rf/S1.fq.gz -o as/
        bwa mem -C ref.fa as/S1.consensus.fq.gz | samtools sort -o S1.bam   # RX/MI carried through
        # then call variants with a STANDARD caller -- depth already is a molecule count
        ```

        `docs/downstream.rst` has the measured tag-survival table, and `docs/variants.rst` covers
        which variant caller composes with this and which replaces a stage of it.
        """
    )
    return


@app.cell
def _(mo, tmp):
    import shutil

    shutil.rmtree(tmp, ignore_errors=True)
    mo.md(f"Cleaned up `{tmp}`.")
    return (shutil,)


if __name__ == "__main__":
    app.run()
