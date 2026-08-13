# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "migec", "altair", "polars"]
# ///
"""UMI diagnostics: the coverage curve, the barcode-rank plot, and where the errors are.

Run with:  marimo edit notebooks/refine_diagnostics.py

Every figure is drawn from a TSV that `migec refine` wrote, never recomputed here, so any of them
can be redrawn from a committed table long after the run.
"""

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # What the UMIs actually look like

        Four questions, each with a table `refine` already wrote:

        1. **How deeply was each molecule sequenced?** The coverage histogram
           (`refine.coverage.tsv`). Below ~3 reads per UMI the consensus buys counting rather
           than error correction, and barcode correction loses most of its evidence.
        2. **Is there a knee?** The barcode-rank curve (`<sample>.rank.tsv`), which is Cell
           Ranger's plot: reads per barcode against rank, both on log axes.
        3. **Where are the errors?** `fraction_erroneous` per size bin (`<sample>.bins.tsv`).
           Error children pile up at one read. Finding them at high counts means correction is
           eating real molecules.
        4. **Is a bin one artefact repeated, or a population?** Payload entropy per bin, in the
           same table.
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

    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run

    ADAPTER = "CAGTGGTATCAACGCAGAGT"
    tmp = Path(tempfile.mkdtemp())

    def build(coverage, tag):
        """Simulate a library at one depth, check it out, refine it. Returns the output dir."""
        cfg = SimConfig(
            n_molecules=20_000, n_clones=100, coverage=coverage, coverage_cv=0.5,
            umi_len=12, umi_error=3e-3, adapter=ADAPTER, seed=5,
        )
        sim = simulate(cfg, tmp / f"sim_{tag}")
        (tmp / f"bc_{tag}.txt").write_text(f"S1\t{sim['pattern']}\n")
        checkout_run(sim["reads"], tmp / f"bc_{tag}.txt", tmp / f"co_{tag}")
        out = tmp / f"ref_{tag}"
        summary = refine_run(tmp / f"co_{tag}" / "S1.fq.gz", out)
        return out, summary

    deep, deep_summary = build(8.0, "deep")
    shallow, shallow_summary = build(1.5, "shallow")
    return Path, deep, deep_summary, shallow, shallow_summary


@app.cell
def _(deep, mo, shallow):
    import polars as pl

    def bins(path, label):
        return pl.read_csv(path / "S1.bins.tsv", separator="\t").with_columns(
            pl.lit(label).alias("library")
        )

    binned = pl.concat([bins(deep, "8 reads/UMI"), bins(shallow, "1.5 reads/UMI")])
    mo.md("### The coverage histogram, and what fraction of each bin is error")
    return binned, pl


@app.cell
def _(binned):
    binned.select(
        ["library", "min_reads", "max_reads", "barcodes", "reads", "merged",
         "fraction_erroneous", "payload_entropy_bits"]
    )
    return


@app.cell
def _(binned, mo):
    import altair as alt

    chart = (
        alt.Chart(binned.to_pandas())
        .mark_bar()
        .encode(
            x=alt.X("min_reads:O", title="reads per barcode (power-of-two bin)"),
            y=alt.Y("fraction_erroneous:Q", title="fraction merged as error",
                    axis=alt.Axis(format="%")),
            color=alt.Color("library:N", title=None),
            xOffset="library:N",
        )
        .properties(width=460, height=260, title="Error children pile up at one read")
    )
    mo.ui.altair_chart(chart)
    return (alt,)


@app.cell
def _(mo):
    mo.md(
        """
        That shape is the diagnostic. A barcode seen **once** is overwhelmingly an error child of
        one seen many times; a barcode seen eight times almost never is. If the curve is flat, or
        rises at high counts, correction is merging real molecules rather than errors — and a
        merge deletes a molecule that nothing downstream can recover.
        """
    )
    return


@app.cell
def _(alt, deep, mo, pl, shallow):
    def rank(path, label):
        return pl.read_csv(path / "S1.rank.tsv", separator="\t").with_columns(
            pl.lit(label).alias("library")
        )

    ranked = pl.concat([rank(deep, "8 reads/UMI"), rank(shallow, "1.5 reads/UMI")])
    mo.md("### The barcode-rank curve")
    mo.ui.altair_chart(
        alt.Chart(ranked.to_pandas())
        .mark_line()
        .encode(
            x=alt.X("rank:Q", scale=alt.Scale(type="log"), title="barcode rank"),
            y=alt.Y("reads:Q", scale=alt.Scale(type="log"), title="reads"),
            color=alt.Color("library:N", title=None),
        )
        .properties(width=460, height=280, title="Reads per barcode against rank")
    )
    return (ranked,)


@app.cell
def _(alt, mo, ranked):
    mo.md("### ...and its CDF, which is where over-sequencing is visible")
    mo.ui.altair_chart(
        alt.Chart(ranked.to_pandas())
        .mark_line()
        .encode(
            x=alt.X("rank:Q", scale=alt.Scale(type="log"), title="barcode rank"),
            y=alt.Y("cumulative_fraction:Q", title="cumulative share of reads",
                    axis=alt.Axis(format="%")),
            color=alt.Color("library:N", title=None),
        )
        .properties(width=460, height=260, title="Cumulative reads by barcode rank")
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        A curve that rises steeply and then flattens means a few barcodes carry most of the reads
        — over-amplification, or a barcode too short to separate molecules. A straight diagonal
        means the reads are spread evenly, which is what a well-mixed library at low depth looks
        like. The rank table is log-spaced, so it is a few hundred rows however large the library.
        """
    )
    return


@app.cell
def _(alt, binned, mo):
    mo.md(
        """
        ### Sequence entropy per bin

        Per-position Shannon entropy of the payload drafts of the barcodes in each bin. Two bits
        is a uniform base. A bin that drops well below the rest is holding **one sequence
        repeated** — a clone, or an artefact like a primer dimer — rather than a population of
        distinct molecules, and it is worth looking at before it becomes a result.
        """
    )
    mo.ui.altair_chart(
        alt.Chart(binned.to_pandas())
        .mark_line(point=True)
        .encode(
            x=alt.X("min_reads:O", title="reads per barcode (power-of-two bin)"),
            y=alt.Y("payload_entropy_bits:Q", title="entropy (bits/base)",
                    scale=alt.Scale(zero=False)),
            color=alt.Color("library:N", title=None),
        )
        .properties(width=460, height=240, title="Payload entropy by MIG size")
    )
    return


@app.cell
def _(deep_summary, mo, shallow_summary):
    mo.md(
        f"""
        ### The two runs, side by side

        | | 8 reads/UMI | 1.5 reads/UMI |
        |---|---|---|
        | barcodes in | {deep_summary['barcodes']:,} | {shallow_summary['barcodes']:,} |
        | merged as error | {deep_summary['merged']:,} | {shallow_summary['merged']:,} |
        | molecules out | {deep_summary['molecules']:,} | {shallow_summary['molecules']:,} |
        | barcode error estimate | {deep_summary['estimated_error']:.2e} |
          {shallow_summary['estimated_error']:.2e} |
        | clonality | {deep_summary['payload_clonality']:.4f} |
          {shallow_summary['payload_clonality']:.4f} |

        20,000 molecules were simulated in both. The deep run recovers them; the shallow one
        cannot, because most of its barcode errors have no observable parent — and it says so
        rather than reporting a corrected count it has no evidence for.
        """
    )
    return


if __name__ == "__main__":
    app.run()
