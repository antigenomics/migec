# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "migec", "altair", "polars"]
# ///
"""Barcode space, collisions and the error budget -- the numbers checkout computes rather than counts.

Run with:  marimo edit notebooks/barcode_space.py

Everything here is drawn from checkout's own TSVs, never recomputed, so a figure can always be
regenerated from a committed table.
"""

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # Barcode space, collisions, and the error budget

        A UMI is a lottery ticket. Two questions decide whether the draw was any good:

        1. **How many tickets are there really?** `4^L` assumes a perfect synthesiser. Real oligo
           mixes are skewed, and the usable space is `1 / prod_j sum_a p_j(a)^2` -- the *collision*
           (Renyi-2) entropy, not Shannon.
        2. **How many molecules drew the same ticket?** That is the birthday problem, and past a
           few percent occupancy it stops being a curiosity and starts merging molecules.

        This notebook builds a library, runs `checkout` on it, and reads both off the output.
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

    from migec.checkout import run

    ADAPTER = "CAGTGGTATCAACGCAGAGT"
    tmp = Path(tempfile.mkdtemp())

    def build(umi_len, n_molecules, freqs=(0.25, 0.25, 0.25, 0.25), umi_error=0.0, tag=""):
        """Simulate a library and check it out. Returns checkout's per-sample summary."""
        name = f"{tag or umi_len}_{n_molecules}"
        cfg = SimConfig(
            n_molecules=n_molecules, n_clones=1, umi_len=umi_len, coverage=6.0,
            umi_base_freqs=freqs, umi_error=umi_error, adapter=ADAPTER, seed=7,
        )
        sim = simulate(cfg, tmp / f"sim_{name}")
        (tmp / f"bc_{name}.txt").write_text(f"S1\t{'N' * umi_len}{ADAPTER.lower()}\n")
        return run(sim["reads"], tmp / f"bc_{name}.txt", tmp / f"out_{name}")["samples"][0]

    return ADAPTER, Path, build, tmp


@app.cell
def _(build, mo):
    # Same molecule count, shrinking barcode. Nothing else changes.
    ladder = {L: build(L, 3000) for L in (12, 10, 8, 7, 6)}
    mo.md("### The same 3,000 molecules, in barcodes of different lengths")
    return (ladder,)


@app.cell
def _(ladder):
    import polars as pl

    space = pl.DataFrame(
        [
            {
                "umi_nt": L,
                "nominal_space": int(s["barcode_space"]["nominal_space"]),
                "effective_space": round(s["barcode_space"]["effective_space"]),
                "observed": s["barcode_space"]["observed"],
                "occupancy": s["barcode_space"]["occupancy"],
                "migs_over_one_molecule": s["barcode_space"]["p_multi"],
                "molecules_implied": round(s["barcode_space"]["molecules"]),
            }
            for L, s in ladder.items()
        ]
    )
    space
    return pl, space


@app.cell
def _(mo, pl, space):
    import altair as alt

    long = space.select(
        ["umi_nt", "occupancy", "migs_over_one_molecule"]
    ).unpivot(index="umi_nt", variable_name="metric", value_name="fraction")

    chart = (
        alt.Chart(long.to_pandas())
        .mark_line(point=True)
        .encode(
            x=alt.X("umi_nt:O", title="UMI length (nt)", sort="descending"),
            y=alt.Y("fraction:Q", title="fraction", axis=alt.Axis(format="%")),
            color=alt.Color("metric:N", title=None),
        )
        .properties(width=420, height=260, title="Occupancy, and the MIGs it ruins")
    )
    mo.ui.altair_chart(chart)
    return (alt,)


@app.cell
def _(mo):
    mo.md(
        """
        `migs_over_one_molecule` is the number to read. It is `P(k > 1 | k >= 1)` for a
        Poisson-occupied space, and it says what fraction of MIGs are a *mixture of templates*
        rather than one molecule. A consensus over those is not a molecule's sequence, and no
        amount of over-sequencing repairs it.

        Note the two curves separate: occupancy grows smoothly, but the damage grows faster,
        because the birthday problem is quadratic before it saturates.
        """
    )
    return


@app.cell
def _(build, mo):
    # A skewed synthesiser mix, at fixed length. This is what `bias_loss` measures.
    even = build(9, 3000, tag="even")
    skew = build(9, 3000, freqs=(0.40, 0.10, 0.40, 0.10), tag="skew")
    mo.md(
        f"""
        ### The synthesiser costs you barcode space

        Same 9 nt barcode, same 3,000 molecules, different oligo mix:

        | | even (25/25/25/25) | skewed (40/10/40/10) |
        |---|---|---|
        | nominal space | {even['barcode_space']['nominal_space']:,.0f} | {skew['barcode_space']['nominal_space']:,.0f} |
        | effective space | {even['barcode_space']['effective_space']:,.0f} | {skew['barcode_space']['effective_space']:,.0f} |
        | effective length | {even['barcode_space']['effective_length']:.2f} nt | {skew['barcode_space']['effective_length']:.2f} nt |
        | space lost to bias | {even['barcode_space']['bias_loss']:.1%} | **{skew['barcode_space']['bias_loss']:.1%}** |
        | MIGs holding >1 molecule | {even['barcode_space']['p_multi']:.2%} | **{skew['barcode_space']['p_multi']:.2%}** |

        The skewed library is nominally the same length and collides
        {skew['barcode_space']['p_multi'] / max(even['barcode_space']['p_multi'], 1e-9):.1f}x more.
        `migec suggest` sees this directly -- it is the shape of the 1/4 PWM trace.
        """
    )
    return even, skew


@app.cell
def _(alt, even, mo, skew):
    import polars as pl2

    comp = pl2.DataFrame(
        [
            {"mix": label, "position": p["position"], "base": b, "frequency": p[b]}
            for label, s in (("even", even), ("skewed", skew))
            for p in s["composition"]
            for b in "ACGT"
        ]
    )
    mo.ui.altair_chart(
        alt.Chart(comp.to_pandas())
        .mark_bar()
        .encode(
            x=alt.X("position:O", title="barcode position"),
            y=alt.Y("frequency:Q", stack="normalize", axis=alt.Axis(format="%")),
            color=alt.Color("base:N", scale=alt.Scale(scheme="category10")),
            row=alt.Row("mix:N", title=None),
        )
        .properties(width=420, height=110, title="The 1/4 PWM trace, even and skewed")
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ### The error budget

        `checkout` estimates the barcode error rate from the excess of barcode pairs at Hamming
        distance 1 -- a measurement. Two things predict what it should be:

        - **sequencing**, from the reported Phred: the mean of `10^(-Q/10)`, *not*
          `10^(-mean Q/10)`. The function is convex, so the low-Q tail carries nearly all the
          error and averaging Q first hides it.
        - **polymerase**, `eps_pol x cycles`.

        When the measurement and the prediction disagree by more than 3x, one of them is wrong.
        """
    )
    return


@app.cell
def _(build, pl):
    # Same injected barcode error, shrinking space: the estimator has a working range.
    sweep = pl.DataFrame(
        [
            {
                "umi_nt": L,
                "occupancy": s["barcode_space"]["occupancy"],
                "neighbourhood_occupied": s["error_budget"]["neighbour_occupancy"],
                "predicted": s["error_budget"]["predicted"],
                "estimated": s["error_budget"]["estimated"],
                "flagged_unreliable": s["error_budget"]["estimate_unreliable"],
            }
            for L in (12, 10, 9, 8)
            for s in [build(L, 20000, umi_error=3e-3, tag=f"err{L}")]
        ]
    )
    sweep
    return (sweep,)


@app.cell
def _(mo):
    mo.md(
        """
        The estimate collapses as the space fills, and always **downward** -- a crowded library
        reports too little barcode error and under-corrects. That is why `estimate_unreliable` is
        set past 5% neighbourhood occupancy rather than the number being quietly reported.

        The fix is not a better estimator, it is a longer barcode.
        """
    )
    return


if __name__ == "__main__":
    app.run()
