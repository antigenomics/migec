# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "migec", "polars"]
# ///
"""How many molecules does a variant caller actually get? -- on real ctDNA reference material.

Run with:  marimo edit notebooks/ctdna_variants.py

Fetches its own data from SRA (PRJNA788522, cell-free DNA reference material at known mutant
allele frequencies), runs the three migec stages, and answers the question that decides which
variant caller matters: is the variant present in enough molecules for any caller to see it?

Nothing here is stored in this repo or on HuggingFace. `scripts/sra_fetch.py` pulls the runs from
NCBI's S3 mirror on demand.
"""

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # Rare variants in ctDNA: the molecule count comes first

        Choosing a variant caller is a second-order decision. The first-order one is arithmetic
        that happens before any software runs:

        ```
        molecules at a site  =  input DNA / 3.3 pg  x  strands recovered  x  efficiency
        variant molecules    =  molecules at a site  x  VAF
        ```

        Deeper sequencing recovers more of the molecules that are in the tube, but only up to the
        number of input molecules -- past that it adds reads per molecule and nothing else. A
        caller cannot recover a molecule that was never sampled, and the supporting count is a
        Poisson draw rather than a guarantee: an expectation of 3 supporting molecules means a
        third of replicates see fewer than 3, whatever the caller.

        This notebook measures the left-hand side on data where the right-hand side is known.
        """
    )
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
        ## The data

        `PRJNA788522` is SiMSen-Seq on commercial cell-free DNA reference material, from the
        UMIErrorCorrect paper (Osterlund et al., Clin Chem 2022). The design crosses:

        * **mutant allele frequency** -- `WT` (0%, a true negative), 0.125%, 0.25%, 1%
        * **DNA input** -- 5, 20, 80 ng
        * **sequencing depth** -- 3.3x, 10x, 30x reads per UMI

        three replicates each. What makes it usable here and not in the published ctDNA caller
        benchmarks: the 12 nt inline UMI **survived deposition**. `migec suggest` finds it from
        base composition alone, with no knowledge of the protocol.
        """
    )
    return


@app.cell
def _():
    import subprocess
    import sys
    import tempfile
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    work = Path(tempfile.mkdtemp(prefix="migec-ctdna-"))

    # One replicate per arm keeps the download to a few hundred MB. The committed table in
    # assets/ has the full 3-replicate design.
    RUNS = {
        "SRR17220921": "20ng_WT_10x",
        "SRR17220931": "20ng_0.125_10x",
        "SRR17220928": "20ng_0.25_10x",
        "SRR17220924": "20ng_1_10x",
        "SRR17220957": "5ng_0.125_10x",
    }

    subprocess.run(
        [sys.executable, str(repo / "scripts" / "sra_fetch.py"), "get", *RUNS, "-o", str(work)],
        check=True,
    )
    sorted(p.name for p in work.glob("*.fastq.gz"))
    return Path, RUNS, repo, subprocess, sys, work


@app.cell
def _(mo):
    mo.md(
        """
        ## Where is the barcode?

        Never ask the protocol -- ask the data. `suggest` segments the per-cycle base composition
        into UMI, constant and payload runs. On these reads it reports a 12 nt near-uniform run at
        cycles 0-11, then the SiMSen-Seq constant spacer.
        """
    )
    return


@app.cell
def _(RUNS, work):
    from migec.suggest import run as suggest_run

    first = work / f"{next(iter(RUNS))}_1.fastq.gz"
    layout = suggest_run(first, reads=50_000)
    layout["segments"]
    return first, layout, suggest_run


@app.cell
def _(mo):
    mo.md(
        """
        ## Run the three stages

        `checkout` extracts and trims the UMI, `refine` corrects barcode errors and counts
        molecules, `assemble` builds one consensus record per molecule.
        """
    )
    return


@app.cell
def _(RUNS, work):
    import polars as pl

    from migec import assemble, checkout, refine
    from migec.cli import _inline_sheet
    from migec.sheet import parse_layout

    pattern, _ = parse_layout("0:12")
    rows = []
    for run_acc, alias in RUNS.items():
        w = work / run_acc
        (w / "co").mkdir(parents=True, exist_ok=True)
        sheet = _inline_sheet(pattern, run_acc, w / "co")
        checkout.run(work / f"{run_acc}_1.fastq.gz", sheet, w / "co", max_offset=0)
        rf = refine.run(w / "co" / f"{run_acc}.fq.gz", w / "rf")
        asm = assemble.run(w / "rf" / f"{run_acc}.fq.gz", w / "as")
        rows.append({
            "run": run_acc, "arm": alias,
            "reads": rf["reads"], "molecules": asm["molecules"],
            "reads_per_molecule": round(rf["reads"] / asm["molecules"], 1),
            "barcode_phred": round(rf["error_phred"], 1),
            "consensus_phred": round(asm["mean_quality"], 1),
        })
    stages = pl.DataFrame(rows)
    stages
    return (
        _inline_sheet, alias, asm, assemble, checkout, parse_layout, pattern,
        pl, refine, rf, rows, run_acc, sheet, stages, w,
    )


@app.cell
def _(mo):
    mo.md(
        """
        ## Molecules per site, not molecules per library

        Never: **the molecule total of a multiplex panel is not the count at any one site.** A
        variant sits on one amplicon, so the evidence a caller gets is `molecules / amplicons`.
        Reporting the total overstates it by exactly the panel size.

        The amplicon count is measured, not assumed: every molecule of one amplicon starts with
        the same bases, so a prefix tally finds them with no reference and no aligner. The share
        threshold has to sit in the **gap** between the smallest real amplicon and the largest
        payload-error prefix -- put it on the slope and the count climbs at shallow depth, exactly
        where the evidence is thinnest.
        """
    )
    return


@app.cell
def _(RUNS, pl, repo, work):
    import sys as _sys

    # `scripts/` itself, not the repo root: it has no `__init__.py`, so `scripts.ctdna_titration`
    # is not an importable path. Putting the directory itself on sys.path is the honest fix --
    # these are runnable scripts, not a package, and making them one for a notebook import would
    # be the tail wagging the dog.
    _sys.path.insert(0, str(repo / "scripts"))
    from ctdna_titration import count_amplicons  # noqa: E402

    shares = []
    for _acc in RUNS:
        cons = work / _acc / "as" / f"{_acc}.consensus.fq.gz"
        n, share = count_amplicons(cons, prefix_len=25, min_share=0.05)
        shares.append({"run": _acc, "amplicons": n, "share_of_molecules": round(share, 3)})
    pl.DataFrame(shares)
    return _acc, _sys, cons, count_amplicons, n, share, shares


@app.cell
def _(mo):
    mo.md(
        """
        ## The answer

        `variant_molecules` is what the caller gets; `p_enough` is the probability that at least
        three of them are actually present, Poisson in the expectation. Below about 0.9 the row is
        a coin flip and the choice of caller stops mattering -- the evidence is not there.

        The `WT` row is the control: its true variant count is zero, so anything a caller reports
        there is its own false-positive rate on real chemistry rather than on a simulation.
        """
    )
    return


@app.cell
def _(count_amplicons, pl, shares, stages, work):
    import math

    VAF = {"20ng_WT_10x": 0.0, "20ng_0.125_10x": 0.00125, "20ng_0.25_10x": 0.0025,
           "20ng_1_10x": 0.01, "5ng_0.125_10x": 0.00125}
    amps = {s["run"]: s["amplicons"] for s in shares}

    out = []
    for r in stages.iter_rows(named=True):
        per_amp = r["molecules"] / amps[r["run"]]
        lam = per_amp * VAF[r["arm"]]
        below = sum(math.exp(-lam) * lam**k / math.factorial(k) for k in range(3))
        out.append({
            "arm": r["arm"],
            "molecules_per_amplicon": round(per_amp),
            "vaf": VAF[r["arm"]],
            "variant_molecules": round(lam, 1),
            "p_at_least_3": round(1 - below, 4),
        })
    pl.DataFrame(out)
    return VAF, amps, below, lam, math, out, per_amp, r


@app.cell
def _(mo):
    mo.md(
        """
        ## What to run next

        The consensus is one record per molecule, so a **standard** caller reads it correctly and
        its depth already is a molecule count:

        ```bash
        minimap2 -ax sr -y ref.fa as/S1.consensus.fq.gz | samtools sort -o S1.bam
        lofreq call -f ref.fa -o S1.vcf S1.bam
        ```

        Never feed it to a UMI-aware caller (`UMI-VarCal`, `UMIErrorCorrect`). Those group and
        consense themselves, so they replace `assemble` rather than following it -- and a
        `--min-family-size 3` filter discards the whole library, because after `assemble` every
        family has size 1 by construction. `docs/variants.rst` has the full table.
        """
    )
    return


if __name__ == "__main__":
    app.run()
