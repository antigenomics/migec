# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "migec", "polars"]
# ///
"""How much of a repertoire is PCR? Error-corrected clonotypes, end to end.

Run with:  uv run marimo edit notebooks/airr_repertoire.py

Immune repertoire sequencing is the case migec was built for, and the one where the answer is a
*count* rather than a genotype: how many distinct clones are there, and how abundant is each. That
makes it the case most damaged by amplification error, because every PCR artefact looks exactly
like a rare new clonotype -- and rare clonotypes are the biology.

This notebook simulates a repertoire with a known clone count, runs the three stages, and shows
what the barcode recovers that read counting cannot.
"""

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # Repertoires: every PCR error is a plausible new clone

        In a repertoire library the deliverable is a list of clonotypes and their abundances.
        Two things make that harder than it sounds:

        1. **A PCR or sequencing error produces a sequence that is a valid clonotype.** There is no
           reference to contradict it, and it will differ from a real clone by one base -- exactly
           like a genuinely related clone does. Counting reads counts the artefacts.
        2. **Abundance is amplification, not biology.** Two clones present in equal numbers can
           differ ten-fold in reads because one primer bound better.

        A UMI answers both: collapsing reads to molecules removes errors made after barcoding, and
        counting *molecules* rather than reads gives an abundance that is not an amplification
        artefact. What it cannot remove is an error made **before** the barcode was attached -- the
        RT/first-cycle floor -- which is why `assemble` caps emitted quality rather than claiming
        Q60. See `docs/quality_floor.rst`.
        """
    )
    return (mo,)


@app.cell
def _(mo):
    mo.md(
        """
        ## A library with a known answer

        The simulator draws a clone size distribution, amplifies it with per-cycle error, attaches
        a UMI, and sequences it. `n_clones` is the truth we will check against. The adapter is the
        SMART sequence from the original MIGEC protocol, so the layout is the real one.
        """
    )
    return


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

    ADAPTER = "CAGTGGTATCAACGCAGAGT"
    tmp = Path(tempfile.mkdtemp(prefix="migec-airr-"))

    # A shallow bulk repertoire: many molecules, few reads each. This is the ORDINARY case for
    # repertoire profiling, not an exotic one -- see `docs/umi_statistics.rst`.
    cfg = SimConfig(
        n_molecules=40_000, n_clones=400, coverage=3.0, coverage_cv=0.9,
        umi_len=12, umi_error=1e-3, pcr_error=2e-5, adapter=ADAPTER, seed=23,
    )
    sim = simulate(cfg, tmp / "sim")
    (tmp / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")

    co = checkout_run(sim["reads"], tmp / "bc.txt", tmp / "co")
    rf = refine_run(tmp / "co" / "S1.fq.gz", tmp / "rf")
    asm = assemble_run(tmp / "rf" / "S1.fq.gz", tmp / "as")
    return (
        ADAPTER, Path, SimConfig, asm, assemble_run, cfg, co, checkout_run,
        refine_run, rf, sim, simulate, sys, tempfile, tmp,
    )


@app.cell
def _(asm, cfg, rf):
    import polars as pl

    pl.DataFrame([
        {"quantity": "reads", "value": rf["reads"]},
        {"quantity": "reads per molecule", "value": round(rf["reads"] / asm["molecules"], 2)},
        {"quantity": "distinct barcodes seen", "value": rf["barcodes"]},
        {"quantity": "barcodes merged as errors", "value": rf["merged"]},
        {"quantity": "molecules called", "value": asm["molecules"]},
        {"quantity": "molecules simulated", "value": cfg.n_molecules},
        {"quantity": "barcode error, Phred", "value": round(rf["error_phred"], 2)},
        {"quantity": "clonality (measured)", "value": rf["payload_clonality"]},
    ])
    return (pl,)


@app.cell
def _(mo):
    mo.md(
        """
        ## Reads are not molecules, and the difference is the point

        Below: distinct sequences counted from the **raw reads**, against distinct sequences
        counted from the **consensus**. The simulated truth is 400 clones. Read counting inflates
        that number by every amplification error that survived; the consensus does not.
        """
    )
    return


@app.cell
def _(ADAPTER, cfg, pl, sim, tmp):
    import gzip
    from collections import Counter

    def distinct_sequences(path, skip):
        """Count distinct payloads, ignoring the leading barcode+adapter region."""
        seen = Counter()
        with gzip.open(path, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 1:
                    seen[line.strip()[skip:skip + 60]] += 1
        return seen

    # The raw read still carries UMI + adapter in front of the payload; the consensus does not,
    # because `checkout` trimmed the pattern off. Note: the simulator's return dict has no
    # `adapter` key -- reading one back from it silently gave a 12 nt offset instead of 32 and
    # compared two different windows of the read.
    skip_raw = cfg.umi_len + len(ADAPTER)
    raw = distinct_sequences(sim["reads"], skip_raw)
    cons = distinct_sequences(tmp / "as" / "S1.consensus.fq.gz", 0)

    pl.DataFrame([
        {"counted from": "raw reads", "distinct sequences": len(raw),
         "inflation vs truth": round(len(raw) / cfg.n_clones, 1)},
        {"counted from": "consensus (molecules)", "distinct sequences": len(cons),
         "inflation vs truth": round(len(cons) / cfg.n_clones, 1)},
        {"counted from": "truth", "distinct sequences": cfg.n_clones, "inflation vs truth": 1.0},
    ])
    return Counter, cons, distinct_sequences, gzip, raw, skip_raw


@app.cell
def _(mo):
    mo.md(
        """
        Note: neither number lands exactly on 400, and both are informative about *why*. Read
        counting inflates because every surviving error is a new sequence. The consensus still
        inflates a little, because an error made **before** barcoding is in every read of that
        molecule and no consensus can remove it -- that is the RT/first-cycle floor, and it is the
        residual a repertoire study actually has to live with.

        This is also why `assemble` never emits a quality above the floor: claiming Q60 for a
        single-molecule consensus would tell a downstream clustering tool that a pre-amplification
        error is a certainty.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## Handing it to AIRR tools

        One consensus record is one molecule, and its identity is in the read **name**
        (`<sample>.<cell>.<umi>`) as well as in SAM tags. That matters here specifically: `dnaio`,
        which most Python AIRR tooling uses, drops FASTQ comments -- so the name has to stand
        alone, and it does.

        ```bash
        migec checkout R1.fq.gz R2.fq.gz --preset migec --sample S1 -o co/
        migec refine   co/S1.fq.gz -o rf/
        migec assemble rf/S1.fq.gz -o as/

        arda amplicon --r1 as/S1.consensus.fq.gz -p S1     # AIRR sequence_id IS the molecule id
        ```

        The resulting AIRR `duplicate_count` is then a **molecule** count, not a read count, which
        is the number a repertoire diversity estimate should be computed from. `docs/downstream.rst`
        has the measured table of which tools preserve the tags.

        Never run a second deduplicator on this output. After `assemble` there is exactly one
        record per molecule, so anything that collapses again deletes real clones -- and in a
        repertoire the rare clones it deletes are usually the ones being looked for.
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
