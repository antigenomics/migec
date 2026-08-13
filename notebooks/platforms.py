# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "migec", "polars", "altair", "huggingface-hub"]
# ///
"""Every barcode layout migec has been tested against, and how to declare it.

Run with:  marimo edit notebooks/platforms.py

The only thing that changes between platforms is where the barcode is. Everything after that -
correction, consensus, the quality cap - is the same three commands. This notebook shows the
declaration for each layout we have run, on real data where the data is public.

Fixtures come from https://huggingface.co/datasets/isalgo/umi_data and are small on purpose.
"""

import marimo

__generated_with = "0.9.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo

    mo.md(
        """
        # Which layout is my library, and what do I type?

        Three commands, always the same:

        ```
        migec checkout  reads -> tagged reads     # find and cut out the barcode
        migec refine    tagged -> corrected       # fix errors IN the barcode
        migec assemble  corrected -> consensus    # collapse each molecule
        ```

        Only `checkout` needs to know anything about the platform, and it needs exactly one thing:
        **where the barcode is**. Four ways to say it, in the order to reach for them.

        | | what it looks like | when |
        |---|---|---|
        | a position | `^NNNNNNNN` or `0:8` | the primary mode - the barcode is at an offset |
        | `--preset` | `10x-v2`, `tso500`, `duplex`, ... | a chemistry with a name |
        | `--read-structure` | `5M5S+T` | fgbio / Picard / TSO500 speak this |
        | barcode table | `S1<TAB>aaACTcagtgg...NNNNtNNNNtNNNN` | many samples in one file |

        In a pattern: **`N`** is a UMI base, **`X`** a cell-barcode base, **uppercase** is matched
        exactly, **lowercase** is the fuzzy adapter, **`.`** is skipped. Slices are half-open and
        0-based like Python's, each a UMI slice unless prefixed `cell:`, so `0:8` is eight bases
        and 10x is `cell:0:16,16:26`.

        A leading `^`, a slice list and a read structure all **anchor the barcode at the first
        base**. That used to be `--max-offset 0`, typed by hand, and forgetting it made 10x fail
        with an error about anchoring. It is now automatic and the flag should not be passed.

        If you do not know the layout, do not guess:

        ```bash
        migec suggest reads.fq.gz
        ```
        """
    )
    return (mo,)


@app.cell
def _():
    from pathlib import Path

    from migec.sheet import from_read_structure

    def show(structure):
        return f"`{structure}` -> `{from_read_structure(structure)}`"

    return Path, from_read_structure, show


@app.cell
def _(mo, show):
    mo.md(
        """
        ## 1. Bulk amplicon with an anchor

        The classic RepSeq / MIGEC case: a UMI followed by a known primer. The primer is what
        places the barcode, so the pattern can be searched for anywhere in the read.

        ```
        S1	aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
        ```

        ```bash
        migec checkout reads.fq.gz -b barcodes.txt -o co/
        ```

        The `NNNN t NNNN t NNNN` is one 12 nt UMI - the `t`s between the runs are matched pattern
        bases, not barcode. That distinction is the whole barcode space: 4^12, not 4^14.

        **Real example.** An HIV-1 Primer ID library, `SRR1763769`. `migec suggest` recovered the
        layout with nothing but the FASTQ:

        ```
        pattern  NNNNNNNNNcagtttaacttttgggccatcca
        ```

        and the assembled consensus places at HXB2 2,328-2,595 - the 3' end of protease into RT.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 2. Droplet single-cell (10x)

        A cell barcode and a UMI, back to back, at the very start of R1 - and **no constant
        sequence anywhere**. Nothing to search for, so the position is the placement. Three
        spellings, all the same run:

        ```bash
        migec checkout R1.fq.gz R2.fq.gz --preset 10x-v2 -o co/
        migec checkout R1.fq.gz R2.fq.gz --bc-pattern 'cell:0:16,16:26' -o co/
        migec checkout R1.fq.gz R2.fq.gz --bc-pattern '^XXXXXXXXXXXXXXXXNNNNNNNNNN' -o co/
        ```

        No `--max-offset`. A pattern with nothing to score is anchored at the first base for you,
        because a free scan over it has no evidence to choose an offset with - and refuses,
        correctly.

        Two things that catch people:

        - **A molecule is cell + UMI**, never the UMI alone. The same UMI in two cells is two
          molecules; migec keys on both.
        - **`refine` and `assemble` take R2.** R1 is 26 nt of barcode and nothing else, so after
          trimming it is empty. Run the later stages on the mate that carries cDNA.

        ```bash
        migec refine   co/PBMC_R2.fq.gz -o ref/ --expect-cells 1000
        migec assemble ref/PBMC.fq.gz   -o asm/
        ```

        **Real example.** `sc5p_v2_hs_PBMC_1k` VDJ-T: 3,155,166 reads, 100% assigned, 221,024
        barcodes at 14.28 reads each, **813 cells** called.
        """
    )
    return


@app.cell
def _(mo, show):
    mo.md(
        f"""
        ## 3. Targeted panel with fgbio read structures (TSO500)

        fgbio, Picard and samtools describe a layout as a *read structure*: runs of `M` (molecular
        barcode), `B` (sample barcode), `S` (skip) and `T` (template). migec takes them directly.

        | structure | translates to | platform |
        |---|---|---|
        | {show("5M5S+T")} | TSO500: 5 nt UMI, 5 nt spacer, then template |
        | {show("16B10M+T")} | 10x 5' |
        | {show("8M+T")} | a plain 8 nt inline UMI |

        TSO500 ctDNA is `5M5S+T +T` -- the UMI is on R1 and R2 is all template:

        ```bash
        migec checkout R1.fq.gz R2.fq.gz --preset tso500 -o co/
        migec checkout R1.fq.gz R2.fq.gz --read-structure 5M5S+T -o co/
        ```

        Warning: **5 nt is 1,024 barcodes, which does not identify a molecule on a real panel.**
        TSO500's own pipeline groups on the UMI *and the mapping position* (`fgbio
        GroupReadsByUmi`, after alignment). migec groups on the barcode, before any alignment
        exists, so it will report the space as saturated -- and on this chemistry that warning is
        the right answer, not a threshold to raise. Extract and tag here; group position-aware,
        downstream.

        When a chemistry really does put a UMI on **both** mates -- duplex sequencing, or MAGERI's
        dual-end design -- the two halves concatenate into one identifier and both must match:

        ```bash
        migec checkout R1.fq.gz R2.fq.gz --preset duplex -o co/
        ```

        Accepting only the first mate would emit 12 nt barcodes beside 24 nt ones, and every
        collision estimate downstream would be computed over two barcode spaces at once.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 4. Shallow libraries: 1-3 reads per molecule

        Not a special case - the common one. Bulk repertoire profiling and shallow 3' single-cell
        both look like this, and nothing about the commands changes. What changes is **what the
        numbers can mean**, and migec says so rather than reporting them flat:

        ```
            MIG size      groups    share
                   1      31,888    79.4%
                 2-3       8,176    20.4%

        warning: 79.4% of molecules were seen once. A consensus over one read is that read --
          the UMI is buying counting here, not error correction
        ```

        Three things stop applying below ~3 reads per molecule, and each is reported rather than
        silently degraded:

        - **Correction gets much harder.** At ~1 read per molecule about 80% of barcode errors
          have no observable parent - the true barcode was never sequenced - so they cannot be
          fixed by any method. migec corrects conservatively and never deletes a molecule.
        - **The split threshold is inert.** Telling two molecules apart inside one barcode needs
          ~30 reads; below that there is no evidence and migec does not invent any.
        - **Nothing is thresholded away.** `--min-reads` defaults to 1. A molecule seen once is
          still a molecule, and cutting it discards real sequence.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## 5. Bulk RNA-seq with an inline UMI

        A 10 nt UMI at the start of a single-end read, then the `GGG` a SMARTer template switch
        leaves behind:

        ```bash
        migec checkout reads.fq.gz --preset smarter-umi -o co/
        ```

        The `GGG` is scored, so it is a real anchor rather than decoration - if it is missing the
        read did not template-switch and the barcode is not where it is claimed to be.

        Warning: check the deposited FASTQ before assuming the UMI is in it. The pipeline this
        layout comes from (`ncgr/UMI-analysis`, SRP150352) moves the UMI **into the read header**,
        and SRA rewrites headers, so the archived copy has no UMI anywhere. `migec suggest` says so
        rather than guessing:

        ```
        warning: the only near-uniform run sits after the last constant sequence, with nothing to
          anchor it. That is what diverse payload looks like.
        ```

        ## 6. Duplex sequencing

        A 12 nt UMI and a 5 nt fixed spacer on **both** mates, 24 nt of identifier together:

        ```bash
        migec checkout R1.fq.gz R2.fq.gz --preset duplex -o co/
        ```

        Warning: this extracts the tags and emits **single-strand** consensuses. Pairing a
        molecule's two strands into a duplex consensus is not implemented, so no duplex error rate
        may be quoted from this output.

        ## Run one, end to end

        The fixtures below come from
        [`isalgo/umi_data`](https://huggingface.co/datasets/isalgo/umi_data). They are small on
        purpose - all the reads of a fraction of the *barcodes*, never a fraction of the reads,
        so the molecule size distribution is the real one.
        """
    )
    return


@app.cell
def _(Path):
    import tempfile

    from huggingface_hub import hf_hub_download

    def fixture(name):
        """One fixture from the public dataset. Cached by huggingface_hub after the first call."""
        return Path(
            hf_hub_download(repo_id="isalgo/umi_data", filename=f"ci/{name}", repo_type="dataset")
        )

    work = Path(tempfile.mkdtemp())
    return fixture, work


@app.cell
def _(fixture, mo, work):
    from migec.assemble import run as assemble_run
    from migec.refine import format_report, run as refine_run

    hiv = fixture("SRR1763769_umi0.5pct.fq.gz")
    hiv_ref = refine_run(hiv, work / "hiv_ref", sample_id="CTRL")
    hiv_asm = assemble_run(work / "hiv_ref" / "CTRL.fq.gz", work / "hiv_asm", sample_id="CTRL")

    mo.md(
        f"""
        ### HIV-1 Primer ID amplicon

        | | |
        |---|---|
        | reads | {hiv_ref['reads']:,} |
        | barcodes | {hiv_ref['barcodes']:,} |
        | molecules after correction | {hiv_ref['molecules']:,} |
        | consensuses | {hiv_asm['molecules']:,} |
        | mean emitted quality | Q{hiv_asm['mean_quality']:.1f} of a Q{hiv_asm['quality_cap']:.0f} cap |
        """
    )
    return assemble_run, format_report, refine_run


@app.cell
def _(assemble_run, fixture, mo, refine_run, work):
    tenx = fixture("sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz")
    tenx_ref = refine_run(tenx, work / "tenx_ref", sample_id="PBMC", expect_cells=10)
    tenx_asm = assemble_run(work / "tenx_ref" / "PBMC.fq.gz", work / "tenx_asm", sample_id="PBMC")

    mo.md(
        f"""
        ### 10x VDJ-T, 1% of cells

        | | |
        |---|---|
        | reads | {tenx_ref['reads']:,} |
        | cell barcodes seen | {tenx_ref['cells_observed']:,} |
        | molecules | {tenx_ref['molecules']:,} |
        | consensuses | {tenx_asm['molecules']:,} |
        | clonality | {tenx_ref['payload_clonality']:.4f} - a diverse repertoire, so payload agreement is worth ~{1 / max(tenx_ref['payload_clonality'], 1e-9):.0f}x odds |

        The cell barcode is carried through: `refine` keys on cell + UMI, and `<sample>.cells.tsv`
        has one row per cell. Cell *calling* needs the whole library rather than 1% of it, which is
        why `--expect-cells` is tiny here.
        """
    )
    return


@app.cell
def _(mo):
    mo.md(
        """
        ## Then what?

        The consensus is ordinary FASTQ, and the molecule identity is carried twice: in the read
        **name** (`<sample>.<cell>.<umi>`, for tools that drop FASTQ comments) and in tab-separated
        **SAM tags** (`RX QX CB CY BC MI cD`, for tools that keep them).

        ```bash
        minimap2 -ax sr -y ref.fa cons/S1.consensus.fq.gz | samtools sort -o S1.bam
        bwa mem -C     ref.fa cons/S1.consensus.fq.gz     | samtools sort -o S1.bam
        arda amplicon --r1 cons/S1.consensus.fq.gz -p S1  # sequence_id IS the molecule id
        salmon quant -i tx.idx -l A -r cons/S1.consensus.fq.gz -o quant/
        ```

        Warning: do **not** run alevin, bustools or STARsolo on a consensus FASTQ. They read the
        barcode out of a raw barcode read and deduplicate themselves - migec already did, and that
        read no longer exists. One consensus is one molecule, so a plain quantifier's count already
        is a molecule count. See `docs/downstream.rst`, where each of these was run.

        ## Where to look next

        - **[`docs/layouts.rst`](../docs/layouts.rst)** - every way to say where the barcode is
        - **[`docs/downstream.rst`](../docs/downstream.rst)** - what each aligner and quantifier did
        - **[`docs/checkout.rst`](../docs/checkout.rst)** - the pattern grammar, dual-end barcodes,
          and what the reported Phred is actually worth
        - **[`docs/refine.rst`](../docs/refine.rst)** - what evidence correction uses, and where it
          cannot work
        - **[`docs/assemble.rst`](../docs/assemble.rst)** - the consensus posterior and the quality
          cap
        - **`notebooks/barcode_space.py`** - is my barcode long enough?
        - **`notebooks/refine_diagnostics.py`** - the coverage curve, the barcode-rank plot, and
          where the errors are
        """
    )
    return


if __name__ == "__main__":
    app.run()
