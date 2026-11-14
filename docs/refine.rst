refine -- correct the errors in the barcode
===========================================

Correct the barcode errors, then hand the reads on. This is the stage that decides **how many
molecules there were**.

.. code-block:: bash

   migec refine out/S1.fq.gz -o ref/
   migec refine out/S1.fq.gz -o ref/ --min-posterior 0.99      # correct less
   migec refine out/S1.fq.gz -o ref/ --no-payload --no-quality # what the count ratio alone does

Input is a per-sample FASTQ from :doc:`checkout` -- or any FASTQ, BAM, SAM or CRAM whose records
carry ``RX``, which is what an fgbio or Picard UMI pipeline writes
(:doc:`Bring your own UMI <byo_umi>`). Output is the same reads with the corrected barcode in
``RX`` and the original preserved in ``OX``, plus the barcode table.

Or the ``.mig`` buckets ``checkout --mig`` wrote:

.. code-block:: bash

   migec checkout reads.fq.gz -b barcodes.txt -o out --mig
   migec refine out/S1.000.mig -o ref/       # buckets in, buckets out
   migec assemble ref/S1.000.mig -o asm/     # ...and assemble skips its partition pass

The output is then buckets too, **re-partitioned on the corrected barcode**. A corrected barcode is
a different key and a key decides its bucket, so copying a bucket through unchanged would stop it
being a partition — and the reads whose barcode was corrected across a bucket boundary would be
grouped with strangers by the next stage. The audit trail moves with it: a ``.mig`` record has no
room for the pre-correction barcode the way a FASTQ comment has ``OX:Z:``, so
``<sample>.barcodes.tsv`` — every barcode with its parent — is the record of what was merged. It is
one row per barcode rather than two ``u64`` per read.

Both routes produce the same numbers, down to the estimated error rate and the consensus that comes
out the far end; ``tests/synthetic/test_mig_chain.py`` asserts exactly that.

What the evidence is
--------------------

A barcode within two substitutions of another is either an **error child** of it or an
**independent molecule** that drew a neighbouring barcode. Three things separate them, and only
the first survives on a shallow library.

**Counts.** A child is much smaller than its parent. This is the whole game on a deeply sequenced
amplicon and it is worth nothing at 1–3 reads per UMI: a parent with 2 reads and a child with 1 is
not an asymmetry, and two singletons are not one either. Two gates in the first version — a parent
must be strictly larger, a child at most half its parent — made a singleton-vs-singleton merge
impossible by construction.

**The barcode's own base quality**, at the position that differs. A sequencing miscall carries a
low Phred *exactly there*. For two mismatches, the likelihood multiplies both position-specific
error probabilities: one poor base cannot license another confidently observed mismatch. Missing
barcode quality falls back to the fitted library rate. A polymerase error from an early PCR cycle carries a high one in every
read that inherits it. ``checkout`` already writes this to ``QX`` and nothing read it. Works at one
read.

**Payload agreement.** A barcode error child is a read of the *parent's molecule*, so its payload
matches. An independent molecule at a nearby barcode has its own payload. Worth ``log(1/clonality)`` —
and the clonality is **measured**, by sampling random barcode pairs and asking how often two
unrelated barcodes carry the same sequence anyway:

.. code-block:: text

   clonality       0.0100 of random barcode pairs carry the same payload anyway
                   -- payload agreement is worth about 100x odds towards the same molecule here

In a diverse repertoire that is decisive. In a clonal library it is worth nothing, and the number
says so rather than the evidence being quietly over-trusted. Agreement also **lifts the count
gates**, which is what makes a singleton merge possible at all; disagreement **refuses** a merge
the count ratio would have made.

.. note::

   The error likelihood is a *rate*, not a conditional. It used to be a zero-truncated Poisson
   weighed against ``a_ind · p_size``, which is an expected count. The truncation divides out
   ``(1 − e^−λ)`` — precisely the term that says whether an error child should exist — so
   ``ZT-Poisson(1, λ) → 1`` for every small λ, and the error rate stopped mattering at exactly the
   coverage where nothing else was available.

How well it works, and where it cannot
--------------------------------------

What one run prints, on a simulated library of 20,000 molecules carrying an injected barcode error
rate of 3.0e-03:

.. code-block:: text

   barcodes    23,910 distinct
     merged    3,855 (16.1%) into a parent, 3,889 reads moved
   molecules   20,055 after correction          <- 20,000 were simulated

   barcode error   2.87e-03 per base            <- 3.0e-03 injected

``scripts/correction_accuracy.py`` scores it against the simulator's truth. Recall is reported
against the **reachable** set, because a child whose parent barcode was never sequenced has nothing
to merge into and correctly stays put:

.. list-table::
   :header-rows: 1

   * - reads/UMI
     - reachable
     - recall of those
     - precision
     - molecules kept
   * - 1.11
     - 0.204
     - 0.108
     - 0.818
     - **1.000**
   * - 2.32
     - 0.904
     - 0.816
     - 0.830
     - 0.987
   * - 3.12
     - 0.975
     - 0.914
     - 0.926
     - 0.991
   * - 7.12
     - 1.000
     - 0.979
     - 0.997
     - 0.999
   * - 13.30
     - 1.000
     - 0.983
     - 0.999
     - 1.000

.. warning::

   At ~1 read per UMI, **80% of barcode errors are unfixable in principle** — the parent barcode
   was never sequenced. Of the rest migec corrects about a tenth, and destroys no real molecule at
   any depth measured.

   That is the side to err on. A wrong merge deletes a molecule and nothing downstream can tell; a
   missed correction only inflates a count, and the count is reported next to the coverage
   histogram that explains it. ``--min-posterior`` is where you move along that trade.

Whitelists: ``--cell-whitelist``
--------------------------------

10x ships the list of cell barcodes it actually synthesised. An observed barcode one substitution
away from a list entry is usually a miscall of it — and *usually* is not *always*, which is the
whole difficulty.

.. code-block:: text

   whitelist   200 of 3,201 cell barcodes were on the list
               1 snapped to it (20 reads), 3,000 left off it
               3,000 are further than one substitution from every entry -- which is what
               measures the off-list prior, 1.59e-05 per barcode

**Never: The posterior needs a background hypothesis.** Without one the model asserts "the true barcode
is one of these 737,000", so every observed barcode is assigned to *some* entry with posterior 1.0
— an index-hopped read, an undeclared sample, free-floating ambient sequence, all silently absorbed
into whichever entry happens to be nearest. With it, "this barcode is not on the list and was read
correctly" competes, and for those it wins.

The prior on that background is **measured**: barcodes at distance ≥2 from every entry cannot be
single substitutions of anything on the list, so the share of reads they carry is a lower bound on
how much of the library is genuinely off-list.

.. warning::

   It is a prior on **this barcode**, not on the library. The whitelist prior is a share of
   ``1 − background`` spread over every entry, so a 737,000-entry list gives each candidate ~10⁻⁶.
   A background quoted as "1% of the library is off-list" would be four orders of magnitude larger
   than any candidate and would win every time. The comparable quantity is the off-list read share
   **divided by** the number of distinct off-list barcodes.

The posterior scales roughly as ``n_parent · e/3``, where ``e`` is the error at the base that
differs. So a snap needs **both** a well-used parent and a base the instrument is unsure of: at
Q30 nothing is overridden however popular the neighbour, and a barcode seen twice never wins an
argument. That is deliberate — a wrong snap moves a molecule into another cell and nothing
downstream can tell.

An ``N`` is **expanded, not discarded**: it is a base the instrument declined to call, consistent
with all four at ``e = 0.75``, so a barcode carrying one is still correctable rather than thrown
away with the molecule it tagged.

What correction left behind
---------------------------

The molecule count is only as good as the errors that were *not* corrected, so ``refine`` estimates
those directly rather than trusting the correction. A surviving barcode that still looks like a
child of a surviving neighbour is one the posterior declined to merge:

.. code-block:: text

   residual        1,294 molecules still look like children of a neighbour -- by count, or by
                   their reads agreeing on the molecule
                   5.25% of 1-read molecules; at >= 2 reads the rate is within the 5% target.
                   REPORTED, not applied -- a molecule seen three times is still a molecule

.. warning::

   **"A much larger neighbour" is not the test.** At 1–3 reads per UMI nothing is 20× anything, so
   a count-ratio criterion reports **zero residual in exactly the regime where the residual is
   worst** — the same trap the correction posterior itself fell into. The payload is what still
   separates them at one read: a neighbour whose reads agree on the molecule is a child whatever
   the counts say. Measured on the same library, the count-only estimator says 0 and the full one
   says 1,294.

Never: The MIG size threshold is **reported, never applied**. Every molecule is in the output whatever
it says: a molecule seen three times with no plausible parent is information, and cutting it
discards real sequence. ``--target-fdr`` sets which size the report points at; filtering on it is
a downstream decision, taken with the coverage histogram in view.

``<sample>.bins.tsv`` carries ``molecules``, ``suspected_residual`` and ``residual_fdr`` per
power-of-two size, next to the fraction that *was* merged.

Cell calling
------------

When the reads carry a cell barcode, ``refine`` counts **molecules per cell** — never reads, since
read depth is amplification — and calls cells with **OrdMag**, Cell Ranger's original rule: take
the 99th percentile of the top ``--expect-cells`` barcodes and keep everything within a tenth of
it.

.. code-block:: text

   cells       500 called of 20,500 barcodes seen, at >= 59 molecules (OrdMag)
               201,463 molecules in called cells (83.5% of all)
               the curve breaks at rank 310 (357 molecules) -- the knee, for comparison

The **knee** — the rank furthest from the chord joining the ends of the log-log curve — is reported
*next to* the call, not instead of it. OrdMag is a rule; the knee is what the data says on its own.
When they disagree by more than a factor of three the report says so, because one of them is
describing a library the other is not.

That distance **is** Kneedle's (`Satopaa 2011 <https://raghavan.usc.edu/papers/kneedle-simplex11.pdf>`_)
normalised ``|y_n - (1 - x_n)|`` up to a constant, taken at its **global maximum**.

.. warning::

   **Kneedle's published selection rule is degenerate without the smoothing spline the paper
   specifies.** It walks the *local* maxima and stops at the first that then falls by a
   sensitivity step. Measured on ``sc5p_v2_hs_PBMC_1k``'s 136,032 barcodes:

   =========================================  =======  ===========
   rule                                       rank     molecules
   =========================================  =======  ===========
   Kneedle, unsmoothed (287 local maxima)     15       1,057
   Kneedle, smoothed over 0.01 log-rank       358      314
   Kneedle, smoothed over 0.10 log-rank       180      435
   **global maximum (what migec computes)**   **376**  **308**
   =========================================  =======  ===========

   The unsmoothed answer would call fifteen cells. Properly smoothed it lands on the global
   maximum to within 5%, and over-smoothing drifts it off again — so the global maximum reaches
   the same place with no smoothing parameter to tune.

.. warning::

   **A global maximum always exists, so the knee is guarded.** An ambient-only library returns a
   rank too, and reporting it beside OrdMag would read as two methods disagreeing about where the
   cells are rather than as a curve with no cells in it. The knee is refused below **ten times the
   mean molecules per observed barcode** — one order of magnitude, the unit a log-log curve is
   read in — and the report then says the curve has no knee.

   The mean alone was tried and is too weak: an ambient-only library of 1-3 molecules per barcode
   is a step curve whose corner sits at 3 against a mean of 2.0 and clears a mean-only test at
   1.5x. On the real library the knee is 308 against a mean of 3.38, a factor of 91. The two cases
   are two orders apart, so the boundary is not a tuned one.

**Never: EmptyDrops-style rescue of low-count cells is deliberately not reproduced.** It is Cell
Ranger's job, and imitating it would make every comparison against their calls unreachable by
construction rather than by measurement. The benchmark gate is written against recall of their
cells, broken down by barcode rank — not a Jaccard we cannot reach.

``<sample>.cells.tsv`` has one row per cell barcode: ``cell``, ``molecules``, ``called``.

What it holds
-------------

The **barcode table**, never the reads: ``(key, count)`` plus this barcode's own evidence — the
mean error at each barcode position and a 32-base payload draft. The reads are streamed twice
instead, and the table size is reported for the same reason ``checkout`` reports its counters.

.. code-block:: text

   0.5 s, two passes over the reads
   peak RSS 184.3 MB of which the barcode table 2.2 MB

**And the table bounds itself.** Past 1 GB it range-partitions to disk, correction follows it into
the partition, and every table refine writes is streamed one bucket at a time. Nothing changes but
the wall clock: a partitioned run and a resident one agree on every number and on every output file
byte for byte, which is what ``tests/synthetic/test_refine_bucketed.py`` asserts. The budget is a
``refine.run()`` argument, never a CLI flag — a run that needs it is one whose library does not fit,
and the stage can see that for itself.

.. warning::

   **Correction is not complete within one plain range partition.** Two substitutions may
   affect two different prefixes. Correction uses three disjoint prefix partitions, rotating
   keys and quality positions together: at least one prefix remains unchanged. Overly wide
   legacy partitions are re-partitioned. Every pass proposes parents; duplicate proposals
   are reconciled using original-coordinate deterministic ties before any merge is applied.
   The distance-one excess estimator and distance-one residual diagnostic retain their
   separately labelled two-pass semantics; their radius is not the correction limit.

   The table budget bounds counter storage, not total process memory. Correction uses
   seqtree's native batch consumer and retains each barcode's best decision rather than all
   radius-two neighbours. For unique canonical DNA tags up to 32 nt, a query has at most 4,561
   radius-two hits: a 32-query primitive block is at most 3,502,848 bytes. Concurrent native
   hit vectors and blocks scale with the worker count, and the current bucket's index also
   occupies memory. The distance-one residual and whitelist paths still use ordinary native
   batches; the counter budget is not a hard limit on the whole process.

The two-substitution sequencing model uses the product of per-position error probabilities.
The multiple-polymerase-event component assumes independent events; the focused synthetic
regressions verify the implementation, not empirical calibration for every chemistry. Candidate
retrieval uses seqtree, including distance one, and acceptance retains the independent-molecule,
count, quality and payload hypotheses. See the seqtree technical appendix's
`UMI and barcode collapse discussion <https://antigenomics.github.io/seqtree/_static/evalue-appendix.pdf>`_.
Homopolymer indel support is an active follow-up and is not implemented by this substitution stage.

Output
------

============================== ====================================================
``<sample>.fq.gz``             the reads, ``RX`` corrected, ``OX`` = what it was
``<sample>.barcodes.tsv``      cell, umi, reads, corrected reads, parent
``<sample>.cells.tsv``         cell, molecules, called — only with cell barcodes
``<sample>.rank.tsv``          the barcode-rank curve and its CDF, log-spaced ranks
``<sample>.bins.tsv``          per MIG size: barcodes, reads, merged as error, entropy
``<sample>.sizes.tsv``         the MIG size spectrum at exact sizes: molecules and their reads
``<sample>.umi_errors.tsv``    per parent depth: error children, their reads, the rate implied
``refine.coverage.tsv``        molecules per power-of-two MIG size, after correction
``refine.json``                all of it, machine-readable
============================== ====================================================

``notebooks/refine_diagnostics.py`` draws all of them. Nothing is plotted inside the C++: a figure
has to be redrawable from a committed TSV long after the run.

``<sample>.umi_errors.tsv`` is the barcode error rate measured at every amplification depth
rather than once for the library, and it is what checks ``estimated_error``: a parent seen
:math:`c` times had :math:`cL` barcode bases to miscall, and its error children are what got
miscalled. Two estimators fall out of the same row and they fail in opposite directions, which is
the point of reporting both — :doc:`umi_errors` works through where each one holds and where it
stops. On a diverse library sequenced 25 deep it lands within 1% of a known injected rate.

The column worth reading first is ``fraction_erroneous`` in ``<sample>.bins.tsv``. About **94% of
singleton barcodes are error children** and ~0.2% at 2–3 reads. A flat curve, or one rising at high
counts, means correction is merging real molecules.

A merged read keeps its original barcode in ``OX:Z:``. A correction nobody can audit is a
correction nobody can check — and merges chain, so a read can end up two substitutions from where
it started while every individual step was one.
