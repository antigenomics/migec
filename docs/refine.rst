refine -- correct the errors in the barcode
===========================================

Correct the barcode errors, then hand the reads on. This is the stage that decides **how many
molecules there were**.

.. code-block:: bash

   migec refine out/S1.fq.gz -o ref/
   migec refine out/S1.fq.gz -o ref/ --min-posterior 0.99      # correct less
   migec refine out/S1.fq.gz -o ref/ --no-payload --no-quality # what the count ratio alone does

Input is a per-sample FASTQ from :doc:`checkout`; output is the same reads with the corrected
barcode in ``RX`` and the original preserved in ``OX``, plus the barcode table.

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

A barcode one substitution away from another is either an **error child** of it or an
**independent molecule** that drew a neighbouring barcode. Three things separate them, and only
the first survives on a shallow library.

**Counts.** A child is much smaller than its parent. This is the whole game on a deeply sequenced
amplicon and it is worth nothing at 1–3 reads per UMI: a parent with 2 reads and a child with 1 is
not an asymmetry, and two singletons are not one either. Two gates in the first version — a parent
must be strictly larger, a child at most half its parent — made a singleton-vs-singleton merge
impossible by construction.

**The barcode's own base quality**, at the position that differs. A sequencing miscall carries a
low Phred *exactly there*; a polymerase error from an early PCR cycle carries a high one in every
read that inherits it. ``checkout`` already writes this to ``QX`` and nothing read it. Works at one
read.

**Payload agreement.** A barcode error child is a read of the *parent's molecule*, so its payload
matches. An independent molecule at distance 1 has its own payload. Worth ``log(1/clonality)`` —
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

   **Correction is not bucketable by a plain range partition.** A range partition on the top *b*
   bits keeps a barcode and its 1-substitution neighbours together for every position except the
   top *b*/2 — and a neighbour that crosses a bucket boundary would never be found, so the
   partition would bound the memory and silently stop correcting. It runs in **two passes**: over
   the buckets as they stand, owning the positions the prefix does not touch, then over a copy with
   the key rotated past it, owning exactly the ones it hid.

   Two things follow, and both were bugs before they were rules. The table carries the
   **evidence** — the barcode's own quality and its payload draft — because a side array indexed
   against the entry list cannot survive a partition, and dropping it would leave the bucketed run
   correcting on the count ratio alone, which reports nothing at 1–3 reads per UMI. And the two
   passes **scan**; they do not merge. A barcode can have a plausible parent on each side of the
   boundary, and merging inside a pass takes the first candidate rather than the best: measured, 2
   barcodes in 6,591 landed on a different parent than the resident run gave them, and every table
   downstream moved with them.

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
