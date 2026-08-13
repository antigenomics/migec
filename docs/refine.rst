refine
======

Correct the barcode errors, then hand the reads on. This is the stage that decides **how many
molecules there were**.

.. code-block:: bash

   migec refine out/S1.fq.gz -o ref/
   migec refine out/S1.fq.gz -o ref/ --min-posterior 0.99      # correct less
   migec refine out/S1.fq.gz -o ref/ --no-payload --no-quality # what the count ratio alone does

Input is a per-sample FASTQ from :doc:`checkout`; output is the same reads with the corrected
barcode in ``RX`` and the original preserved in ``OX``, plus the barcode table.

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

⛔ **EmptyDrops-style rescue of low-count cells is deliberately not reproduced.** It is Cell
Ranger's job, and imitating it would make every comparison against their calls unreachable by
construction rather than by measurement. The benchmark gate is written against recall of their
cells, broken down by barcode rank — not a Jaccard we cannot reach.

``<sample>.cells.tsv`` has one row per cell barcode: ``cell``, ``molecules``, ``called``.

What it holds
-------------

The **barcode table**, never the reads: ``(key, count)`` plus this barcode's own evidence — the
mean error at each barcode position and a 32-base payload draft. The reads are streamed three
times instead, and the table size is reported for the same reason ``checkout`` reports its
counters.

.. code-block:: text

   0.5 s, three passes over the reads
   peak RSS 184.3 MB of which the barcode table 2.2 MB

.. warning::

   **Correction is not bucketable by a plain range partition.** A range partition on the top *b*
   bits keeps a barcode and its 1-substitution neighbours together for every position except the
   top *b*/2 — and a neighbour that crosses a bucket boundary can never be found. Doing this in
   buckets needs two passes with the key rotated, so that every pair shares a bucket in at least
   one of them. Until that lands the table is held whole.

Output
------

============================== ====================================================
``<sample>.fq.gz``             the reads, ``RX`` corrected, ``OX`` = what it was
``<sample>.barcodes.tsv``      cell, umi, reads, corrected reads, parent
``<sample>.cells.tsv``         cell, molecules, called — only with cell barcodes
``<sample>.rank.tsv``          the barcode-rank curve and its CDF, log-spaced ranks
``<sample>.bins.tsv``          per MIG size: barcodes, reads, merged as error, entropy
``refine.coverage.tsv``        molecules per power-of-two MIG size, after correction
``refine.json``                all of it, machine-readable
============================== ====================================================

``notebooks/refine_diagnostics.py`` draws all of them. Nothing is plotted inside the C++: a figure
has to be redrawable from a committed TSV long after the run.

The column worth reading first is ``fraction_erroneous`` in ``<sample>.bins.tsv``. About **94% of
singleton barcodes are error children** and ~0.2% at 2–3 reads. A flat curve, or one rising at high
counts, means correction is merging real molecules.

A merged read keeps its original barcode in ``OX:Z:``. A correction nobody can audit is a
correction nobody can check — and merges chain, so a read can end up two substitutions from where
it started while every individual step was one.
