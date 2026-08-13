suggest
=======

Where is the barcode? Read it off the data rather than off the protocol.

.. code-block:: bash

   migec suggest reads.fq.gz
   migec suggest reads.fq.gz -o out/ --cycles 40

A protocol description is written for the bench, not for the file: it says what was ordered, not
what arrived, and it does not know that the run was demultiplexed with the tag already trimmed, or
that the barcode is on the other mate.

How it works
------------

A UMI cycle is one the synthesiser **mixed**: all four bases near 1/4, about 2 bits. A constant
cycle is one base near 100%, about 0 bits. Everything else is payload. Segmenting the per-cycle
base composition on that gives a pattern that pastes straight into a barcode table.

.. code-block:: text

    cycle      A      C      G      T  1/4 dev     Q  layout
        0  0.271  0.205  0.257  0.267    0.045    33  N  UMI
        9  0.020  0.971  0.004  0.006    0.721    37  |  constant

   segments:
       0-8   umi         9 nt  (mean 1/4 deviation 0.038)
       9-31  constant   23 nt  (mean 1/4 deviation 0.718)  CAGTTTAACTTTTGGGCCATCCA

   pattern  NNNNNNNNNcagtttaacttttgggccatcca

That is a real HIV Primer ID library (``SRR1763769``) with nothing supplied but the FASTQ. Checking
that pattern out assigns 94.8% of reads, and :doc:`the consensus places at HXB2 2,328-2,595 <refine>`.

The measure is **total-variation distance from uniform**, not entropy: it is bounded (0 for a flat
cycle, 0.75 for a fixed one), scale-free, and does not need a log per cycle.

What it cannot do
-----------------

.. warning::

   **The pattern stops at the last constant run.** Composition alone cannot tell a UMI from diverse
   payload — both are four flat lines at 25%. What separates them is that a barcode is *anchored*
   and payload is not, so a uniform run with nothing constant after it is reported in the note and
   left out of the pattern. Claiming it would produce a pattern that matches everywhere.

This is the guard that makes ``suggest`` usable as an answer rather than a suggestion, and it fires
on real data: on the ctDNA runs of `Maruzani et al. 2024
<https://doi.org/10.1186/s12864-024-10737-w>`_ it correctly reports **no pattern** in either mate of
either accession — those submissions are aligned BAMs whose UMIs were generated in silico, so there
is no barcode in the reads to find.

If the barcode really is 3', raise ``--cycles`` so the primer past it is profiled too. For paired
data with no UMI in R1, try R2 — the note says so.

Options
-------

===================== ==========================================================================
``--cycles``          leading cycles profiled (default 60)
``--max-reads``       reads used; the composition converges long before the default 200,000
``--umi-deviation``   how far from flat a cycle may sit and still be UMI. Real synthesiser mixes
                      are routinely 20/30/30/20, which is 0.05
``-o/--out``          write ``suggest.cycles.tsv``, ``suggest.segments.tsv``, ``suggest.json``
===================== ==========================================================================
