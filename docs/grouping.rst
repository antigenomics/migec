Grouping accuracy: Calib, UMI-tools, fgbio
==========================================

"Which reads came from the same original molecule" is the question migec,
`Calib <https://github.com/vpc-ccg/calib>`_, `UMI-tools <https://github.com/CGATOxford/UMI-tools>`_
and `fgbio <https://fulcrumgenomics.github.io/fgbio/>`_ all answer, so comparing them is a
clustering comparison: score each tool's partition of the reads against a known truth with the
adjusted Rand index.

.. code-block:: bash

   python scripts/compare_calib.py --truth truth_reads.tsv \
       --migec out/S1.fq.gz --calib calib_out.cluster

``truth_reads.tsv`` is ``read_id``/``molecule_id``; ``tests/synthetic/_sim.py`` writes one. Any
tool can be scored by handing it a two-column partition through ``--partition name=file.tsv``.

ARI alone is not enough
-----------------------

A single number hides the direction of the error, and the two directions have opposite costs:

* **Splitting** one molecule across several clusters inflates the molecule count. Recoverable — the
  reads are still there and still correct.
* **Merging** several molecules into one cluster mixes their sequences, which is what destroys a
  real variant. Not recoverable by anything downstream.

So the script reports the fraction of reads in split molecules and the fraction in merged clusters
alongside the ARI, and the tests assert on the *direction*, not just the magnitude.

Where each tool wins
--------------------

Calib clusters on the barcode **and** the read sequence, with a locality-sensitive index over
minimizers. migec today groups on the barcode alone — ``assemble`` is what will split a group by
its sequence, and it lands in M1. That difference is exactly measurable:

.. list-table::
   :header-rows: 1
   :widths: 14 14 14 18 20 20

   * - UMI length
     - UMI error
     - ARI
     - reads split
     - reads merged
     - clusters / molecules
   * - 12 nt
     - 0
     - **1.0000**
     - 0.0000
     - 0.0000
     - 2000 / 2000
   * - 12 nt
     - 5·10⁻³
     - 0.9348
     - 0.5165
     - 0.0004
     - 2928 / 2000
   * - 8 nt
     - 0
     - 0.9917
     - 0.0000
     - 0.0267
     - 1974 / 2000
   * - 6 nt
     - 0
     - 0.8877
     - 0.0000
     - 0.3982
     - 1575 / 2000

2000 molecules, 8 reads each, simulated. Read the rows as three separate statements:

* **A clean 12 nt barcode needs nothing cleverer.** 4¹² is 16.8 million; collisions are negligible
  and barcode-only grouping is exact. Calib cannot beat 1.0.
* **UMI errors split, and only split.** Every read whose barcode picked up a substitution starts a
  cluster of its own — 52% of reads at a 5·10⁻³ per-base rate — while merging stays at 0.04%. This
  is what ``migec refine`` corrects (M3) and what Calib avoids by clustering barcodes at an edit
  distance.
* **A short barcode merges, and no amount of barcode cleverness fixes it.** At 6 nt the birthday
  bound guarantees collisions: 40% of reads land in a cluster holding more than one molecule. Two
  molecules that drew the same barcode are separable only by their *sequence*, which is precisely
  what Calib uses and what migec will use at ``assemble``.

The point of tabulating it is that the gap has a known size and a known cause. It is the collision
rate, and ``effective_length`` in ``checkout.summary.tsv`` predicts it before any clustering runs.

.. note::

   ``tests/synthetic/test_grouping_accuracy.py`` asserts the migec column on every test run, so the
   number does not rot between the occasions when someone has Calib installed. The Calib column
   needs Calib; see ``SOURCES.md`` for how to get it.

UMI-tools and fgbio: what the alignment is worth
------------------------------------------------

`UMI-tools <https://github.com/CGATOxford/UMI-tools>`_ ``group`` and
`fgbio <https://fulcrumgenomics.github.io/fgbio/>`_ ``GroupReadsByUmi`` are the map-first tools:
they align the raw reads and group on *(position, UMI)*. migec groups on *(sample, cell, UMI)* and
aligns once, afterwards. That is one difference and it has one consequence, which is what this
comparison measures — **the position is only evidence when reads land in different places.**

.. code-block:: bash

   python scripts/compare_grouping.py --out /tmp/cmp --molecules 20000 --clones 200 --coverage 5

The table is ``assets/grouping_tools.tsv``; 20,000 molecules, a 12 nt barcode at a 3·10⁻³ per-base
error rate, on one laptop core. ``clones`` is how many distinct sequences the molecules were drawn
from, so it is exactly how much the aligner has to work with.

.. list-table:: reads per molecule 5, varying how much the reference tells the aligner
   :header-rows: 1
   :widths: 12 16 12 14 14 12 12

   * - clones
     - tool
     - ARI
     - reads split
     - reads merged
     - seconds
     - peak RSS
   * - 1
     - **migec**
     - **0.9967**
     - 0.0111
     - **0.0065**
     - **0.17**
     - 234 MB
   * - 1
     - UMI-tools
     - 0.9864
     - **0.0056**
     - 0.0298
     - 4.61
     - **233 MB**
   * - 1
     - fgbio
     - 0.9817
     - **0.0056**
     - 0.0389
     - 7.98
     - 588 MB
   * - 200
     - migec
     - 0.9985
     - 0.0038
     - 0.0016
     - **0.18**
     - 231 MB
   * - 200
     - **UMI-tools**
     - **0.9994**
     - **0.0034**
     - **0.0000**
     - 1.49
     - **132 MB**
   * - 200
     - fgbio
     - **0.9994**
     - **0.0034**
     - **0.0000**
     - 4.65
     - 514 MB
   * - 20000
     - migec
     - 0.9987
     - 0.0034
     - 0.0015
     - **0.19**
     - 229 MB
   * - 20000
     - **UMI-tools**
     - **0.9995**
     - 0.0029
     - **0.0000**
     - 1.82
     - **136 MB**
   * - 20000
     - fgbio
     - **0.9995**
     - **0.0028**
     - **0.0000**
     - 4.59
     - 570 MB

Three statements, in the order they matter:

* **On one reference, migec wins, and it wins on the direction that cannot be undone.** A single
  amplicon, a clonal control, a targeted ctDNA panel: every read maps to the same place, the
  position carries nothing, and the map-first tools are left grouping on the barcode alone with no
  error model for it. They put **3.0%** (UMI-tools) and **3.9%** (fgbio) of reads into clusters that
  mix molecules, against migec's **0.65%** — 4.6× and 6× fewer molecules destroyed. migec pays for
  it in splitting (1.1% against 0.56%), which inflates a count and is recoverable.
* **On a diverse reference, the map-first tools win, by 0.001 ARI.** With 200 or 20,000 distinct
  sequences the aligner separates molecules that collided on a barcode, which barcode-only grouping
  cannot do by construction. The gap is the collision rate and nothing else, and it is small
  because a 12 nt barcode is large. ``assemble``'s linkage sub-clustering recovers that
  discriminating power from the payload with no aligner, but it does it after grouping, so it is
  not in this column.
* **migec is 8–48× faster and does not need the alignment at all.** 0.17–0.26 s against 0.98–4.61 s
  (UMI-tools) and 3.90–7.98 s (fgbio), *including* the aligner run the other two cannot skip.
  fgbio's memory is a JVM heap; UMI-tools streams a BAM and stays flat, while migec's grows with
  the barcode count until the table partitions itself.

.. note::

   Depth changes nothing about the ranking: over 1.2, 2.5, 5 and 10 reads per molecule the ARI gap
   holds at ~0.001 and the speed ratio at 8–28×. Those rows are in the same TSV.

.. note::

   The dividing line for a downstream tool is **transport vs deduplicate** — a tool that carries
   ``RX`` composes with migec, a tool that dedups on it replaces a stage of it. UMI-tools and fgbio
   are the second kind, which is why they are compared here rather than in
   :doc:`downstream`.
