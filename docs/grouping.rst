Grouping accuracy, and Calib
============================

"Which reads came from the same original molecule" is the question both migec and
`Calib <https://github.com/vpc-ccg/calib>`_ answer, so comparing them is a clustering comparison:
score each tool's partition of the reads against a known truth with the adjusted Rand index.

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
