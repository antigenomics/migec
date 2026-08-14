Post-processing: what to run after migec, and what it is worth
==============================================================

``migec assemble`` writes ordinary FASTQ, one record per molecule. Everything after that is someone
else's tool, and this section is about which one, run how, and how much the collapsing was actually
worth -- measured against certified reference material rather than argued.

The headline
------------

On commercial cfDNA reference material with **certified** mutant allele frequencies, including a
0%-certified arm that is a true negative by construction, three replicates per arm, one panel, one
aligner, matched molecule-support thresholds:

.. list-table::
   :header-rows: 1
   :widths: 30 16 14 14 14 12

   * - pipeline
     - false calls per sample, 0% arm
     - 0.125%
     - 0.25%
     - 1%
     - measured VAF at 1%
   * - **migec + Mutect2**
     - **0.67**
     - 0 of 3
     - 1 of 3
     - **3 of 3**
     - 0.0103
   * - **migec + LoFreq**
     - **2.00**
     - **1 of 3**
     - **3 of 3**
     - **3 of 3**
     - **0.0102**
   * - UMIErrorCorrect (own consensus + own caller)
     - 7.67
     - 1 of 3
     - 3 of 3
     - 2 of 2
     - 0.0094

Substitutions only, because migec emits no indels by design and 56% of UMIErrorCorrect's calls are
deletions; both totals are in ``assets/ctdna_callers.tsv``. Detection is at the certified hotspot,
PIK3CA H1047R.

**migec + LoFreq matches the best sensitivity at every arm and reports 3.8x fewer false positives
on the true negative.** migec + Mutect2 gives the fewest false positives of anything measured here
and pays for it at 0.25%. Which of the two to run is a decision about which error is expensive,
and :doc:`Detection limits <detection>` is where that arithmetic lives.

Two things this does *not* say. It does not say Mutect2 is a better caller than LoFreq in general
-- both were run at defaults plus the one non-default each needs on consensus input, and that flag
matters more than the choice between them (see :doc:`Variant calling <variants>`). And it does not
say UMIErrorCorrect is a worse tool -- it is a different *pipeline*, aligning first and grouping on
*(position, UMI)*, and what the table measures is the end-to-end result of that order against this
one.

Read next
---------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - page
     - the question it answers
   * - :doc:`Downstream tools <downstream>`
     - which aligner, quantifier or QC tool takes the consensus, and what survives the trip
   * - :doc:`Variant calling <variants>`
     - which caller, which non-default flags it needs on consensus input, and which callers
       *replace* migec instead of following it
   * - :doc:`Detection limits <detection>`
     - how low the frequency can go before the answer is set by molecules, by the chemistry, or
       by an artifact -- and which of the three you are in
   * - :doc:`Grouping accuracy <grouping>`
     - map-first against collapse-first, scored against a known truth

.. toctree::
   :maxdepth: 1
   :hidden:

   downstream
   variants
   detection
