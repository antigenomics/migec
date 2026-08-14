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
   :widths: 26 14 11 11 11 13 14

   * - pipeline
     - false calls per sample, 0% arm
     - 0.125%
     - 0.25%
     - 1%
     - measured VAF at 1%
     - median depth
   * - **migec + Mutect2**
     - **0.67**
     - 0 of 3
     - 1 of 3
     - **3 of 3**
     - 0.0103
     - 2,811 molecules
   * - **migec + LoFreq**
     - **2.00**
     - **1 of 3**
     - **3 of 3**
     - **3 of 3**
     - **0.0102**
     - 2,832 molecules
   * - no consensus + LoFreq
     - 5.67
     - 0 of 3
     - 3 of 3
     - 3 of 3
     - 0.0127
     - **52,628 reads**
   * - UMIErrorCorrect (own consensus + own caller)
     - 7.67
     - 1 of 3
     - 3 of 3
     - 2 of 2
     - 0.0094
     - 5,010

Substitutions only, because migec emits no indels by design and 56% of UMIErrorCorrect's calls are
deletions; both totals are in ``assets/ctdna_callers.tsv``. Detection is at the certified hotspot,
PIK3CA H1047R.

**migec + LoFreq matches the best sensitivity at every arm and reports 3.8x fewer false positives
on the true negative.** migec + Mutect2 gives the fewest false positives of anything measured here
and pays for it at 0.25%. Which of the two to run is a decision about which error is expensive,
and :doc:`Detection limits <detection>` is where that arithmetic lives.

Note: **UMIErrorCorrect is not losing on evidence.** The median depth at the sites it called is
7,475 molecules against migec + LoFreq's 4,046 at the matched threshold -- it has *more* molecules
and still reports 3.8x more calls on a sample where the right answer is none. The difference is in
what each pipeline does with them, not in how many it has.

What the collapsing itself is worth
-----------------------------------

The **no consensus** row is the same reads, the same trimming, the same barcode correction, the
same aligner and the same caller as the migec + LoFreq row. One thing differs: a record is a read
rather than a molecule. It is the row that says what ``assemble`` is for, and it moves everything:

* **2.8x the false positives** on the certified true negative -- 5.67 calls per sample against
  2.00.
* **The measured frequency stops being right.** 0.0127 against a certified 1% (1.27x) and 0.0038
  against a certified 0.25% (1.52x), where the consensus reads 0.0102 and 0.0022. The gap is an
  additive floor of 0.1-0.3%, which is what an uncollapsed pileup carries and what a consensus
  removes.
* **It detects less, from 38x more depth.** At 0.125% the consensus finds the hotspot in 1 of 3
  replicates and the read pileup in **0 of 3**, with 197,772x read coverage against 5,903
  molecules.

Never: **a read count is not a molecule count, and no amount of the first substitutes for the
second.** 153,675x coverage on a sample whose right answer is "nothing" produced more wrong answers
than 4,046 molecules did. Depth buys statistical power over a noise process the depth itself does
not reduce; collapsing reduces the noise process.

Note: the surviving artifact is also *more* concentrated after collapsing -- 77% of migec's
0%-arm substitutions are ``-> G`` against 41% without the consensus. That is the expected
direction: the consensus removes the random per-read errors and leaves the common-mode ones, which
is exactly the class :doc:`detection` says needs a per-position background model rather than a
threshold.

Never: **UMI-VarCal is missing from that table because it could not be run, not because it lost.**
It requires paired-end reads -- ``Extract.py`` pairs reads *by read id*, holding each one until the
same id appears a second time -- and these arms are single-end 151 nt, so both of its intermediate
FASTQs come out empty and the extracted BAM holds 0 of 901,938 reads after its log has reported
"Working 100 %". That is an input-shape mismatch, not a configuration error, and it is recorded
with the evidence in ``SOURCES.md`` rather than as a poor score.

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
