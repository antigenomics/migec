How low can you go: exome, ctDNA and MRD
========================================

Three applications, one question. *What is the lowest allele frequency this library can detect?*
It is answered by two numbers, and the variant caller is neither of them:

.. code-block:: text

    N   molecules covering the site       -- what `migec assemble` counts
    p   per-MOLECULE error floor          -- the RT/first-cycle floor, `docs/quality_floor.rst`

``N`` sets how low you *could* go; ``p`` sets how low you *can*. Sequencing deeper raises reads per
molecule, not ``N``. Only more input DNA, or more tracked sites, raises the evidence.

.. code-block:: bash

    python scripts/detection_limit.py --input-ng 20 --sites 5              # a ctDNA panel
    python scripts/detection_limit.py --input-ng 50 --sites 30 --rt-error duplex   # MRD
    python scripts/detection_limit.py --from-json asm/assemble.json --sites 5

The three regimes
-----------------

Everything below is one of three situations, and knowing which you are in tells you what to buy.
The third is the one that surprises people, because the usual lever makes it worse.

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * -
     - molecule-limited
     - floor-limited
     - artifact-limited
   * - what binds
     - too few molecules cover the site
     - the chemistry errs at random, per molecule
     - the chemistry errs *systematically*, at this base
   * - symptom
     - the variant is absent from the library
     - the variant is as common as background
     - a reproducible false call at a fixed position
   * - fix
     - more input DNA, or more tracked sites
     - a lower floor: proofreading enzyme, or **duplex**
     - a **per-position background model** from known negatives
   * - what does **not** help
     - deeper sequencing, a better caller
     - deeper sequencing, more input, a better caller
     - **more molecules** -- see below

The molecule/floor crossover is at ``VAF = p/3``, the frequency at which a true variant molecule is
as rare as the chemistry's own false ones. At the default RT floor of 1e-4 that is **3.3e-5**, and
no amount of input DNA reaches below it.

Never: **an assay designed past its floor spends money on sequencing that cannot work.** A 50 ng,
30-site MRD panel has enough molecules for 6.9e-6 -- but on a single-strand protocol the floor sits
at 3.3e-5, five times higher. The molecules promise something the chemistry cannot deliver.

Never: **more molecules makes an artifact easier to call, not harder.** Measured on the
0%-certified arm at fixed 20 ng input, varying only sequencing depth:

.. list-table::
   :header-rows: 1
   :widths: 18 26 28 28

   * - depth
     - molecules per site
     - artifact called in
     - VAF reported
   * - 3.3x
     - ~8,000
     - 0 of 3 replicates
     - --
   * - 10x
     - ~12,000
     - **3 of 3 replicates**
     - **0.66%**

The bias is present at both depths. What changes is the statistical power to call it: a systematic
error does not average out with more molecules, so the extra evidence that makes a real variant
significant makes the artifact significant too. At 20 ng and 10x the 0.125% arm read 0.17%
(detected in 1 of 3) while the **true negative read 0.66% in 3 of 3** -- the negative outscoring
the positive. Below the artifact level, the ranking is meaningless without a background model.

ctDNA
-----

Molecule-limited almost always, because the input is a blood draw and cell-free DNA is scarce:
5-30 ng from 10 mL of plasma is typical, and 20 ng is only ~6,000 haploid genomes.

:doc:`variants` has the measurement: over 100 runs of cfDNA reference material at certified
frequencies, **the input mass decided the outcome and the caller did not**. At 0.125% VAF, 5 ng
gave a coin flip and 20 ng gave a certainty.

Two things a total molecule count hides, both measured on that panel by aligning to GRCh38:

* **Coverage is not uniform.** The weakest target held **0.31-0.61x** of the on-target mean, so a
  variant sitting on it has up to 3x fewer molecules than the average implies.
* **Off-target product is invisible without a reference, and it grows as input falls.** One locus
  outside any coding sequence took a share that tracked the DNA input almost perfectly:

  .. list-table::
     :header-rows: 1
     :widths: 30 30 40

     * - input
       - off-target share
       - on-target molecules
     * - 80 ng
       - 5-7%
       - 110,000-122,000
     * - 20 ng
       - 24%
       - 45,000-57,000
     * - 5 ng
       - **47-58%**
       - 11,500-12,700

  The *absolute* off-target count barely moves (~15,000 molecules at both 5 and 20 ng at the same
  depth) while on-target scales with input -- which is what template-independent product looks
  like. At 5 ng more than half the library is not evidence about anything.

So for ctDNA: count molecules **per target**, quote the weakest target rather than the mean, and
never read a library total as on-target depth. Precisely when input is scarce -- the case that
matters -- the total is most misleading.

How much plasma DNA do you need?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combining the measured per-target molecule counts with the arithmetic above answers the question
people actually ask, for this panel, at 95% detection and three supporting molecules:

.. list-table::
   :header-rows: 1
   :widths: 16 30 26 28

   * - input
     - molecules at the **weakest** target
     - limit of detection
     - detects 0.125%?
   * - 5 ng
     - 796
     - **0.79%**
     - no -- six times too high
   * - 20 ng
     - 3,529
     - **0.18%**
     - no -- marginally too high
   * - 80 ng
     - ~24,000
     - **0.026%**
     - yes, with room to spare

The reference material is certified at 0.125%, 0.25% and 1%, so **only the 80 ng arm can call the
lowest of those at every target in the panel**. Solving for the threshold: ~5,000 molecules at the
weakest target are needed for 0.125%, which on this panel means roughly **30 ng of input**.

Never: quoting the panel *average* would have said 20 ng was sufficient. It is not, for a variant
that happens to sit on the weakest amplicon -- and which amplicon a patient's variant sits on is
not something you get to choose.

What input actually buys: precision, not accuracy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Running the full chain -- ``assemble`` consensus, ``minimap2 -y``, LoFreq on the inferred panel --
over the undiluted arm recovers **PIK3CA H1047R** (``3:179234297 A>G``), and the certified 1%
dilution comes back at **0.92%** across three replicates. Across a 16x range of DNA input the point
estimate does not move; only its scatter does:

.. list-table::
   :header-rows: 1
   :widths: 12 8 16 14 16 18 16

   * - input
     - n
     - mean VAF
     - SD
     - CV observed
     - CV if Poisson
     - excess
   * - 5 ng
     - 8
     - 3.72%
     - 0.73 pp
     - 0.196
     - 0.119
     - **1.6x**
   * - 20 ng
     - 9
     - 3.54%
     - 0.38 pp
     - 0.108
     - 0.055
     - **2.0x**
   * - 80 ng
     - 6
     - 3.60%
     - 0.27 pp
     - 0.074
     - 0.033
     - **2.2x**

The assay is **unbiased** -- 3.5-3.7% at every input -- and precision improves roughly as
``1/sqrt(N)``. That is the practical meaning of a molecule count: more DNA does not give you a
different answer, it gives you a more certain one.

Never: **the observed scatter is consistently about twice the Poisson prediction.** Molecule
sampling explains only half of it; the rest is library preparation and PCR efficiency. So the limit
of detection computed by ``scripts/detection_limit.py`` is a *floor*, not a field estimate -- a
real assay will do worse, and validating against a dilution series is the only way to know by how
much.

The true negative is not empty
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The arm that matters most is the one certified at **0% mutant**, and running the same
consensus-then-LoFreq chain over it does not return nothing. It returns **9-11 calls per sample at
0.4-1.4% VAF**, and they are not noise:

* **94% of them are** ``-> G`` (14 A>G, 9 T>G, 6 C>G, against one C>A and one T>A).
* **Eight positions recur in 3 of 3 replicates**, with VAF reproducible to the third decimal --
  ``3:179234288 A>G`` at 0.0137 / 0.0133 / 0.0140, ``4:54733163 A>G`` at 0.0117 / 0.0127 / 0.0140.

A ``-> G`` bias is the signature of **2-colour chemistry**, where G is the base call for *no
signal*: any position that loses fluorescence reads as G. These runs are MiniSeq, which is
2-colour. Consensus does not remove it, because it is not a random sequencing error -- it is a
systematic, position-specific, sequence-context-driven bias that most reads of a molecule share.

Never: **the artifact lands on the hotspot too.** ``3:179234297 A>G`` is PIK3CA H1047R, and the
true-negative arm calls it at **0.58-0.79%** while the certified 1% arm reads 0.92%.

.. list-table::
   :header-rows: 1
   :widths: 46 27 27

   * - arm
     - measured VAF at H1047R
     - truth
   * - undiluted ``cell_line``
     - 3.6%
     - positive
   * - certified 1% (5 ng, 3.3x)
     - 0.92%
     - 1%
   * - **certified 0%** (20 ng, 10x)
     - **0.66%**
     - **0%**

Note: those two arms differ in input and depth as well as in truth, so this is not a matched
comparison -- but the artifact is at the same base, is reproducible across replicates, and is 72%
the size of the certified signal. A pipeline that reports 0.92% as a detection must explain why
0.66% is not one.

What follows for the caller
^^^^^^^^^^^^^^^^^^^^^^^^^^^

**A standard caller on consensus reads is not sufficient on its own.** The molecule count is
right, the consensus is right, and the caller still reports systematic false positives, because
nothing in that chain knows that *this base on this strand in this context* reads high.

What fixes it is a **per-position background model** built from samples known not to carry the
variant -- and that, rather than UMI handling as such, is what the UMI-aware and panel-of-normals
callers actually contribute:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - approach
     - what it models
   * - ``Mutect2`` + panel of normals
     - per-site artifact rate across a normal cohort
   * - ``Shearwater``, ``DREAMS-vc``
     - per-position beta-binomial / learned error model over a panel of normals
   * - ``UMIErrorCorrect``
     - per-position beta-binomial background
   * - ``UMI-VarCal``
     - per-position Poisson background
   * - ``LoFreq`` alone on consensus
     - the base qualities only -- which the artifact does not violate

So the recommendation in :doc:`variants` stands for *which* caller, and gains a condition: run it
against a background model. On this data a panel of normals is not optional, and the WT arm of a
reference material series is exactly the cohort to build one from.

.. code-block:: bash

    migec assemble rf/S1.fq.gz -o as/
    minimap2 -ax sr -y ref.fa as/S1.consensus.fq.gz | samtools sort -o S1.bam
    samtools index S1.bam
    # molecules at one target: one consensus record is one molecule, so count distinct MI
    samtools view S1.bam 17:7673727-7673832 | grep -o 'MI:Z:[^\t]*' | sort -u | wc -l

Exome
-----

Also molecule-limited, but for a different reason: an exome spreads its molecules over ~200,000
targets, so **the mean is meaningless and the distribution is everything**. A "100x mean" exome
routinely has thousands of targets in the single digits, and a variant in one of those is
undetectable no matter how good the caller is.

migec's contribution here is that its depth *is* a molecule count, so the per-target number you
compute is the real evidence rather than a duplicate-inflated read count. Capture panels also make
coordinate deduplication actively wrong -- probes pile reads on identical start positions whether
or not they came from one molecule, which is what ``notebooks/exome_capture.py`` demonstrates.

Target definitions do not need a vendor login: ``AstraZeneca-NGS/reference_data`` ships hg38 BEDs
for Agilent V2-V6, IDT V1, MedExome, NGv3 and a canonical CDS set (``SOURCES.md``).

Never: those BEDs are ``chr1``-style. An Ensembl reference is ``1``, and intersecting the two
returns **zero** rather than an error -- a silent clean-looking negative. Strip the prefix first.

MRD
---

The case migec was built for. Minimal residual disease means tracking a *known* variant set --
originally the leukaemic clone's **IGH rearrangement**, which is why the original MIGEC paper is a
repertoire paper -- at frequencies far below anything a discovery assay reaches.

Two things make that possible, and both are arithmetic rather than software:

**You know where to look.** No multiple-testing burden across the genome, so a single supporting
molecule can be meaningful where a discovery assay would need many.

**You track many sites at once.** Evidence pools. Thirty patient-specific variants is thirty times
the molecules that can carry a signal, and thirty times lower a reachable frequency. This is why an
MRD panel follows tens of mutations rather than one:

.. code-block:: text

    50 ng input, single site      LOD 2.1e-04
    50 ng input, 30 sites pooled  LOD 6.9e-06     <- 30x lower, same blood draw

**And then the floor stops you.** At the RT floor of 1e-4 the crossover is 3.3e-5, so the pooled
6.9e-6 above is unreachable on a single-strand protocol: 30 background false molecules against the
3 you are trying to call. Duplex sequencing -- requiring both strands of the original duplex to
agree -- moves the floor by orders of magnitude and makes the molecule count binding again.

.. list-table::
   :header-rows: 1
   :widths: 22 20 20 38

   * - protocol
     - floor ``p``
     - crossover
     - lowest useful VAF
   * - RT / cDNA (``rt``)
     - 1e-4
     - 3.3e-5
     - ~1e-4; below this, duplex or nothing
   * - ordinary polymerase (``medium``)
     - 1e-5
     - 3.3e-6
     - ~1e-5
   * - proofreading, no RT (``high``)
     - 1e-6
     - 3.3e-7
     - ~1e-6
   * - duplex (both strands agree)
     - ~1e-9
     - ~3e-10
     - molecule-limited, not floor-limited

Never: **migec v2 extracts duplex tags but emits single-strand consensuses.** It does not yet build
duplex consensus (``ROADMAP.md``), so no error-suppression claim here rests on duplex data. For IGH
MRD the clonotype half of the problem is `arda <https://github.com/antigenomics/arda>`_'s: migec
gives it one record per molecule, and its AIRR ``duplicate_count`` is then a molecule count, which
is the number a residual-disease burden should be computed from.

What to report
--------------

For any of the three, the honest read-out is three numbers, and migec prints all of them:

.. list-table::
   :header-rows: 1
   :widths: 26 30 44

   * - number
     - where
     - why it matters
   * - molecules **per target**
     - the consensus BAM, distinct ``MI``
     - the evidence; never the library total, never a read count
   * - barcode error, as a Phred
     - ``refine`` ``error_phred``
     - whether grouping is trustworthy at this depth (:doc:`umi_errors`)
   * - the emitted quality cap
     - ``assemble`` ``quality_cap``
     - the floor; anything below ``p/3`` is not measurable (:doc:`quality_floor`)

Never: do not quote a limit of detection without saying which of the two regimes produced it. "We
detect 0.1%" means one thing when 12,000 molecules cover the site and another when the chemistry
floor is 1e-3, and only the first is improved by a bigger blood draw.
