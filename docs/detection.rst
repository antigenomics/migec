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

Never: **more molecules makes an existing artifact easier to call, not harder** -- but only where
the artifact is there to begin with. Two factors decide it, and they are separable. Splitting 72
runs by preparation and by molecule count:

.. list-table::
   :header-rows: 1
   :widths: 26 24 24 26

   * - preparation
     - molecules per site
     - calls per sample
     - n
   * - diluted, fewer molecules
     - 2,579
     - 1.0
     - 24
   * - **diluted, more molecules**
     - 7,354
     - **7.1**
     - 24
   * - undiluted, fewer molecules
     - 3,686
     - 1.1
     - 12
   * - undiluted, more molecules
     - 13,611
     - 1.4
     - 12

Within the **diluted** material, 2.9x the molecules gives **7.1x** the calls. Within **undiluted**
material at *more* molecules than that, 3.7x the molecules gives 1.3x. So the two factors do
different jobs:

* **Preparation decides whether the artifact exists.** The dilution series was made by mixing, and
  the extra handling is what the artifact tracks; the raw material barely shows it.
* **Molecule count decides whether you can see it.** A systematic error does not average out, so
  the evidence that makes a real variant significant makes a *present* artifact significant too.

That is why the usual lever backfires here. At 20 ng and 10x the 0.125% arm read 0.17% while the
**true negative read 0.66%** -- the negative outscoring the positive. Below the artifact level the
ranking carries no information without a background model.

ctDNA
-----

Molecule-limited almost always, because the input is a blood draw and cell-free DNA is scarce:
5-30 ng from 10 mL of plasma is typical, and 20 ng is only ~6,000 haploid genomes.

Note: **ctDNA and cfDNA are not synonyms and this page uses both deliberately.** Cell-free DNA is
the input, all of it, from a blood draw. Circulating tumour DNA is the tumour-derived fraction of
that input, and the VAF is what measures the fraction. An assay is run *on* cfDNA *for* ctDNA.

:doc:`variants` has the measurement over 100 runs of cfDNA reference material at certified
frequencies, and it is a two-part answer. **Above ~1% VAF the input mass decides the outcome and
the caller does not** -- 1% was called in 11 of 12 runs at an accurate 0.93%. **Below 1% neither
does**, because the assay becomes artifact-limited: 0.125% was called in 1 of 12 runs while the
0%-certified arm was called in 3 of 12.

Two things a total molecule count hides, both measured on that panel by aligning to GRCh38:

* **Coverage is not uniform.** Across 72 runs the weakest target held **0.09-0.64x** of the panel
  mean (median 0.36), so an average overstates the thinnest target by up to elevenfold.
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

  The *absolute* count barely moves (~15,000 molecules at both 5 and 20 ng at the same depth)
  while on-target scales with input, so at 5 ng more than half the library is not evidence about
  anything.

  Never: **it is adapter read-through, not a PCR product and not poly-G reads.** The sequences
  there are 85-95 bases soft-clipped with ~45 aligned (``87S52M``, ``94S45M``, ``95S44M``) at
  **MAPQ 4-16**, against MAPQ 60 for 1,844 of 1,907 molecules on real TP53, and they carry the
  TruSeq adapter ``GATCGGAAGAGCACACGTCTGAACTCCAGTCAC``. Short inserts let the read run past the
  fragment into the adapter; the leftover mismaps into a low-complexity locus which is **97% G with
  an 81 bp pure-G run**, against 19-33% G for every real amplicon.

  So it is removed by things that cost nothing:

  .. list-table::
     :header-rows: 1
     :widths: 26 22 52

     * - fix
       - where
       - effect
     * - ``-q 20`` MAPQ filter
       - alignment
       - removes the mismapped *alignments* -- MAPQ 4-16 against 60 -- and **does not reduce the
         call burden**. Measured below: it raises it
     * - adapter trimming
       - before ``checkout``
       - removes the cause. Diagnosed, still not measured
     * - ``--min-reads 3``
       - ``assemble``
       - a *different* artifact class, and the one that works; see below

  Note: the G-rich locus is a red herring worth naming, because it is the sort of thing that
  invites a poly-G explanation. The barcodes of those molecules have ordinary G content (mode 3 of
  12, same as a real amplicon) and under 5% of consensus records are even 50% G. The reference is
  G-rich; the reads are not.

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

Never: quoting the panel *average* would have said 20 ng was sufficient. It is not, for a variant
that happens to sit on the weakest amplicon -- and which amplicon a patient's variant sits on is
not something you get to choose. The weakest target holds as little as **0.09x** the panel mean, so
an average can overstate it elevenfold.

Never: **that table is the molecule-limited answer only, and it is optimistic.** It says 80 ng
reaches 0.026% and therefore calls 0.125% comfortably. Scoring actual calls against the certified
frequencies says otherwise -- 0.125% was detected in **1 of 12** runs, and the 0%-certified arm was
called in 3 of 12. The molecules are there; what is missing is a background model, because below
1% this panel is artifact-limited rather than molecule-limited. Use the table to rule inputs
*out*, never to rule one *in*.

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

``--min-reads``, and the artifact it removes
--------------------------------------------

Which floor applies to your chemistry, and what a consensus is worth there, is the other axis:
:doc:`assays` has the per-assay recipe and the RNA/DNA split of the pre-amplification floor. What
follows is the measurement that sets ``ctdna``'s ``--min-reads 3``.

Never: **``--min-reads`` defaults to 1, which is right for counting molecules and wrong for
calling variants.** A consensus over one read *is* that read -- no error correction at all, just
counting. Measured on certified cfDNA reference material at 20 ng and 10x:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - variant
     - ``--min-reads 1``
     - ``--min-reads 3``
     - ``--min-reads 5``
   * - ``3:179234297 A>G`` PIK3CA H1047R
     - 0.0118
     - 0.0118
     - 0.0125
   * - ``17:7674220 C>T``
     - 0.0139
     - 0.0137
     - 0.0138
   * - ``17:7673768 T>G``
     - 0.0041
     - **gone**
     - gone
   * - ``17:7674219 C>G``
     - 0.0072
     - **gone**
     - gone
   * - ``3:179199161 C>G``
     - 0.0067
     - **gone**
     - gone
   * - ``3:179234288 A>G``
     - 0.0077
     - **gone**
     - gone
   * - ``4:54733163 A>G``
     - 0.0068
     - **gone**
     - gone

Scored against truth across three replicates per arm at 20 ng / 10x, the effect is categorical:

.. list-table::
   :header-rows: 1
   :widths: 20 14 16 14 18 18

   * - arm
     - ``--min-reads``
     - calls/sample
     - ``-> G`` calls
     - H1047R called
     - mean VAF
   * - **0%, truth: absent**
     - 1
     - 10.0
     - **29**
     - **3/3 wrong**
     - 0.0066
   * - **0%**
     - **3**
     - 2.0
     - **0**
     - **0/3 correct**
     - --
   * - 0%
     - 5
     - 1.0
     - 0
     - 0/3 correct
     - --
   * - **1%, truth: present**
     - 1
     - 9.0
     - 17
     - 3/3 correct
     - 0.0101
   * - **1%**
     - **3**
     - 6.3
     - 3
     - **3/3 correct**
     - **0.0102**
   * - 1%
     - 5
     - 5.0
     - 3
     - 3/3 correct
     - 0.0104

**Specificity goes from 0% to 100% with no loss of sensitivity, and the measured frequency does not
move** (1.01% -> 1.02% against a certified 1%). Every ``-> G`` artifact in the true negative is
gone. Requiring three reads discards the molecules that
were never error-corrected, which is exactly the population the dark-G bias rides on.

It is not free, and the trade is worth stating in full. Measured on the 20 ng / 10x arm:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * -
     - ``--min-reads 1``
     - ``--min-reads 3``
     - change
   * - molecules at the site
     - 12,471
     - 7,393
     - **-41%**
   * - molecule-limited LOD
     - 5.1e-4
     - 8.5e-4
     - 1.7x worse
   * - artifact floor (measured)
     - **6.6e-3**
     - none observed
     - removed
   * - **what actually binds**
     - **6.6e-3** (the artifact)
     - **8.5e-4** (the molecules)
     - **8x better**

So the filter throws away 41% of molecules and makes the molecule-limited limit 1.7x worse -- while
removing a floor that was **13x higher than the molecule limit to begin with**. Never: judge
``--min-reads`` on molecules retained and it looks like a loss. Judge it on what binds and it is an
eightfold gain, because at ``--min-reads 1`` the molecules were never the constraint.

Note: retention across the arm was 63% at ``--min-reads 3`` and 54% at 5. Going to 5 costs another
9% of molecules and removed nothing further here, so 3 is where the curve flattens on this
chemistry. On a library with more reads per molecule the same threshold costs less; on a shallower
one it costs more, which is why it is a recommendation per assay rather than a new default.

The intermediate arms, and what the artifact was adding
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two arms between the extremes are where the artifact's *additivity* shows, because there the
true frequency and the artifact floor are the same order of magnitude. Three replicates per arm,
LoFreq on the same consensus BAM, ``assets/ctdna_callers.tsv``:

.. list-table::
   :header-rows: 1
   :widths: 16 14 18 18 18 16

   * - arm
     - certified
     - ``--min-reads 1``
     - ``--min-reads 3``
     - ``--min-reads 5``
     - H1047R seen
   * - 0.125%
     - 0.00125
     - 0.0052
     - **0.0014**
     - 0.0012
     - 1 of 3
   * - 0.25%
     - 0.0025
     - 0.0079
     - **0.0022**
     - 0.0023
     - 3 of 3
   * - 1%
     - 0.01
     - 0.0101
     - **0.0102**
     - 0.0104
     - 3 of 3

At ``--min-reads 1`` the measured frequency is the certified one **plus** a floor of 0.4-0.6%: the
0.125% arm reads 4.2x its truth and the 0.25% arm 3.2x, while the 1% arm -- where the floor is a
twentieth of the signal -- reads correctly. At ``--min-reads 3`` all three track the certified
value to within 12%. That is the same additivity the 0%-certified arm shows directly, seen from the
other side.

Never: **the 0.125% arm is detected in 1 replicate of 3 at every threshold.** That is a molecule
limit, not a filter effect, and no threshold moves it -- see the plasma-DNA table above.

The MAPQ floor is not the fix
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The adapter read-through above is real and a ``-q 20`` filter does remove those *alignments*. It
does not remove the calls. Applied to the same consensus BAMs before calling, three replicates per
arm:

.. list-table::
   :header-rows: 1
   :widths: 22 20 20 20 18

   * - arm
     - ``--min-reads``
     - calls/sample, MAPQ 0
     - calls/sample, MAPQ >= 20
     - change
   * - **0%, truth: absent**
     - 1
     - 10.0
     - **25.7**
     - **2.6x worse**
   * - 0%
     - 3
     - 2.0
     - 2.0
     - none
   * - 0%
     - 5
     - 1.0
     - 1.0
     - none
   * - 0.25%
     - 1
     - 10.3
     - 26.3
     - 2.5x worse
   * - 1%
     - 1
     - 9.0
     - 10.0
     - +1

The call set at MAPQ >= 20 is a strict **superset** of the one at MAPQ 0 on the 0% arm at
``--min-reads 1``: 13 positions shared, 18 added, none removed. The 18 additions sit at 0.22-0.38%
VAF and **16 of the 18 are** ``-> G`` -- the same dark-G artifact, previously just under the
threshold. Depth at those loci barely moved (6,434 against 6,473), so the filter did not take
evidence away; it moved where the caller's threshold fell.

Never: **a filter that removes the diagnosed cause is not the same thing as a filter that removes
the calls.** The mismapping is real, the MAPQ floor removes it, and the false-positive burden goes
*up*. At ``--min-reads 3`` and 5 the MAPQ floor changes nothing at all, because the artifact it
would have to remove is already gone. Adapter trimming -- removing the cause before ``checkout``
rather than its symptom after alignment -- is still diagnosed and still unmeasured.

What you can do blind
---------------------

That table is also a **test you can run without a panel of normals**, and it is the answer to
"what if I have no matched controls". migec knows something a caller does not: how many reads
built each molecule.

* A **real** variant at VAF ``f`` sits in ~``f`` of molecules regardless of how many reads built
  them, so its frequency does not move as the threshold rises.
* A **context artifact** is carried disproportionately by molecules made from few reads, because a
  singleton consensus is one raw read carrying the raw per-base error rate.

So call the same sample at several thresholds and keep what holds still. Note: **assemble does
not need rerunning** -- every consensus record already carries ``cD:i:N``, its true molecule depth,
and ``minimap2 -y`` carries it into the BAM. One assemble, one alignment, then filter on the tag:

.. code-block:: bash

   migec assemble rf/S1.fq.gz -o as/
   minimap2 -ax sr -y ref.fa as/S1.consensus.fq.gz | samtools sort -o S1.bam
   samtools index S1.bam

   for mr in 1 3 5; do
       samtools view -e "[cD]>=$mr" -b S1.bam > S1.mr$mr.bam && samtools index S1.mr$mr.bam
       lofreq call -f ref.fa -l panel.bed -o S1.mr$mr.vcf S1.mr$mr.bam
   done
   python scripts/blind_artifact_filter.py \
       --vcf 1=S1.mr1.vcf 3=S1.mr3.vcf 5=S1.mr5.vcf

That is a third of the work, and the subsets are guaranteed nested because they are the same
records. ``--min-reads`` on ``assemble`` is still right once you have *chosen* a threshold; this is
for when you want the **trend**, which is the thing that discriminates.

On the sample above it returns 4 real and 5 artifact, and **all five artifacts are** ``-> G`` --
which the script says out loud, because a spectrum that lopsided is a platform signature rather
than biology.

Note: this is weaker than a real per-position background model, and it is not a substitute for one
where normals exist. What it does is convert "I have no controls" from *nothing* into *one
orthogonal axis of evidence* -- and it costs two extra ``assemble`` runs, which are the cheapest
stage in the pipeline.

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
