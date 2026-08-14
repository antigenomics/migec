Variant calling: which caller, and what it can possibly see
===========================================================

migec is not a variant caller and will not become one. This page answers the question people
actually arrive with -- *I want rare variants out of ctDNA or a tumour, what do I run after
assemble?* -- and the question underneath it, which is the one that decides the answer: **how many
molecules does the caller get, and is that enough for the frequency you are chasing?**

Two rules, then the numbers.

Compose or replace, never both
------------------------------

:doc:`downstream` draws the line for aligners and quantifiers: a tool that *transports* the
barcode composes with migec, a tool that *deduplicates* on it replaces a stage of migec. Variant
callers fall on both sides of that line, and running one from the wrong side is the most common
way to lose molecules silently.

.. list-table::
   :header-rows: 1
   :widths: 24 18 58

   * - caller
     - after ``assemble``?
     - why
   * - ``Mutect2``, ``LoFreq``, ``FreeBayes``, ``VarDict``, ``bcftools``
     - **yes**
     - they call from a BAM and never look at ``RX``. One consensus record is one molecule, so
       their depth *is* a molecule count and their allele fraction *is* a molecule fraction
   * - ``UMI-VarCal``
     - **no** -- it replaces ``assemble``
     - it does its own UMI pileup and its own consensus. Feed it raw reads with the UMI, or feed
       migec's consensus to a standard caller. Not both
   * - ``UMIErrorCorrect``
     - **no** -- it replaces ``checkout`` + ``refine`` + ``assemble``
     - it aligns first, then groups on *(position, UMI)* at edit distance <= 1 and consenses.
       It is an alternative pipeline, not a stage
   * - ``DREAMS-vc``, ``Shearwater``
     - **no** as normally run
     - both fit a per-position error model across a panel of normals built from *reads*. Run on a
       consensus the model is fitted to a different noise process than the one it will see

Running a UMI-aware caller on a consensus counts each molecule once and then collapses the result
again. Nothing errors; the molecule count just quietly drops.

.. warning::

   Never: do not derive an allele fraction from a consensus BAM's read depth *and* apply a
   UMI-aware caller's family-size filter. The family size is 1 by construction after ``assemble``
   -- every record already is a family -- so a ``--min-family-size 3`` filter discards the entire
   library and reports zero variants without an error.

What to run
-----------

.. warning::

   **A standard caller on consensus reads needs a background model, not just a threshold.** Run
   over the 0%-certified arm of a reference series, ``LoFreq`` on migec consensus returns 9-11
   calls per sample at 0.4-1.4% VAF; **94% are** ``-> G``, eight positions recur in 3 of 3
   replicates, and one of them is the PIK3CA H1047R hotspot itself. That is 2-colour chemistry's
   dark-G bias, and no consensus removes it because it is systematic rather than random. See
   :doc:`detection`. What removes it is a per-position background built from normals -- which is
   the real contribution of the UMI-aware and panel-of-normals callers below.

**If you have UMIs and you are running migec**, collapse first and use a standard caller:

.. code-block:: bash

    migec checkout reads.fq.gz --bc-pattern '0:12' --sample S1 -o co/
    migec refine   co/S1.fq.gz -o rf/
    migec assemble rf/S1.fq.gz -o as/
    minimap2 -ax sr -y ref.fa as/S1.consensus.fq.gz | samtools sort -o S1.bam
    samtools index S1.bam

    # then, in decreasing order of how much the choice matters:
    lofreq call   -f ref.fa -o S1.vcf S1.bam                    # balanced, and it is the
                                                                # specificity end of the trade
    gatk Mutect2  -R ref.fa -I S1.bam -O S1.vcf                 # more sensitive, needs
                                                                # FilterMutectCalls after it

The published comparison behind that ordering is Maruzani et al. 2024, which benchmarked six
callers on ctDNA at 0.5-7.5% VAF (:doc:`sources <reference>`; full citation below). Their finding,
in one line each:

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - caller
     - what the benchmark found
   * - ``LoFreq``
     - fewest putative false positives of any standard caller at every depth tested, and second
       only to Mutect2 on sensitivity. The balanced default
   * - ``Mutect2``
     - highest sensitivity; the most privately-called variants of any caller, which is the
       false-positive signature. Balanced once ``FilterMutectCalls`` runs -- which their
       comparison deliberately did not run
   * - ``UMI-VarCal``
     - fewest putative false positives of *all* callers on UMI-encoded data. Also the fewest
       calls overall -- 24 against 234-1,728 for the others on one sample
   * - ``UMIErrorCorrect``
     - most sensitive at the lowest VAF at every depth; false positives climb steeply with depth
   * - ``FreeBayes``
     - false-positive rate rises with allele frequency, unlike every other caller tested
   * - ``bcftools``
     - called nothing below ~8% VAF at any depth. Not a low-frequency caller and not written as one

.. note::

   That benchmark ran every caller at **default parameters and without base quality score
   recalibration**, and skipped ``FilterMutectCalls`` on purpose, to keep the comparison even. It
   is a fair ranking of defaults, not a ceiling for any single tool. Mutect2 in particular is
   being scored without the filter its own authors require.

**If you do not have UMIs**, none of this applies and no consensus is possible: use LoFreq for a
balanced call set or Mutect2 for sensitivity, and accept that the floor is the PCR error rate.

The number that decides it
--------------------------

Caller choice is a second-order effect. The first-order effect is whether the variant is present
in enough *molecules* to be called by anything, and that is set by the wet lab:

.. code-block:: text

    molecules at the site  =  input DNA / 3.3 pg  x  strands recovered  x  efficiency
    variant molecules      =  molecules at the site  x  VAF

A caller needs some minimum number of supporting molecules -- three is typical -- and the
supporting count is a *Poisson draw*, not a guarantee. If the expectation is 3, a third of
replicates see fewer than 3, and no threshold setting recovers a molecule that was never sampled.

This is why ``assemble`` reports molecules rather than reads, and why the RT floor
(:doc:`quality_floor`) is a per-*molecule* rate: an error made before amplification is in every
read of that molecule and no consensus removes it. At a floor of 1e-4, 10,000 molecules at a site
expect one false variant molecule from the chemistry alone -- and that, not the sequencer's Q30,
is what a caller's specificity is fighting.

Measured on cell-free DNA reference material
--------------------------------------------

``scripts/ctdna_titration.py`` runs the three migec stages over 100 runs of two SiMSen-Seq studies
on commercial cfDNA reference material with **certified** mutant allele frequencies
(:doc:`sources <reference>`; both are in ``SOURCES.md``):

* **PRJNA788522** -- 0% (``WT``), 0.125%, 0.25% and 1% VAF, crossed with 5/20/80 ng input and
  3.3/10/30x reads per UMI, three replicates each.
* **PRJNA507366** -- six polymerases on the same material, plus 0.031% and 0.0625% VAF.

The ``WT`` arms are true negatives: **33 runs and 6.1 million molecules in which the variant
frequency is zero by construction.** That is a false-positive floor measured on real chemistry
rather than on a simulation.

.. code-block:: bash

    python scripts/sra_fetch.py get SRR17220921 SRR17220924 ... -o simsen/
    python scripts/ctdna_titration.py --reads simsen/ --out ctdna/

.. csv-table:: migec over both titrations, averaged over replicates
   :file: ../assets/ctdna_titration.tsv
   :delim: tab
   :header-rows: 1

``molecules_per_amplicon`` is the evidence a caller gets at one site; ``variant_molecules`` is that
times the certified VAF; ``p_enough`` is the probability that at least three of them are really
there. ``assets/ctdna_titration_runs.tsv`` beside it has all 100 runs individually.

Where detection stops being about the caller
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Calling variants on the consensus and scoring against the certified frequencies gives the honest
performance of this panel and pipeline. The locus is PIK3CA H1047R (``3:179234297 A>G``), twelve
runs per arm:

.. list-table::
   :header-rows: 1
   :widths: 22 18 20 20 20

   * - certified VAF
     - detected
     - sensitivity
     - mean measured
     - verdict
   * - undiluted (~3.6%)
     - 24/24
     - 100%
     - 3.61%
     - reliable
   * - 1%
     - 11/12
     - **92%**
     - 0.93%
     - reliable, and accurate
   * - 0.25%
     - 4/12
     - **33%**
     - 0.23%
     - accurate when seen, mostly missed
   * - 0.125%
     - 1/12
     - **8%**
     - 0.04%
     - effectively undetectable
   * - **0% (true negative)**
     - **3/12**
     - -- (**25% false positive**)
     - 0.16%
     - **calls a variant that is not there**

Read the last two rows together. Below 1% VAF **sensitivity collapses faster than specificity
does**: at 0.125% the assay finds 8% of real variants while still calling 25% of the true
negatives. That is the worst possible shape, and it is not fixed by choosing a different caller --
:doc:`detection` shows the artifact is systematic, reproducible and ``-> G`` biased.

Never: **that table is ``--min-reads 1``, which is the wrong setting for variant calling**, and it
is why the numbers looked so poor. Re-run at ``--min-reads 3`` (:doc:`detection`), three replicates
per arm at 20 ng / 10x:

.. list-table::
   :header-rows: 1
   :widths: 16 20 20 20 24

   * - truth
     - ``mr=1`` called
     - ``mr=1`` VAF
     - ``mr=3`` called
     - ``mr=3`` VAF
   * - **0%**
     - 3/3 **wrong**
     - 0.66%
     - **0/3 correct**
     - --
   * - 0.125%
     - 1/3
     - 0.17%
     - 1/3
     - 0.05%
   * - **0.25%**
     - 3/3
     - **0.79%** (3.2x too high)
     - **3/3**
     - **0.22%**
   * - 1%
     - 3/3
     - 1.01%
     - 3/3
     - 1.02%

Note: **the artifact was inflating true positives, not only inventing false ones.** At 0.25% it
read 0.79%, and 0.25% + the 0.57% artifact floor is 0.82% -- the contamination is additive, so a
quantitative result was wrong by threefold at a frequency where the call itself looked fine.

**The reliable limit is 0.25%, quantified accurately**, once singleton molecules are excluded. That
is fourfold better than the same pipeline at ``--min-reads 1`` and in the range Illumina specify
for TruSight Oncology 500 ctDNA v2 (0.2% for SNVs at 20 ng). 0.125% remains out of reach here at
1 of 3 replicates -- that arm is molecule-limited, and the published <0.1% claims for this
chemistry assume both more input and a per-position background model.

Note: **depth does buy molecules, until it does not.** Deeper sequencing recovers more of the
molecules that are in the tube -- 20 ng of undiluted material gives 6,310 / 10,299 / 16,809
molecules per amplicon at 3.3 / 10 / 30x -- but the ceiling is the number of input molecules, and
past it further reads only raise reads-per-molecule. The honest summary is that depth and input are
*both* worth spending on until the molecule count stops rising.

Never: **the molecule total of a multiplex panel is not the count at a site**, and dividing by the
panel size is not enough either. Aligning to GRCh38 and counting per target shows the weakest one
holds **0.09-0.64x** the panel mean (median 0.36) -- so an average overstates the weakest target by
up to elevenfold. The panel is TP53 x2, PIK3CA x2 and KIT, inferred from coverage;
:doc:`detection` has the method and the off-target share that comes with it.

The barcode error tracks polymerase fidelity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

PRJNA507366 varies the polymerase while holding template and protocol fixed, which makes it an
independent check on ``refine``'s barcode-error estimate (:doc:`umi_errors`): a higher-fidelity
enzyme should miscall fewer barcode bases, and the estimator should say so without being told
which enzyme it is looking at.

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - polymerase
     - barcode Phred
     - error per base
     - runs
   * - Accuprime
     - 30.82
     - 8.3e-4
     - 3
   * - Accuprime HiFi
     - 30.84
     - 8.2e-4
     - 3
   * - Platinum
     - 31.01
     - 7.9e-4
     - 3
   * - Phusion
     - 31.40
     - 7.2e-4
     - 3
   * - Platinum HiFi
     - **33.99**
     - **4.0e-4**
     - 3
   * - Platinum SuperFi
     - **34.78-35.11**
     - **3.1-3.3e-4**
     - 7

The estimator separates the high-fidelity enzymes from the standard ones by **3 Phred, a factor of
two in error rate**, and puts Platinum HiFi and Platinum SuperFi together at the top -- from the
barcodes alone, with no reference, no alignment and no knowledge of the experiment. Note: these
are *within-study* comparisons; the 80 ng paired-enzyme arm used a different library prep and its
absolute numbers are not comparable with the rows above.

Why this dataset and not the published ctDNA benchmark
------------------------------------------------------

Maruzani et al. benchmarked on ``SRR10296599`` and eight metastatic breast cancer runs, and had to
**generate the UMIs in silico** -- 9 nt, Phred fixed at 37, assigned by Poisson to reads sharing
start and end positions -- because neither deposit kept its real ones. We confirmed that
independently: both runs report ``nreads=2`` with no index read, ``migec suggest`` finds no
barcode in either mate, and ``vdb-dump`` shows an empty linkage group. ``scripts/sra_fetch.py
probe`` is the one-line version of that check.

Two consequences for how far their UMI results carry:

* **Their synthetic UMIs cannot be wrong.** A random 12-mer with a fixed Phred has no error rate,
  no collisions beyond chance, and no ambiguous parent. On the real reference material above,
  migec measures the barcode error at Q30 -- so roughly 1.3% of 12 nt barcodes carry at least one
  miscalled base, and correcting them is most of what separates one UMI-aware caller from another.
  The benchmark removed the variable it was measuring.
* **Their assignment rule assumes reads of a molecule are co-terminal**, which
  :doc:`X1 measured as false <fragmented>` (7.8% of 10x groups overall). cfDNA has preferred cut
  sites so it is less wrong for a capture panel than for 3' GEX, but it is not free.

This is not a criticism of the ranking -- it is the best available comparison of caller *defaults*,
and we quote it above. It is a statement about which half of it transfers: the standard-caller
ordering rests on real reads and carries over; the UMI-aware arm rests on UMIs that were invented
after the fact.

.. note::

   Public ctDNA data with real, recoverable UMIs does exist -- it just was not the data that
   benchmark used. ``PRJNA507366`` and ``PRJNA788522`` both carry a 12 nt inline UMI that survived
   deposition, and ``migec suggest`` recovers it from base composition alone with no prior
   knowledge of the protocol. ``SOURCES.md`` records both.

Citations
---------

* Maruzani R, Brierley L, Jorgensen A, Fowler A. *Benchmarking UMI-aware and standard variant
  callers for low frequency ctDNA variant detection.* BMC Genomics 2024;25(1):827.
  `doi:10.1186/s12864-024-10737-w <https://doi.org/10.1186/s12864-024-10737-w>`_, PMID 39227777.
* Sater V, Viailly PJ, Lecroq T, Prieur-Gaston E, Bohers E, Viennot M, Ruminy P, Dauchel H, Vera P,
  Jardin F. *UMI-VarCal: a new UMI-based variant caller that efficiently improves low-frequency
  variant detection in paired-end sequencing NGS libraries.* Bioinformatics 2020;36(9):2718-2724.
  `doi:10.1093/bioinformatics/btaa053 <https://doi.org/10.1093/bioinformatics/btaa053>`_,
  PMID 31985795.
* Sater V, Viailly PJ, Lecroq T, Ruminy P, Berard C, Prieur-Gaston E, Jardin F. *UMI-Gen: A
  UMI-based read simulator for variant calling evaluation in paired-end sequencing NGS libraries.*
  Comput Struct Biotechnol J 2020;18:2270-2280.
  `doi:10.1016/j.csbj.2020.08.011 <https://doi.org/10.1016/j.csbj.2020.08.011>`_, PMID 32952940.
* Osterlund T, Filges S, Johansson G, Stahlberg A. *UMIErrorCorrect and UMIAnalyzer: Software for
  Consensus Read Generation, Error Correction, and Visualization Using Unique Molecular
  Identifiers.* Clin Chem 2022;68(11):1425-1435.
  `doi:10.1093/clinchem/hvac136 <https://doi.org/10.1093/clinchem/hvac136>`_, PMID 36031761.
  The source of ``PRJNA788522``, the titration above.
* Filges S, Yamada E, Stahlberg A, Godfrey TE. *Impact of Polymerase Fidelity on Background Error
  Rates in Next-Generation Sequencing with Unique Molecular Identifiers/Barcodes.* Sci Rep
  2019;9(1):3503. `doi:10.1038/s41598-019-39762-6 <https://doi.org/10.1038/s41598-019-39762-6>`_,
  PMID 30837525. The source of ``PRJNA507366``.
* `BesenbacherLab/UMIseq_variant_calling <https://github.com/BesenbacherLab/UMIseq_variant_calling>`_
  -- a gwf workflow comparing Shearwater, Mutect2, VarScan2 and DREAMS-vc on cfDNA, with a panel of
  normals built first. Read as the design reference for the panel-of-normals approach; not run here.
