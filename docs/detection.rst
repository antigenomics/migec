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

The two regimes
---------------

Everything below is one of two situations, and knowing which you are in tells you what to buy:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * -
     - molecule-limited
     - floor-limited
   * - what binds
     - too few molecules cover the site
     - the chemistry makes the same error
   * - symptom
     - the variant is absent from the library
     - the variant is indistinguishable from background
   * - fix
     - more input DNA, or track more sites
     - a lower floor: proofreading enzyme, or **duplex**
   * - what does **not** help
     - sequencing deeper, a better caller
     - sequencing deeper, more input DNA, a better caller

The crossover is at ``VAF = p/3`` -- the frequency at which a true variant molecule is as rare as
the chemistry's own false ones. At the default RT floor of 1e-4 that is **3.3e-5**, and no amount
of input DNA reaches below it.

Never: **an assay designed past its floor spends money on sequencing that cannot work.** A 50 ng,
30-site MRD panel has enough molecules for 6.9e-6 -- but on a single-strand protocol the floor sits
at 3.3e-5, five times higher. The molecules promise something the chemistry cannot deliver.

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
