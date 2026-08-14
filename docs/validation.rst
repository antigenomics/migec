Validation — the spike-in metric
================================

The question a UMI pipeline has to answer is not "can you remove errors". Anything removes errors
by discarding rare sequences. The question is:

   **can you remove PCR and sequencing error while keeping a real variant that is rarer than the
   error cloud around it?**

Shugay *et al.* 2014 built the instrument that answers it. Three known IGH clones were spiked into
a real B-cell library — ``EHEB`` and two deliberate variants at one and two substitutions from it:

============  ==========================================  =========  ===========
clone         junction (48 nt)                            subs       % plasmid
============  ==========================================  =========  ===========
EHEB          ``TGTGCGAGAGATGATGGCGGGGG…GACTTTGG``        0          ~99
EHEB-V1       ``TGTG``\ **G**\ ``GAGAGATGATGGCGGGGG…``    1          ~0.5
EHEB-V2       ``TGTGCGA``\ **CA**\ ``CATGATGGCGGGGG…``    2          ~0.05
============  ==========================================  =========  ===========

The metric is the ratio of a real spike-in to the **worst error at the same substitution
distance**:

* ``Err1`` — the most abundant junction exactly 1 substitution from EHEB, *excluding* V1
* ``Err2`` — the most abundant junction exactly 2 substitutions from EHEB, *excluding* V2

.. list-table::
   :header-rows: 1
   :widths: 25 20 27 28

   * - quantity
     - raw reads
     - standard processing
     - MIGEC (UMI consensus)
   * - EHEB / Err1
     - 362
     - 1041–1085
     - 9007–24696
   * - **V1 / Err1**
     - **1.35**
     - 3.1–3.8
     - **26.5–75.9**
   * - **V2 / Err2**
     - **0.28**
     - 1.7–2.0
     - **4.6–6.2**

Read the second row of numbers carefully. **V2 is** *less* **abundant than the worst
2-substitution PCR error** — 0.28× — so no abundance threshold can keep V2 and drop that error,
because the error is bigger. V1 is only 1.35× its worst competitor. This is the regime where
abundance-based error correction is provably unsolvable, and it is precisely why molecular
barcodes exist: consensus assembly moves V1/Err1 from ~1.4 to 26–76, which is a change in the
*evidence*, not in a threshold.

.. note::

   Measured independently on ``SRR1200517`` by the arda benchmark project, which found the same
   ordering on a 1.95 M-read subsample and concluded that the missing capability was UMI consensus
   itself. Those numbers are a subsample and are not like-for-like against the paper's; the
   ordering and the conclusion do not depend on that.

Running it
----------

.. code-block:: bash

   python scripts/spikein_ratio.py reads.fq.gz                     # baseline, raw reads
   python scripts/spikein_ratio.py consensus.fq.gz --label migec   # after assembly

The junction is located by its conserved 3′ anchor, so this runs on raw reads with no reference
and no V/J calling, and it searches both orientations.

.. warning::

   Anchor on the 3′ end only. Requiring the junction's *first* twelve bases as well looks more
   robust and is catastrophically wrong here: V1 differs from EHEB at position 4 and V2 at
   positions 7–8, i.e. inside a 5′ anchor. Both variants then count as zero and the metric looks
   perfect. This is covered by a test.

Acceptance gate
---------------

``V1/Err1`` in **26.5–75.9** and ``V2/Err2`` in **4.6–6.2** after consensus assembly, on
``SRR1200517``. That is the headline claim of the rewrite, and it is not yet met — consensus
assembly lands in M1. Until then the script reports the raw baseline, which reproduces the
published ~1.4 and ~0.3.

Data
----

``PRJNA239303``, runs ``SRR1200517``–``SRR1200520``. Only Experiment 2 is public; the
12-clonotype TRA/TRB truth of Supplementary Table 1a has no public raw reads and lives on the
cluster. See ``SOURCES.md``.

The shallow regime, on a real repertoire library
------------------------------------------------

Every claim migec makes about **1–3 reads per UMI** — that the count ratio carries nothing there,
that payload agreement is what remains, that singleton filtering is unaffordable — was measured on
simulated libraries and on a deep amplicon. Neither is the shape the tool is actually used on. This
is: one HiSeq lane of 5′-RACE bulk TCR β from an ageing cohort, ten donors multiplexed by a 4 nt
sample tag, 149,588,907 read pairs, a 16 nt UMI, published in Britanova *et al.*, *J Immunol*
2014;192(6):2689–98 (`10.4049/jimmunol.1302064
<https://doi.org/10.4049/jimmunol.1302064>`_) — retrieved from PubMed.

Demultiplexing recovered **36.0%** of the lane into the ten declared samples at **2.35–2.62 reads
per UMI**, which is the regime. Four donors then went through correction and consensus:

.. list-table::
   :header-rows: 1
   :widths: 12 13 13 11 13 13 11

   * - sample
     - reads
     - barcodes
     - merged
     - by payload
     - molecules
     - ε at depth
   * - A2-i129
     - 4,832,212
     - 1,940,846
     - 6.12%
     - 10,819
     - 1,822,070
     - 8.1e-4
   * - A2-i131
     - 4,398,345
     - 1,802,598
     - 6.21%
     - 10,144
     - 1,690,723
     - 1.1e-3
   * - A2-i132
     - 6,190,949
     - 2,411,942
     - 6.23%
     - 13,716
     - 2,261,640
     - 1.6e-3
   * - A2-i133
     - 4,921,328
     - 1,997,079
     - 7.18%
     - 13,083
     - 1,853,650
     - 1.4e-3

``assets/shallow_repertoire.tsv`` is the full table.

What it says, in order of how much it matters:

**The payload term is not decoration.** Between 10,144 and 13,716 merges per donor — about 9% of
all merges — were ones *the count ratio alone would have refused*. At 2.5 reads per UMI a parent
with two reads and a child with one is not an asymmetry, and without the reads agreeing on the
molecule those barcodes stay separate and inflate the count. See :doc:`refine`.

**The distance-1 estimator correctly refused to answer.** It read 1.3e-2, sixteen times the
children-based estimate, and flagged itself: 10% of every barcode's 1-substitution neighbourhood is
occupied at these densities, so the excess it measures is a small difference of two large numbers.
The estimate from the children of well-covered parents — 8.1e-4 to 1.6e-3, Q28–Q31 — is the one to
read, and `err_unreliable` is what says so. See :doc:`umi_errors`.

**A 16 nt UMI was worth 12.1–13.0.** The effective length is three to four bases below the nominal
one, which is a 100–1000× smaller barcode space than the sheet implies, and it is why the collision
correction declined: at 8.7–11.9% occupancy the space is saturated and the inversion would collapse
onto the observed count. See :doc:`barcode_space`.

**The residual FDR is 5.0–5.7% of one-read molecules**, and at ≥2 reads it is inside the 5% target.
Reported, never applied.

Note: the numbers above are from the published wheel of the previous release, run on the cluster
that holds the reads — so they are also a check that the artefact on PyPI works on a machine nobody
developed on. Peak RSS was 932 MB for the demultiplex, of which 568 MB was the UMI counters, on a
lane that is 21 GB compressed: the allocation that :doc:`performance` bounds, measured on the
library that produces it.

MIGEC 1.2.9, the implementation this replaced
---------------------------------------------

The rewrite's claim is that the *algorithms* are the specification and the *code* is not, so the
comparison worth running is against the Groovy MIGEC this repo replaced — the same barcode
dialect, the same sheet, the same simulated library, both pipelines end to end.

.. code-block:: bash

   gh release download 1.2.9 --repo antigenomics/migec -p 'migec-1.2.9.zip'
   python scripts/compare_migec_v1.py --out /tmp/v1 --jar migec-1.2.9.jar \
       --molecules 20000 --clones 200 --coverage 8 --min-count 1

20,000 molecules over 200 clones, 149,103 reads, a 12 nt barcode at a 2·10⁻³ per-base error rate.
v1 is ``Checkout`` → ``Assemble --filter-collisions``; migec 2 is ``checkout`` → ``refine`` →
``assemble``. ``assets/migec_v1.tsv``:

.. list-table::
   :header-rows: 1
   :widths: 12 18 16 14 14 12 14

   * - min reads
     - tool
     - consensuses
     - exactly a template
     - precision
     - seconds
     - peak RSS
   * - 1
     - MIGEC 1.2.9
     - 22,717
     - 21,338
     - 0.9393
     - 2.95
     - 974 MB
   * - 1
     - **migec 2**
     - **20,017**
     - **19,977**
     - **0.9980**
     - **0.31**
     - **267 MB**
   * - 5
     - MIGEC 1.2.9
     - 16,292
     - 15,490
     - 0.9508
     - 3.54
     - 941 MB
   * - 5
     - **migec 2**
     - **16,636**
     - **16,635**
     - **0.9999**
     - **0.30**
     - **266 MB**

* **9.4–11.9× the wall clock**, against an M5 gate of 3×, and 3.5–3.7× less memory. Single core in
  both cases; migec 2 threads and v1 does not, so this understates it.
* **The molecule count is right.** 20,017 consensuses for 20,000 molecules — 0.09% over. v1 emits
  22,717, 13.6% over, which is the barcode errors ``--filter-collisions`` did not catch: its rule
  is a count ratio, and :doc:`refine` is where the measurement lives that says a count ratio
  carries nothing below ~3 reads per UMI.
* **99.8–99.99% of consensuses are exactly a template**, against 93.9–95.1%.

.. warning::

   ``--min-count`` is applied to **both**. v1 defaults to 5 and migec 2 to 1 — v1 names its output
   ``.t5.`` for exactly this reason — so leaving each at its own default compares defaults rather
   than implementations, and would have credited us with recovering molecules v1 was told to throw
   away.

MAGERI 1.1.1, the other descendant
----------------------------------

MAGERI is the closer comparator of the two: same author, same UMI model, and the assembler this
repo's consensus is a rewrite of. It is also a superset — it assembles, maps and calls variants in
one run — so what is compared is the part they share, and the extra work is named rather than
divided out.

.. code-block:: bash

   gh release download 1.1.1 --repo mikessh/mageri -p mageri.zip
   python scripts/compare_mageri.py --out /tmp/mageri --jar mageri.jar \
       --molecules 20000 --clones 200 --min-count 1

The same library as above: 20,000 molecules over 200 clones, 149,103 reads, a 12 nt barcode at a
2·10⁻³ per-base error rate. ``assets/mageri.tsv``:

.. list-table::
   :header-rows: 1
   :widths: 10 24 15 16 12 11 12

   * - min reads
     - tool
     - consensuses
     - exactly a template
     - precision
     - seconds
     - peak RSS
   * - 1
     - MAGERI 1.1.1
     - 22,707
     - 22,377
     - 0.9855
     - 2.02
     - 2,246 MB
   * - 1
     - **migec 2 + minimap2**
     - **20,017**
     - **19,977**
     - **0.9980**
     - **0.41**
     - **280 MB**
   * - 2
     - MAGERI 1.1.1
     - 19,966
     - 19,930
     - 0.9982
     - 1.98
     - 1,171 MB
   * - 2
     - **migec 2 + minimap2**
     - **19,974**
     - **19,940**
     - **0.9983**
     - **0.39**
     - **279 MB**
   * - 5
     - MAGERI 1.1.1
     - 16,282
     - 16,282
     - **1.0000**
     - 1.98
     - 1,146 MB
   * - 5
     - **migec 2 + minimap2**
     - **16,636**
     - **16,635**
     - 0.9999
     - **0.39**
     - **281 MB**

* **4.9–5.1× the wall clock**, with the alignment MAGERI does folded into migec's row:
  ``minimap2 -ax sr -y`` onto the same reference plus a ``samtools sort`` costs 0.08 s of the 0.39.
  MAGERI's clock still also covers variant calling, which migec does not do at all.
* **4.1–8.0× less memory.** The range is the JVM's, not ours: MAGERI's peak RSS moves between
  1,146 and 2,246 MB on the same input across runs, where migec 2 sits at 279–281 MB across all
  six.
* **At one read per MIG the molecule count separates them.** MAGERI emits 22,707 consensuses for
  20,000 molecules — **13.5% over**, within a tenth of a point of MIGEC 1.2.9's 13.6%, which is
  what a shared count-ratio correction lineage looks like. migec 2 is 0.09% over. At ``--min-count
  2`` the two agree to 0.04%, because a threshold of 2 discards most of what a count ratio cannot
  correct.
* **At five reads per MIG migec keeps 2.2% more molecules at identical recall** (16,636 against
  16,282, both at 0.9873 of the recoverable templates), and one of its 16,636 is not exactly a
  template against none of MAGERI's 16,282. That is the trade the whole of :doc:`refine` is about:
  err on precision when *merging*, and do not throw a molecule away to buy a fourth decimal place.

.. warning::

   The MIG size threshold is applied to **both**, and it is checked afterwards rather than assumed.
   MAGERI's preset carries ``forceOverseq=true`` with ``defaultOverseq=5``, so out of the box it
   assembles only MIGs of five reads or more; ``compare_mageri.py`` rewrites that value from the
   exported preset and then reads back the threshold MAGERI *reports* it used, refusing to score
   the run if the two differ.

.. note::

   MAGERI's sub-clustering test runs at ``pcrMinorTestPValue = 0.01`` — the nominal threshold
   :doc:`nulls` measured as over-calling by 19× once both margins of the reads × positions matrix
   are preserved. It is left at its default here on purpose: a head to head measures what the other
   tool does, not what it would do if it were this one.

The same comparison at the variant level
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MAGERI calls variants; migec stops at the consensus and hands it to a caller. Stopping the
comparison at the consensus therefore stops it one step short of what either tool is *for*. The
simulator's ``variant_af`` turns the clone set into one reference plus five point variants of it at
a known allele fraction — the ctDNA shape, and the only shape in which a call set has an answer —
and the two pipelines are then scored on the calls.

.. code-block:: bash

   python scripts/compare_mageri.py --out /tmp/mageri --jar mageri.jar \
       --molecules 20000 --clones 6 --coverage 1.5 --min-count 1 \
       --variant-af 0.01 --caller lofreq

``assets/mageri_variants.tsv``, 20,000 molecules, 5 variants, one replicate per cell:

.. list-table::
   :header-rows: 1
   :widths: 12 20 12 10 10 10 13 13

   * - reads/UMI
     - pipeline
     - variant AF
     - called
     - true
     - false
     - sensitivity
     - precision
   * - 8
     - MAGERI
     - 1%
     - 5
     - 5
     - 0
     - 1.000
     - 1.000
   * - 8
     - migec 2 + LoFreq
     - 1%
     - 5
     - 5
     - 0
     - 1.000
     - 1.000
   * - 8
     - MAGERI
     - 0.1%
     - 5
     - 5
     - 0
     - 1.000
     - 1.000
   * - 8
     - migec 2 + LoFreq
     - 0.1%
     - 5
     - 5
     - 0
     - 1.000
     - 1.000
   * - 1.5
     - MAGERI
     - 1%
     - 142
     - 5
     - **137**
     - 1.000
     - 0.035
   * - 1.5
     - **migec 2 + LoFreq**
     - 1%
     - **5**
     - **5**
     - **0**
     - **1.000**
     - **1.000**
   * - 1.5
     - MAGERI
     - 0.1%
     - 142
     - 5
     - **137**
     - 1.000
     - 0.035
   * - 1.5
     - **migec 2 + LoFreq**
     - 0.1%
     - 1
     - 1
     - **0**
     - 0.200
     - **1.000**

* **At 8 reads per UMI the two are indistinguishable**: 5 of 5 with no false positive, at 1% and at
  0.1%, and the called allele fraction is within 1–4·10⁻⁴ of the injected one for both. A consensus
  over 8 reads removes essentially all sequencing error, and what is left is set by the molecule
  count, which is the same number for both tools.
* **At 1.5 reads per UMI they separate by 28×, and it is the calling that separates them, not the
  collapsing.** Both tools emit the *same* 20,170 consensuses at the same accuracy — 0.8990 of
  MAGERI's are exactly a template against 0.8985 of migec's — so the input to the two callers is
  matched to within 5·10⁻⁴. On that input MAGERI reports 142 variants of which 137 are at positions
  nothing was injected at; LoFreq on migec's consensus reports 5 and is right about all of them.
* **At 0.1% and 1.5 reads per UMI the trade is explicit**: MAGERI buys 4 extra true positives with
  137 false ones, a positive predictive value of 3.5%. migec 2 finds 1 of 5 and is right about it.
  Which of those a study wants is the study's choice, and :doc:`detection` is where the arithmetic
  for making it lives — but a caller that returns 137 false positives per sample cannot be run
  without a per-position background model, which is exactly what that page concludes.

.. note::

   The variant arm ran on the cluster rather than on this laptop, because LoFreq is what migec
   would hand its consensus to and that is where it is installed (``scripts/compare_mageri.sbatch``,
   SLURM, 8 cores). Its wall clocks are therefore not comparable with the consensus table above,
   which was measured on one machine end to end; on the cluster the ratio is 2.3–3.7× including
   LoFreq, which migec's row pays and MAGERI's does not.
