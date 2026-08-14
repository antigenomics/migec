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
