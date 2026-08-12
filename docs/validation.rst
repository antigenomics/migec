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
