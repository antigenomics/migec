Fragmented libraries — why 10x needs a different consensus
===========================================================

MIGEC assumes a MIG is one molecule amplified from fixed primers, so its reads start at the same
base and one ungapped consensus is the right model. A 3' GEX library is not that: the molecule is
captured, then fragmented, so reads sharing a ``(CB, UMI)`` tile it at different offsets and may
share no sequence at all.

That assumption is cheap to falsify, and it had to be falsified before ``assemble`` was designed
rather than after — hence experiment **X1**.

What was measured
-----------------

``scripts/read_start_dispersion.py`` groups aligned reads by ``(CB, UB, contig, strand)`` and asks
two questions of every group holding more than one read: do the reads start at the same base, and
do they form one connected component under "shares at least one aligned base"?

.. code-block:: bash

   uv pip install pysam    # not a migec dependency; X1 is a one-off
   python scripts/read_start_dispersion.py \
       --bam https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam \
       --region 11:65497688-65508073 --region 7:5527151-5530601 --region 15:44711477-44718877

The BAM is read over HTTPS by byte range, so this costs a few hundred megabytes rather than the
4.8 GB the file weighs.

.. important::

   Overlap is computed on **aligned blocks**, not on ``reference_start``–``reference_end``. These
   are spliced alignments: a read's genomic span includes its introns, so two reads either side of
   a junction can span the same megabase while sharing no aligned base. The first pass of this
   experiment used genomic spans and reported reads starting 500 kb outside a 10 kb window, which
   is what exposed the mistake.

The answer
----------

MALAT1, ACTB and B2M in ``pbmc_1k_v3`` — 1.94 M reads, 560 k carrying both ``CB`` and ``UB``,
540,619 distinct ``(CB, UMI)`` groups of which 8,326 hold more than one read:

.. list-table::
   :header-rows: 1
   :widths: 16 14 20 22 22

   * - reads/UMI
     - groups
     - co-terminal
     - one component
     - more than one
   * - 2
     - 3,582
     - 14.5%
     - 71.9%
     - 28.1%
   * - 3
     - 1,778
     - 5.5%
     - 69.3%
     - 30.7%
   * - 4
     - 1,272
     - 2.0%
     - 71.9%
     - 28.1%
   * - 5
     - 810
     - 0.6%
     - 75.8%
     - 24.2%
   * - 6+
     - 884
     - 0.3%
     - 81.3%
     - 18.7%
   * - **all**
     - **8,326**
     - **7.8%**
     - **72.7%**
     - **27.3%**

Three things follow, and they are the design of ``assemble`` for 10x:

**Co-terminality is not merely rare, it is an artefact of small groups.** 14.5% at two reads,
0.3% at six or more. That is the signature of independently placed reads: with two reads a
coincidence is possible, with six it is not. Any offset-scan consensus keyed on a shared start —
MIGEC's ``--max-offset`` — has essentially nothing to work with here. 92% of groups have a
footprint wider than a single 91 nt read.

**But most groups are still one overlap component**, and the fraction *rises* with depth (72% at
two reads, 81% at six or more) because more reads mean more chances to bridge a gap. So the
fragmented case is not hopeless: partition by overlap first and roughly three quarters of groups
reduce to exactly the ungapped consensus problem MIGEC already solves.

**The remaining quarter is why the partition is mandatory rather than optional.** 27.3% of
multi-read groups split into two or three components. Handing those to a single ungapped consensus
does not produce a slightly worse answer — it produces a sequence asserted across a gap that no
read covers. That is fabricated sequence with a consensus quality attached to it, which is worse
than no output.

.. note::

   Only **1.5%** of ``(CB, UMI)`` groups hold more than one read at this depth. Consensus barely
   applies to shallow 3' GEX at all: the overwhelming majority of molecules are seen once, and for
   those the "consensus" is the read. The value of UMIs here is counting, not error correction —
   which is a different claim from the one the amplicon literature makes, and it should not be
   blurred in anything migec publishes.

Consequences for M1
-------------------

* ``--mode {amplicon,fragmented}`` is explicit, not inferred. The two make opposite assumptions and
  guessing wrong is silent in both directions.
* ``fragmented`` partitions each group into overlap components by union-find and emits one
  consensus per component. The contig path — "assemble overlapping read groups sharing a UMI" —
  is that same code, not a separate feature.
* A component must never be extended across a gap. Emit the components; let the caller decide
  whether to scaffold them.
