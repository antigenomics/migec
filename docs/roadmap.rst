Roadmap
=======

Milestones are ordered by risk, not by pipeline order: the consensus quality model is the
scientific claim, so it is validated before any of the throughput work is done.

.. list-table::
   :header-rows: 1
   :widths: 12 60 28

   * - milestone
     - content
     - state
   * - M0
     - repo, ``.mig`` format, FASTQ IO, read simulator
     - done
   * - M2
     - ``checkout``: patterns, trimming, header transfer, UMI statistics
     - single-end done
   * - M1
     - ``assemble``: consensus, sub-clustering, quality model
     - planned
   * - M3
     - ``refine``: error model, barcode correction, QC
     - planned
   * - M4
     - end-to-end, ``suggest``, ``sort``, ``subsample``, notebooks
     - planned
   * - M5
     - benchmarks, ``isalgo/umi_data``, release
     - planned

Explicitly out of scope for v2.0
--------------------------------

* **Alignment and variant calling.** MAGERI's job. The pipeline ends at consensus FASTQ and hands
  off to arda, minimap2 or bwa-meme.
* **Indels.** Illumina indel rates are around :math:`10^{-6}` per base and there is no dataset in
  the benchmark set that would let us verify indel handling. Substitutions only, everywhere.
* **Duplex consensus (DCS).** v2.0 extracts duplex *tags* and emits single-strand consensuses.
  Pairing the two strand families is a later addition, and until it exists no error-suppression
  claim here is based on duplex data.
* **EmptyDrops-style cell rescue.** Cell calling is OrdMag plus a knee. Reproducing Cell Ranger's
  second pass is its job, not ours.
