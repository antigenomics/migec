Roadmap
=======

Milestones are ordered by risk, not by pipeline order: the consensus quality model is the
scientific claim, so it is validated before any of the throughput work is done.

This page is the summary. ``ROADMAP.md`` in the repository root carries the per-item state and is
the one kept current as work lands.

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
     - ``checkout``: patterns, trimming, header transfer, UMI statistics, paired input,
       strand normalisation, multi-core, whitelists, dual-end barcodes
     - done bar ``.mig`` bucket output, i7xi5 and the bit-parallel matcher
   * - M1
     - ``assemble``: consensus, sub-clustering, quality model, ``--contig``, ``--fast``
     - done bar ``--rt-error auto`` and R1/R2 overlap merge
   * - M3
     - ``refine``: error model, barcode correction, cell calling, QC
     - done bar the template's own error split
   * - M4
     - end-to-end, ``suggest``, ``subsample``, ``plot``, notebooks
     - done, downstream contract measured
   * - M5
     - benchmarks, ``isalgo/umi_data``, release
     - in progress

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
