migec
=====

UMI barcode extraction, correction and consensus assembly for barcoded sequencing data — a
complete C++20 rewrite of `MIGEC <https://doi.org/10.1038/nmeth.2960>`_ and
`MAGERI <https://doi.org/10.1371/journal.pcbi.1005480>`_.

.. warning::

   **Version 2 is under construction.** ``checkout`` works today — barcode extraction, trimming,
   header transfer and the UMI statistics. ``refine`` and ``assemble`` land over the following
   milestones; see :doc:`roadmap`. The Groovy MIGEC 1.2.9 that this replaces is archived on the
   ``legacy-v1`` branch and at tag ``v1-final``.

What it does
------------

Molecules are tagged with a random barcode (a UMI) before amplification, so every read carrying
the same UMI descends from one original molecule. Collapsing those reads into a consensus removes
essentially all sequencing error, which is what makes rare-variant and repertoire work possible.
Doing it correctly is harder than it looks:

* barcodes themselves acquire errors, and telling an error-child barcode from a genuine collision
  needs the birthday bound, the base qualities, and the fact that a polymerase error in an early
  PCR cycle carries *high* quality in every read that inherits it;
* a molecule seen 3–5 times is still information. Throwing it away because it is under a coverage
  threshold discards real sequence; migec keeps it and reports the uncertainty instead;
* consensus cannot fix an error made during reverse transcription or the first PCR cycle, because
  it is present in every read. Any quality score above that floor is a fiction, so migec measures
  the floor from the data and refuses to emit a quality above it.

Pipeline
--------

.. code-block:: text

    FASTQ ──checkout──▶ .mig ──refine──▶ .mig + .pumi ──assemble──▶ consensus FASTQ
              │                   │                                      │
         suggest              QC tables, plots                    per-molecule tables

The output is ordinary FASTQ with the sample, cell barcode and UMI in the read name and in
SAM-style tags, so ``minimap2``, ``bwa``, ``arda``, ``salmon`` and ``kallisto`` consume it
directly -- see :doc:`downstream`, where each of those was run against it.

The only thing that changes between platforms is **where the barcode is**, and the primary way to
say that is a position:

.. code-block:: bash

    migec checkout reads.fq.gz --bc-pattern '^NNNNNNNN' -o out/
    migec checkout reads.fq.gz --bc-pattern '0:8'       -o out/
    migec checkout R1.fq.gz R2.fq.gz --preset 10x-v2    -o out/

See :doc:`layouts` for the presets, the slice grammar, fgbio read structures and barcode tables.

.. toctree::
   :maxdepth: 2

   installation
   examples
   layouts
   suggest
   checkout
   downstream
   refine
   assemble
   subsample
   umi_statistics
   barcode_space
   performance
   grouping
   fragmented
   quality_floor
   nulls
   validation
   formats
   roadmap
