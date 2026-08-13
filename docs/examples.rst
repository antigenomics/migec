Examples
========

Three commands, whatever the platform. Only ``checkout`` needs to know where the barcode is.

.. code-block:: bash

   migec checkout  reads.fq.gz -b barcodes.txt -o co/   # find and cut out the barcode
   migec refine    co/S1.fq.gz -o ref/                  # fix errors IN the barcode
   migec assemble  ref/S1.fq.gz -o asm/                 # collapse each molecule

Declaring the layout
--------------------

Three equivalent ways to say where the barcode is:

.. list-table::
   :header-rows: 1

   * - form
     - looks like
     - when
   * - barcode table
     - ``S1<TAB>aaACTcagtgg...NNNNtNNNNtNNNN``
     - many samples in one file; MIGEC's own format, read verbatim
   * - ``--bc-pattern``
     - ``XXXXXXXXXXXXXXXXNNNNNNNNNN``
     - one sample, inline; what ``umi_tools``, ``umitools`` and ``mgatk`` take
   * - ``--read-structure``
     - ``5M5S+T``
     - fgbio, Picard, samtools and TSO500 all speak this

In a pattern, ``N`` is a UMI base, ``X`` a cell-barcode base, uppercase is matched exactly,
lowercase is the fuzzy adapter region, and ``.`` is skipped.

Do not guess the layout. :doc:`suggest` reads it off the data.

By platform
-----------

.. list-table::
   :header-rows: 1

   * - platform
     - how to declare it
     - notes
   * - Bulk amplicon (MIGEC, RepSeq)
     - ``aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN``
     - the primer anchors it, so a free scan places it
   * - HIV-1 Primer ID
     - ``NNNNNNNNNcagtttaacttttgggccatcca``
     - recovered from the data by ``migec suggest``
   * - 10x droplet
     - ``--bc-pattern XXXXXXXXXXXXXXXXNNNNNNNNNN --max-offset 0``
     - positional; ``refine``/``assemble`` take **R2**
   * - TSO500
     - ``--read-structure 5M5S+T --read-structure2 5M5S+T``
     - a UMI on both mates, concatenated into one identifier
   * - Dual-end (MAGERI)
     - column 3 of the sheet, ``NNNNNNNNNNNNtgact`` / ``agtcaNNNNNNNNNNNN``
     - both halves must match
   * - Shallow bulk (1-3 reads/UMI)
     - any of the above
     - nothing changes; what the numbers can mean does

Runnable notebooks
------------------

`marimo <https://marimo.io>`_ notebooks, each self-contained. The fixtures they use come from
`isalgo/umi_data <https://huggingface.co/datasets/isalgo/umi_data>`_ and download on first run.

.. code-block:: bash

   marimo edit notebooks/platforms.py

.. list-table::
   :header-rows: 1

   * - notebook
     - what it answers
   * - ``notebooks/platforms.py``
     - every layout above, with a real end-to-end run of two of them
   * - ``notebooks/barcode_space.py``
     - is my barcode long enough? collisions, occupancy, the error budget
   * - ``notebooks/refine_diagnostics.py``
     - the coverage curve, the barcode-rank plot, and where the errors are

Pipelines
---------

``integrations/nextflow/migec/`` is an nf-core-style local module for
`nf-core/airrflow <https://nf-co.re/airrflow>`_ or any pipeline handing you FASTQ pairs. SLURM is
the pipeline's business; the module declares ``label`` and ``task.cpus`` and nothing more.

Note: only ``checkout`` threads. ``refine`` and ``assemble`` are single-threaded by construction,
so ask for the cores ``checkout`` can use and no more.
