Examples
========

Three commands, whatever the platform. Only ``checkout`` needs to know where the barcode is.

.. code-block:: bash

   migec checkout  reads.fq.gz -b barcodes.txt -o co/   # find and cut out the barcode
   migec refine    co/S1.fq.gz -o ref/                  # fix errors IN the barcode
   migec assemble  ref/S1.fq.gz -o asm/                 # collapse each molecule

Declaring the layout
--------------------

Four equivalent ways to say where the barcode is, in the order to reach for them:

.. list-table::
   :header-rows: 1

   * - form
     - looks like
     - when
   * - a position
     - ``^NNNNNNNN`` or ``0:8``
     - the primary mode -- the chemistry fixes the barcode at an offset
   * - ``--preset``
     - ``10x-v2``, ``tso500``, ``duplex``, ...
     - a chemistry with a name; ``migec sheet --presets`` lists them
   * - ``--read-structure``
     - ``5M5S+T``
     - fgbio, Picard, samtools and TSO500 all speak this
   * - barcode table
     - ``S1<TAB>aaACTcagtgg...NNNNtNNNNtNNNN``
     - many samples in one file; MIGEC's own format, read verbatim

In a pattern, ``N`` is a UMI base, ``X`` a cell-barcode base, uppercase is matched exactly,
lowercase is the fuzzy adapter region, and ``.`` is skipped. Slices are half-open and 0-based like
Python's, each a UMI slice unless prefixed ``cell:``. A leading ``^``, a slice list and a read
structure all anchor the barcode at the first base, so ``--max-offset`` never has to be passed.

:doc:`layouts` has all of it in one place. Do not guess the layout -- :doc:`suggest` reads it off
the data, and :doc:`downstream` is what happens after.

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
     - ``--preset 10x-v2``, or ``--bc-pattern 'cell:0:16,16:26'``
     - positional; ``refine``/``assemble`` take **R2**
   * - TSO500
     - ``--preset tso500``, or ``--read-structure 5M5S+T``
     - ``5M5S+T +T``: a 5 nt UMI on R1 only. Warning: 1,024 barcodes does not identify
       a molecule -- TSO500 groups position-aware, downstream
   * - UMI RNA-seq (SMARTer)
     - ``--preset smarter-umi``
     - 10 nt inline UMI, then the ``GGG`` the template switch leaves
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
