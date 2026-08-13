Commands
========

Five commands move reads, and three read no reads at all. That count is fixed: a sixth pipeline
command needs a failing benchmark that the existing five cannot pass, because every flag and every
stage is a thing that can disagree with another one.

The pipeline
------------

.. code-block:: bash

   migec checkout  reads.fq.gz --bc-pattern '^NNNNNNNN' -o co/   # find and cut out the barcode
   migec refine    co/S1.fq.gz                          -o ref/  # fix the errors IN the barcode
   migec assemble  ref/S1.fq.gz                         -o asm/  # collapse each molecule

.. list-table::
   :header-rows: 1
   :widths: 18 34 48

   * - command
     - what it does
     - the number it decides
   * - :doc:`suggest <suggest>`
     - reads the barcode layout off the data
     - the pattern you paste into ``checkout``
   * - :doc:`checkout <checkout>`
     - finds the pattern, cuts the barcode out, demultiplexes
     - how many reads carry a usable barcode
   * - :doc:`refine <refine>`
     - corrects errors in the barcode itself, calls cells
     - **how many molecules there were**
   * - :doc:`assemble <assemble>`
     - collapses each molecule's reads into one consensus
     - the sequence, and the quality it is allowed to claim
   * - :doc:`subsample <subsample>`
     - keeps every read of a fraction of the barcodes
     - the size of a fixture that still behaves like the library

``suggest`` and ``subsample`` are outside the pipeline in different directions: one runs before it
to tell you what to type, the other cuts a library down to something a laptop can iterate on.

Reading the output
------------------

.. list-table::
   :header-rows: 1
   :widths: 18 34 48

   * - command
     - what it does
     - reads no reads because
   * - :doc:`plot <plots>`
     - draws sixteen QC panels from the TSVs the stages wrote
     - it computes nothing; a figure that cannot be redrawn from a committed table will
       eventually disagree with the report
   * - ``sheet``
     - writes or explains a barcode table, lists the presets
     - it only ever touches the layout
   * - ``info``
     - prints what a ``.mig`` bucket or a tagged FASTQ contains
     - it reads headers

Every stage writes plain TSV next to its output and a ``.json`` summary beside that, so the
:doc:`figures <plots>` and any downstream script see the same numbers the report printed.

.. toctree::
   :maxdepth: 1

   suggest
   checkout
   refine
   assemble
   subsample
   plots
