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
   * - Capture, exome, ctDNA, MRD
     - nothing -- the UMI is already in ``RX``
     - the kit puts it in the index read, so there is no layout and no ``checkout``. See
       :doc:`Bring your own UMI <byo_umi>`

Starting from a BAM
-------------------

A hybrid-capture kit reads the UMI on the index, so what arrives is an fgbio, Picard or vendor BAM
with the UMI in ``RX`` — never a FASTQ with a barcode inside the read. Two commands, no layout,
no ``checkout``:

.. code-block:: bash

   migec refine   tagged.bam  -o ref/   # correct errors in the UMI
   migec assemble ref/S1.fq.gz -o asm/  # one consensus per molecule

BAM, SAM and CRAM are recognised from the file, not the name; ``samtools`` does the conversion and
the temporary FASTQ is deleted when the stage returns. To see the round trip on data you already
have:

.. code-block:: bash

   migec checkout reads.fq.gz --bc-pattern '^NNNNNNNNNNNN' -o co/
   samtools import -T '*' -s co/S1.fq.gz -o S1.bam     # the tags become real BAM tags
   migec refine   S1.bam -o from_bam/                  # identical to refining co/S1.fq.gz

Runnable notebooks
------------------

Six `marimo <https://marimo.io>`_ notebooks, each a plain Python file with its own PEP 723
dependency header, so ``uv`` builds the environment and nothing has to be installed first. Two of
them need the network: ``platforms.py`` downloads its fixtures from
`isalgo/umi_data <https://huggingface.co/datasets/isalgo/umi_data>`_ on first run, and
``ctdna_variants.py`` fetches its runs straight from SRA with ``scripts/sra_fetch.py``, because
anything with a public accession is regenerated rather than mirrored. The other four simulate a
library whose true molecule and clone counts are known, so every number they print can be checked
rather than admired.

.. code-block:: bash

   uv run marimo edit notebooks/platforms.py

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
   * - ``notebooks/exome_capture.py``
     - duplicates or real molecules? why coordinate deduplication undercounts a capture panel
   * - ``notebooks/airr_repertoire.py``
     - how much of a repertoire is PCR? clonotype counts from reads against molecules
   * - ``notebooks/ctdna_variants.py``
     - how many molecules a variant caller actually gets, on cell-free DNA reference material at
       known allele frequency (:doc:`variants`)

Looking at the run
------------------

.. code-block:: bash

   migec plot co/ -o figs/          # every panel whose table is in co/
   migec plot ref/ -o figs/
   migec plot asm/ -o figs/ --format pdf

Four of the twenty panels are figures you have already read somewhere else, and they are the ones
to look at first:

.. list-table::
   :header-rows: 1
   :widths: 24 34 42

   * - panel
     - the question
     - the failure it shows
   * - ``cell_rank``
     - is my cell calling right?
     - Cell Ranger's barcode rank plot, on unique UMIs. No knee means no cells.
   * - ``mig_size_spectrum``
     - is the library over-sequenced?
     - most molecules shallow while most reads sit in the deep ones
   * - ``mig_size_zipf``
     - is amplification even?
     - a bent rank curve where Zipf would be straight
   * - ``umi_error_rate``
     - is the barcode error rate believable?
     - the two estimators part company, which means the barcode space filled
   * - ``consensus_quality``
     - what quality am I allowed to claim?
     - the boxes flatten at the RT floor, not at the instrument

Every panel is a gnuplot script over a TSV the stage already wrote, so a figure can be redrawn
without the FASTQ. Without gnuplot installed the ``.gp`` scripts are still written. See
:doc:`plot <plots>`.

Then what
---------

The consensus is ordinary FASTQ, so an aligner or a quantifier takes it directly. All four were run
against real ``assemble`` output; :doc:`downstream` has the table and the record counts.

.. code-block:: bash

   minimap2 -ax sr -y  ref.fa asm/S1.consensus.fq.gz | samtools sort -o S1.bam
   minibwa map -y -t8  ref.fa asm/S1.consensus.fq.gz | samtools sort -o S1.bam
   bwa mem -C          ref.fa asm/S1.consensus.fq.gz | samtools sort -o S1.bam
   salmon quant -i tx.idx -l A -r asm/S1.consensus.fq.gz -o quant/   # NumReads = molecules

Never: not ``alevin``, ``bustools`` or ``STARsolo``. They deduplicate from a raw barcode read that
no longer exists, so running them on a consensus collapses the library twice. :doc:`downstream`
also covers when to align *before* collapsing instead.

Pipelines
---------

Two, for two different situations (:doc:`nextflow` has both in full):

* **Nextflow** -- ``integrations/nextflow/`` runs the three stages and then continues:
  ``--mode ctdna`` aligns and calls variants, ``--mode airr`` calls clonotypes with arda,
  ``--mode consensus`` stops at the consensus. It drops into
  `nf-core/airrflow <https://nf-co.re/airrflow>`_ or any pipeline that hands you FASTQ pairs.
* **SLURM** -- ``integrations/slurm/`` is two sbatch templates and a sample sheet, for a cohort
  where the deliverable is the consensus. Both run as ordinary bash without SLURM, which is how
  to check a layout before queueing anything.

All three stages are byte-identical at any thread count, so a retry with more cores cannot change
a result -- which is what makes an escalating retry safe in either.

.. code-block:: bash

   nextflow run integrations/nextflow --mode ctdna --input 'd/*_R{1,2}.fq.gz' --fasta ref.fa
   sbatch --array=1-11 integrations/slurm/migec_array.sbatch samples.tsv
