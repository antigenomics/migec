Nextflow
========

``integrations/nextflow/`` is an nf-core-style local module set: three processes, a plotting
process, a subworkflow that chains them, and a config of defaults.

.. code-block:: text

   modules/migec/checkout/main.nf     reads     -> tagged FASTQ + QC tables
   modules/migec/refine/main.nf       tagged    -> corrected FASTQ + barcode table + cell calls
   modules/migec/assemble/main.nf     corrected -> one consensus per molecule
   modules/migec/plot/main.nf         the tables -> SVG figures (needs gnuplot)
   subworkflows/migec/main.nf         the three chained, with the payload mate handled
   nextflow.config                    defaults, every one read with getOrDefault

.. code-block:: groovy

   include { MIGEC } from './integrations/nextflow/subworkflows/migec/main'

   workflow {
       ch_reads = Channel.fromFilePairs(params.input)
           .map { id, files -> [ [ id: id, preset: '10x-v2', payload_mate: 2 ], files ] }

       MIGEC(ch_reads, file(params.cell_whitelist ?: 'NO_FILE'))

       MIGEC.out.consensus.view()
   }

Per-sample keys in ``meta`` win over the ``params.migec_*`` defaults, so one run can mix
chemistries: ``bc_pattern``, ``preset``, ``read_structure``, ``read_structure2``, ``barcodes``,
``max_offset``, ``payload_mate``, ``expect_cells``, ``rt_error``, ``contig``, ``fast``.

Why three processes
-------------------

The stages have different shapes, and one process forces one answer for all three:

.. list-table::
   :header-rows: 1
   :widths: 16 26 28 30

   * - stage
     - scales with
     - threads
     - memory
   * - ``checkout``
     - reads
     - 1,056,472 reads/s at 16
     - chunk-bounded, plus the UMI counters
   * - ``refine``
     - **distinct barcodes**
     - 1,012,368 reads/s at 16
     - the barcode table, ~96 B each
   * - ``assemble``
     - reads, then buckets
     - 2,051,937 reads/s at 16
     - one bucket per worker

Splitting them means a failed ``assemble`` resumes without re-running the demultiplex, and each
stage gets the label and the retry that fits it. ``refine`` carries ``process_high_memory`` because
its memory is set by the number of distinct barcodes and by nothing else: a 200 GB shallow run and
a 200 GB deep one need wildly different amounts of it, and FASTQ size predicts neither.

Retries are safe. Every stage's output is byte-identical at any thread count
(:doc:`performance`), so an attempt that gets 4 cpus and an attempt that gets 32 produce the same
bytes -- which is what makes ``errorStrategy 'retry'`` with escalating resources sound here rather
than merely convenient.

.. code-block:: groovy

   process {
       withLabel: process_high        { cpus = 16; memory = 32.GB }
       withLabel: process_medium      { cpus = 8;  memory = 16.GB }
       withLabel: process_high_memory { memory = { 64.GB * task.attempt } }
   }

The three things that go wrong
------------------------------

**The barcode read is not always the payload read.** On 10x, R1 is 26 nt of cell barcode and UMI
and nothing else. ``checkout`` is given both mates; the later stages must then run on the mate that
carries cDNA. ``payload_mate: 2`` says so. Assuming R1 produces empty consensuses on every droplet
chemistry, and nothing in the run reports it as an error.

**A per-sample `false` is not an absent value.** Groovy's ``?:`` treats ``false`` as absent, so
``[ id: 's1', contig: false ]`` against ``params.migec_contig = true`` would silently mean its
opposite. The modules use ``meta.containsKey('contig') ? meta.contig : params...`` for every
boolean, which is the one direction a per-sample override exists to make possible.

**``--rt-error`` names a chemistry, not a number.** ``rt`` (1e-4, caps at Q40) for anything with a
reverse transcription step, ``medium`` (1e-5) for an ordinary polymerase and no RT, ``high``
(1e-6) for a proofreading one. It is the *one-molecule* floor -- 10x's Q60 requires two UMIs to
agree, and combining molecules is `arda <https://github.com/antigenomics/arda>`_'s job, downstream
of this. See :doc:`quality_floor`.

Smoke tests
-----------

Every process has a ``stub:`` block, so ``nextflow run -stub-run`` walks the whole graph without
migec installed -- the test to run before the real one.

``params.migec_limit_read`` and ``params.migec_limit_umi`` pass through to the stages that take
them (``--limit-read`` on all three, ``--limit-umi`` on ``refine`` and ``assemble``), which gets an
answer out of a 400 GB run in a minute.

.. warning::

   A limit is not a sample, and every limited run says so in its own report. Use
   :doc:`subsample` when the output has to behave like the library.

.. note::

   The ``-stub-run`` and the real run have **not** been executed on the machine these docs were
   written on -- nextflow is not installed there. The modules are reviewed against the nf-core
   module spec, not verified by a pipeline run; treat the first run in your own pipeline as the
   verification.
