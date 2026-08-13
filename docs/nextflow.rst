Pipelines: Nextflow and SLURM
=============================

Two ways to run migec over more than one sample. They answer different questions: Nextflow when
the pipeline continues past the consensus into alignment, calling or AIRR; SLURM when you want
three commands over a cohort and nothing else.

Nextflow
--------

``integrations/nextflow/`` is an nf-core-style local module set with a runnable entry point.

.. code-block:: text

   main.nf                                  --mode consensus | ctdna | airr

   modules/migec/checkout/main.nf           reads     -> tagged FASTQ + QC tables
   modules/migec/refine/main.nf             tagged    -> corrected FASTQ + barcode table + cells
   modules/migec/assemble/main.nf           corrected -> one consensus per molecule
   modules/migec/plot/main.nf               the tables -> SVG figures (needs gnuplot)

   modules/downstream/align/main.nf         consensus -> BAM, tags carried and checked
   modules/downstream/callvariants/main.nf  BAM       -> VCF (LoFreq or Mutect2)
   modules/downstream/arda/main.nf          consensus -> AIRR clonotypes

   subworkflows/migec/main.nf               the three migec stages chained
   subworkflows/migec_ctdna/main.nf         + align + call: rare somatic variants
   subworkflows/migec_airr/main.nf          + arda: immune repertoires
   nextflow.config                          defaults, every one read with getOrDefault

.. code-block:: bash

   nextflow run integrations/nextflow -profile docker \
       --mode ctdna --input 'data/*_R{1,2}.fq.gz' --preset tso500 \
       --fasta ref.fa --outdir results/

   nextflow run integrations/nextflow --mode airr --input 'data/*_R{1,2}.fq.gz' --preset migec

Or include a subworkflow in a pipeline of your own:

.. code-block:: groovy

   include { MIGEC } from './integrations/nextflow/subworkflows/migec/main'

   workflow {
       ch_reads = Channel.fromFilePairs(params.input)
           .map { id, files -> [ [ id: id, preset: '10x-v2', payload_mate: 2 ], files ] }

       MIGEC(ch_reads, file(params.cell_whitelist ?: 'NO_FILE'))

       MIGEC.out.consensus.view()
   }

Per-sample keys in ``meta`` win over the ``params.*`` defaults, so one run can mix chemistries:
``bc_pattern``, ``preset``, ``read_structure``, ``read_structure2``, ``barcodes``, ``max_offset``,
``payload_mate``, ``expect_cells``, ``rt_error``, ``contig``, ``fast``, ``aligner``, ``caller``,
``species``.

What the downstream modules encode
----------------------------------

Two rules, both from :doc:`variants`:

**Collapse first, then align once.** Aligning raw reads and grouping on *(position, UMI)* is the
other order in use, and it costs one alignment per *read* rather than per *molecule*, with the
aligner seeing uncorrected sequence.

**A standard variant caller, never a UMI-aware one.** After ``assemble`` a caller's depth already
is a molecule count. ``UMI-VarCal`` and ``UMIErrorCorrect`` group and consense themselves, so they
replace ``assemble`` rather than following it.

Never: **do not set a family-size filter downstream of** ``assemble``. Every family has size 1 by
construction, so ``--min-family-size 3`` discards the entire library and reports zero variants
without an error.

The align module checks its own output: if no ``MI:Z:`` tag reaches the BAM it exits non-zero and
names the flag each aligner needs (``-y`` for minimap2 and ``minibwa map``, ``-C`` for
``bwa mem``). Without that check the failure is silent and surfaces much later as an untagged BAM.

SLURM
-----

``integrations/slurm/`` is two sbatch templates and a sample sheet, for the case where the
deliverable is the consensus and a workflow engine would be ceremony.

.. code-block:: bash

   # one sample
   sbatch --export=ALL,R1=s1_R1.fq.gz,R2=s1_R2.fq.gz,SAMPLE=s1,PRESET=10x-v2,PAYLOAD_MATE=2 \
          integrations/slurm/migec_sample.sbatch

   # a cohort, one array task per row
   sbatch --array=1-$(($(wc -l < samples.tsv) - 1)) \
          integrations/slurm/migec_array.sbatch samples.tsv

Both run **without** SLURM as ordinary bash scripts -- every SLURM variable has a fallback -- which
is how they are tested and how a layout should be checked before a cohort is queued:

.. code-block:: bash

   R1=s1.fq.gz SAMPLE=s1 BC_PATTERN='0:12' bash integrations/slurm/migec_sample.sbatch

Note: **array task 1 is the first data row.** The header is skipped rather than counted, so the
range is ``1-(rows - 1)``. An ``--array=0-N`` runs a task that reads the header as a sample and
fails somewhere confusing instead of at the sheet.

Sizing the request: ``checkout`` scales with reads (give it cores), ``refine`` with *distinct
barcodes* (give it memory; ``table_bytes`` in its JSON sizes the next run), ``assemble`` with one
bucket at a time (cores; peak memory is set by the bucket count, not the library). 16 cores and
32 GB covers a typical targeted or single-cell library.

Never: **-t changes the wall clock and nothing else.** Every stage is byte-identical at any thread
count, so a retry on a different node or an escalating ``--requeue`` cannot produce a result that
disagrees with the first attempt. That is what makes automatic retries safe here, and
``tests/benchmark/`` asserts it.

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
     - 1,548,835 reads/s at 16
     - chunk-bounded, plus the UMI counters
   * - ``refine``
     - **distinct barcodes**
     - 1,554,156 reads/s at 16
     - the barcode table, ~96 B each
   * - ``assemble``
     - reads, then buckets
     - 2,470,928 reads/s at 16
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

   The SLURM templates are the opposite case: they *have* been run end to end here, as ordinary
   bash, producing a consensus and its figures, and their error paths were exercised too. Between
   the two, the migec commands themselves are covered either way -- both call the same three.
