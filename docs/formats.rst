File formats
============

This page is the contract between stages. It is frozen before the stages that use it are written,
and ``tests/cpp/test_mig_record.cpp`` fails if any of it changes by accident.

What a stage will read
----------------------

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - format
     - recognised by
   * - FASTQ
     - anything else. Plain or gzipped, decided by the first two bytes rather than the suffix,
       because half the world writes gzipped data to a name ending ``.fastq``
   * - ``.mig``
     - the ``MIGB`` magic; ``refine`` and ``assemble`` only, as a file, a directory or a glob
   * - BAM, SAM, CRAM
     - ``BAM\1`` after inflating the first block, ``CRAM`` at byte 0, or a leading ``@HD``/``@SQ``.
       Converted once with ``samtools`` into a temporary FASTQ inside the output directory, then
       read as one — see :doc:`Bring your own UMI <byo_umi>`

There is no ``--format`` flag and there will not be one: the file says what it is, and a flag that
can disagree with the file is a flag that will.

.. _mig-format:

The ``.mig`` intermediate
-------------------------

One format between all stages. ``checkout --mig`` writes it, ``assemble`` reads it, and
``assemble`` also writes it as the temporaries of its own partition pass when the input was FASTQ.
The file name is ``<sample>.<bbb>.mig``: one file per sample per range-partition bucket, with the
bucket index zero-padded so that a directory listing is in key order.

Layout::

    [FileHeader] [Block]* [Terminator][u64 n_records]["MIGB"]

FileHeader
~~~~~~~~~~

======================  =========  ======================================================
field                   type       meaning
======================  =========  ======================================================
magic                   char[4]    ``MIGB``
format_version          u16        ``2``; a reader takes 1 and 2 and refuses what it cannot read
umi_len                 u8         UMI length in bases, 0 if there is no UMI
cell_len                u8         cell barcode length, 0 if there is none
bucket_index            u8         which range partition this file is
bucket_bits             u8         number of key bits used to partition; 0 = one bucket
paired                  u8         1 if mate 2 is present
barcode_quality         u8         v2: 1 if every record carries the barcode's own quality
sample_id               str        length-prefixed (u32 + bytes)
provenance              str        length-prefixed JSON: command line, version, pattern
quality_calibration     f32[]      length-prefixed; measured error rate per reported Phred
======================  =========  ======================================================

``quality_calibration`` being empty means "not measured, fall back to :math:`10^{-q/10}`". It is
carried in the file rather than recomputed because it is estimated once, by ``checkout``, from
mismatches against the constant segments of the barcode pattern — and on a 2-colour instrument
that emits only four distinct quality values, the nominal Phred is wrong by an order of magnitude
and every downstream likelihood inherits the error.

Block
~~~~~

A block header in plaintext, then a compressed payload::

    n_records u32 | raw_bytes u32 | stored_bytes u32 | crc32 u32 | codec u8 | reserved u8[3]

``codec`` is 0 for stored and 1 for zlib deflate level 1. The CRC is over the *uncompressed*
payload. The payload is **column-major**:

1. ``n_records`` fixed records: ``cell u64, umi u64, src_index u64, flags u16, umi_minq u8,
   cell_minq u8, len1 u32, len2 u32``
2. all of ``seq1``, concatenated
3. all of ``seq2``
4. all of ``qual1``
5. all of ``qual2``
6. v2, and only when ``barcode_quality`` is set: all of the UMI's own quality, ``umi_len`` bytes
   per record
7. ...and all of the cell barcode's, ``cell_len`` bytes per record

Both barcode quality columns are fixed width, because both lengths are file constants — so they
cost no length field and the writer refuses a record that disagrees with the header rather than
shifting every column after it.

.. note::

   **v2 exists because v1 dropped the evidence ``refine`` needs.** v1 stored ``umi_minq`` and
   ``cell_minq``, the minimum over the barcode. The correction posterior does not want a minimum:
   it weighs the reported quality **at the position that differs**, and a minimum says every
   position is as bad as the worst one. That overstates the error everywhere, which makes merges
   easier — the wrong direction, since a wrong merge destroys a molecule and a missed one only
   inflates a count. A v1 file still reads: the quality comes back empty and ``refine`` falls back
   to the library's global rate, exactly as it does for a FASTQ with no ``QX:Z:`` tag.

Three decisions worth knowing, because they look wrong until you measure them:

**Sequence is raw ASCII, not 2-bit packed.** Packing saves 0.75 bytes per base, but the quality
string is the same length and is near-incompressible, so packing touches only about an eighth of
the record — and it destroys the cross-read redundancy that a compressor finds in amplicon data.
Measured on a 2×150 amplicon block: 197 B/pair raw+deflate versus 227 B/pair packed+deflate.
Packing came out *worse*.

**Column-major, not per-record.** Sequence and quality have very different symbol distributions;
interleaving them costs the compressor 10–20% on the same data.

**``src_index`` is a u64.** It is the sort tiebreak, so it is what makes output byte-identical at
one thread and at eight. A u32 caps at 4.29·10⁹ read pairs, which a NovaSeq X run exceeds, and on
overflow the guarantee fails silently and nondeterministically.

Buckets and ordering
~~~~~~~~~~~~~~~~~~~~

Files are **range** partitions of the sort key, not hash partitions: bucket
:math:`b = \mathrm{key} \gg (64 - \mathrm{bucket\_bits})`. Two consequences, both load-bearing:

* a barcode and its 1-mismatch neighbours mostly land in the same bucket, so correction can be
  applied locally. A hash sends them to uncorrelated files and permanently splits the molecule.
* bucket order *is* key order, so the on-disk sort by sample/cell/UMI is a property of the layout
  rather than a separate pass over the data.

Barcodes are 2-bit packed with **base 0 in the high bits**, so that the packed integer order
equals the lexicographic order of the barcode string. An ``N`` is stored as ``A`` with
``kUmiHasN``/``kCellHasN`` set, keeping the key a plain integer and the ambiguity out of band.

Flags describe what has **already been applied**, never what remains to be done — in particular
``kRevComp1``/``kRevComp2`` mean the stored mate is already reverse-complemented, so ``assemble``
must never re-orient anything.

Consensus FASTQ
---------------

The pipeline output, and the contract with everything downstream::

    @<sample>[.<cell>].<umi>[.c<k>][.<m>] RX:Z:<umi>\tBC:Z:<sample>\tCB:Z:<cell>\tMI:Z:<name>\tcD:i:<reads>

Every field is separated by a **dot**, and there is no colon anywhere in the name. ``.c<k>`` is the
overlap component in ``--contig`` mode and ``.<m>`` the linkage split index; each is written only
when there is more than one, so the name carries three, four or five fields. A parser must
therefore read it **from the right** rather than indexing from the left -- and a sample id may
itself contain a dot, since ``validate_sample_id`` does not forbid one.

Tags are separated by **TAB**, not space: ``bwa mem -C`` and ``minimap2 -y`` append the FASTQ
comment verbatim into the SAM record, so it has to be SAM-conformant or the resulting BAM is
malformed. The UMI comes last in the read name because that is the convention ``fgbio``'s
``CopyUmiFromReadName`` and ``umi_tools`` both assume.

.. warning::

   ``dnaio`` — used by arda's rnaseq module — **drops the comment entirely**. Anything a
   downstream Python tool must see has to be in the read *name*, which is why the name is
   self-sufficient rather than a bare integer.

QC tables
---------

Every stage writes plain TSV beside its output. These are the contract with :doc:`plot <plots>` and
with any script of your own -- a figure is never computed from the FASTQ, only from one of these,
which is what stops a figure and a report disagreeing.

.. list-table::
   :header-rows: 1
   :widths: 32 16 52

   * - file
     - written by
     - columns
   * - ``checkout.summary.tsv``
     - ``checkout``
     - one row per **sample**: ``reads``, ``umis``, effective length, effective space, error rate,
       saturation. The per-sample UMI count the ``sample_umis`` panel draws.
   * - ``checkout.coverage.tsv``
     - ``checkout``
     - ``sample_id``, ``mig_size`` (power-of-two bin start), ``reads``, ``units``
   * - ``<sample>.sizes.tsv``
     - ``refine``
     - ``size``, ``log1p_size``, ``molecules``, ``reads`` -- the **exact** MIG size spectrum, one
       row per distinct depth
   * - ``<sample>.rank.tsv``
     - ``refine``
     - ``rank``, ``reads``, ``cumulative_reads``, ``cumulative_fraction``; molecules, log-spaced
   * - ``<sample>.cell_rank.tsv``
     - ``refine``
     - ``rank``, ``umis``, ``called``, ``cumulative_umis``, ``cumulative_fraction``; cell barcodes
       sorted by **distinct UMIs**, log-spaced. Written only when the reads carry a cell barcode.
   * - ``<sample>.cells.tsv``
     - ``refine``
     - ``cell``, ``molecules``, ``called`` -- one row per cell barcode, unsorted
   * - ``<sample>.bins.tsv``
     - ``refine``
     - per size bin: barcodes, reads, merged, erroneous fraction, molecules, residual FDR, payload
       entropy
   * - ``<sample>.barcodes.tsv``
     - ``refine``
     - ``cell``, ``umi``, ``reads``, ``corrected_reads``, ``parent`` -- the whole barcode table
   * - ``assemble.coverage.tsv``
     - ``assemble``
     - ``sample_id``, ``min_reads``, ``max_reads``, ``groups``
   * - ``assemble.quality_by_depth.tsv``
     - ``assemble``
     - per power-of-two depth bin: ``molecules``, ``q_min``, ``q_p25``, ``q_median``, ``q_p75``,
       ``q_max``, ``q_mean``
   * - ``<sample>.mig.tsv``
     - ``assemble``
     - one row per molecule: barcode, contigs, reads, support, length, quality, error, linkage

.. note::

   ``<sample>.sizes.tsv`` is at **exact** sizes, not power-of-two bins, and that is deliberate: the
   rank/Zipf curve is its cumulative count, and four bins make four steps. It costs one row per
   distinct depth -- a few thousand on a real library -- rather than one row per molecule.

.. note::

   ``<sample>.cell_rank.tsv`` and ``<sample>.rank.tsv`` are both **log-spaced**: consecutive rows
   step by about 5% of the rank, with the first and last always emitted so the ends of the curve
   are exact. One row per barcode would be hundreds of millions of rows for a figure that is read
   on a log axis anyway.

``assemble.quality_by_depth.tsv`` holds order statistics rather than a sample. ``assemble``
accumulates the exact joint distribution of (depth bin, rounded Phred) -- both are small integers,
so it is 61 counters per bin -- and the quantiles are read off that. Nothing is thinned, which is
why the quality panel can be a box instead of a scatter.
