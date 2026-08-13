Layouts: where the barcode is
=============================

There is exactly one thing migec has to be told, and everything after it -- correction, consensus,
the quality cap -- is the same three commands whatever the answer. Four ways to say it, in the
order you should reach for them.

1. A position
-------------

Most libraries put the barcode at a fixed offset in one read. That is the primary mode, and it
needs neither a sample sheet nor an anchor:

.. code-block:: bash

    migec checkout reads.fq.gz --bc-pattern '^NNNNNNNN' -o out/     # 8 nt UMI at the read start
    migec checkout reads.fq.gz --bc-pattern '0:8'       -o out/     # the same, as a slice

Two spellings, one meaning:

.. list-table::
   :header-rows: 1
   :widths: 30 30 40

   * - pattern
     - slice
     - what it is
   * - ``^NNNNNNNN``
     - ``0:8``
     - an 8 nt UMI at the first base
   * - ``^NNNN.NNNNN``
     - ``0:4,5:10``
     - a 9 nt UMI split by one skipped base
   * - ``^XXXXXXXXXXXXXXXXNNNNNNNNNN``
     - ``cell:0:16,16:26``
     - a 16 nt cell barcode then a 10 nt UMI (10x)
   * - ``^NNNNXNNN``
     - ``0:4,cell:4:5,5:8``
     - a UMI interrupted by one cell-barcode base

In a pattern, ``N`` is a UMI base, ``X`` a cell-barcode base, ``.`` a base that is skipped --
neither scored nor captured -- and anything else (``ACGT`` or any IUPAC symbol) is constant
sequence that gets scored. Slices are **half-open and 0-based**, like Python's: ``0:8`` is eight
bases and the next slice may start at 8. Each is a UMI slice unless it is prefixed ``cell:``.
Slices must be in increasing order and must not overlap, because one base belongs to one barcode.

Anchoring
~~~~~~~~~

A leading ``^`` says the barcode starts at the first base. Every slice list says the same thing by
construction, since a position is only a position if it is measured from somewhere. Both set
``--max-offset 0``, and a layout with nothing to score is anchored automatically even without the
caret -- which is why the flag no longer appears in any of the examples.

Never: this is not a convenience. Placement is a hypothesis test, and a pattern with no constant
sequence supplies no evidence for it. Asked to scan freely, ``compile()`` refuses rather than
picking an offset; asked to scan a 5 nt dual-end handle, it refuses too, because ``TGACT`` occurs
by chance about every kilobase and the bar for a free scan is ``log2(offsets/alpha)`` bits, which
five bases cannot pay. Anchored, there is only one place to be and the bar does not apply. See
:doc:`checkout` for the arithmetic.

2. A named preset
-----------------

.. code-block:: bash

    migec sheet --presets                                       # all of them, and their sources
    migec checkout R1.fq.gz R2.fq.gz --preset 10x-v2 -o out/

.. list-table::
   :header-rows: 1
   :widths: 14 34 52

   * - preset
     - layout
     - what it is, and where the layout is written down
   * - ``umi``
     - ``^NNNNNNNN``
     - generic inline UMI. Change the run length, or write the slice.
   * - ``migec``
     - ``cagtggtatcaacgcagagtNNNNtNNNNtNNNN``
     - MIGEC 5'-RACE RepSeq: the SMART adapter then a 12 nt UMI split by two spacers.
       ``misc/barcodes.txt`` of MIGEC 1.2.9, tag ``v1-final``. Prefix a sample tag per row to
       demultiplex.
   * - ``primerid``
     - ``NNNNNNNNNcagtttaacttttgggccatcca``
     - HIV-1 Primer ID amplicon as used by MAGERI. Recovered by ``migec suggest`` from
       ``SRR1763769``; the primer places it, so the scan stays free.
   * - ``duplex``
     - ``^NNNNNNNNNNNN.....`` on both mates
     - duplex sequencing: a 12 nt UMI and a 5 nt spacer per mate, 24 nt together.
   * - ``10x``
     - ``^XXXXXXXXXXXXXXXXNNNNNNNNNNNN``
     - 10x Chromium 3' v3/v3.1: 16 nt cell barcode, 12 nt UMI, on R1.
   * - ``10x-v2``
     - ``^XXXXXXXXXXXXXXXXNNNNNNNNNN``
     - 10x Chromium 3' v2 and 5' v1/v2: 16 nt cell barcode, 10 nt UMI.
   * - ``tso500``
     - ``^NNNNN.....`` on R1 only
     - Illumina TSO500 ctDNA. The fgbio read structure is ``5M5S+T +T`` -- R2 is all template.
       See the warning below.
   * - ``smarter-umi``
     - ``^NNNNNNNNNNGGG``
     - SMARTer template-switching RNA-seq: a 10 nt inline UMI, then the ``GGG`` the template switch
       leaves behind. ``ncgr/UMI-analysis``, whose quality filter reads offset 0 length 10.

Warning: the ``duplex`` preset extracts the tags and emits **single-strand** consensuses. Pairing
the two strands of a molecule into a duplex consensus is not implemented, so no duplex error rate
should be quoted from this output.

Never: **a 5 nt UMI does not identify a molecule, and TSO500's does not pretend to.** 4^5 is 1,024
barcodes against the tens of thousands of fragments a ctDNA panel region carries, so the space is
saturated by construction and the birthday bound says most barcodes are shared. TSO500's own
pipeline resolves that by grouping on the UMI **and the mapping position** (``fgbio
GroupReadsByUmi``, which runs after alignment); migec groups on the barcode, before any alignment
exists, so it cannot. It will report the space as saturated, set ``err_unreliable``, and warn --
and on this chemistry that warning is the correct answer, not a threshold to raise. Use migec here
to extract and tag; do the grouping position-aware, downstream.

3. A read structure
-------------------

fgbio, Picard, samtools and the TSO500 pipelines describe a layout as a *read structure*, and
migec takes them verbatim: ``M`` a molecular barcode, ``B`` a sample/cell barcode, ``S`` a skip,
``T`` template.

.. code-block:: bash

    migec checkout R1.fq.gz R2.fq.gz --read-structure 5M5S+T -o out/    # TSO500: `5M5S+T +T`
    migec checkout R1.fq.gz R2.fq.gz --read-structure 12M5S+T --read-structure2 12M5S+T -o out/

.. list-table::
   :header-rows: 1
   :widths: 20 44 36

   * - structure
     - pattern
     - platform
   * - ``5M5S+T``
     - ``NNNNN.....``
     - TSO500
   * - ``16B10M+T``
     - ``XXXXXXXXXXXXXXXXNNNNNNNNNN``
     - 10x 5'
   * - ``8M+T``
     - ``NNNNNNNN``
     - a plain inline UMI

The pattern stops at the first template segment, because everything after it is payload and migec
trims to exactly that point. ``+`` means "the rest of the read" and is valid only on the last
segment; ``+M`` is refused, since an unbounded barcode has no length for the collision arithmetic
to use.

A read structure is positional by definition, so it carries its own anchor.

4. A barcode table
------------------

For many samples in one file. This is MIGEC's own ``barcodes.txt``, read verbatim, which is why
the published tables run unchanged:

.. code-block:: text

    S1	aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
    S2	aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN

Uppercase is matched exactly, lowercase is the fuzzy adapter (scored at half weight), and the
``t``\ s between UMI runs are pattern bases, not barcode -- the UMI is 12 nt, so the barcode space
is ``4^12`` and not ``4^14``. Column 3 is the *slave* pattern, on the other mate, whose captured
positions **extend** the UMI:

.. code-block:: text

    S1	NNNNNNNNNNNNtgact	agtcaNNNNNNNNNNNN

Never: both halves must match or the read is unmatched. Accepting the master alone would emit 12 nt
UMIs beside 24 nt ones, and every collision estimate downstream would then be computed over two
barcode spaces at once.

Rows may share a sample id -- that is how a sample sequenced with more than one tag is declared,
and one output file is written per *sample id*, never per row.

.. code-block:: bash

    migec sheet barcodes.txt       # what will each row extract, before anything runs

If you do not know the layout
-----------------------------

Do not guess:

.. code-block:: bash

    migec suggest reads.fq.gz

It segments the per-cycle base composition into UMI, constant and payload runs and prints a
paste-ready pattern. It recovered the 9 nt + ``CAGTTTAACTTTTGGGCCAT`` layout of ``SRR1763769``
unaided. Note: it stops at the last *constant* run -- composition alone cannot tell a UMI from
diverse payload, only the anchor can, and it says so when that is all it found. See :doc:`suggest`.
