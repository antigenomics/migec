Bring your own UMI
==================

``migec checkout`` is one way of producing migec's input. It is not the only one, and for a capture
assay it is usually the wrong one: on an exome, a hybrid-capture panel, a ctDNA or an MRD library
the UMI is in the **index read**, so it never appears inside R1 or R2 and there is no pattern for
``checkout`` to find. What the sequencing facility hands over is a BAM from fgbio, Picard or the
kit vendor's pipeline, with the UMI already in the ``RX`` tag.

That file is a migec input. ``refine``, ``assemble`` and ``subsample`` group on ``RX``, and where
that tag came from is not their business.

The contract
------------

A record is migec input when it carries these tags. Only ``RX`` is required.

.. list-table::
   :header-rows: 1
   :widths: 10 30 60

   * - tag
     - what it is
     - if it is missing
   * - ``RX``
     - the UMI
     - ``refine`` and ``assemble`` refuse the file by name. Nothing to group on
   * - ``QX``
     - the UMI's own base qualities, one Phred per base
     - the correction posterior falls back to the library's global error rate, exactly as it does
       for a v1 ``.mig`` file. Correction still runs; it just has one fewer piece of evidence
   * - ``CB``, ``CY``
     - cell barcode and its qualities
     - the library is treated as one cell, which is what a bulk library is
   * - ``BC``
     - sample id
     - the sample is named ``sample``, or whatever ``--sample`` says

In a FASTQ those tags live in the **comment** — everything after the first space or tab in the
header — separated by tabs, which is what makes the line a valid SAM record once an aligner appends
it (:doc:`formats`). In a BAM they are ordinary tags. The two are the same input.

From a BAM
----------

.. code-block:: bash

   migec refine   umi_tagged.bam -o ref/
   migec assemble ref/S1.fq.gz   -o asm/

BAM, SAM and CRAM are recognised from the file itself, not from the name — there is no flag. Behind
the scenes migec converts once with ``samtools`` into a temporary FASTQ inside the output
directory, which is deleted when the stage returns:

.. code-block:: bash

   samtools collate -u -O in.bam | \
     samtools fastq -n -T RX,QX,CB,CY,BC -1 R1.fq -2 R2.fq -0 R0.fq -s S.fq -

Two things about that command are load-bearing.

.. warning::

   **Collate is not optional on an aligned file.** ``samtools fastq -1/-2`` pairs by *adjacency*,
   so on a coordinate-sorted BAM it writes one molecule's mate 1 beside another molecule's mate 2.
   ``assemble`` matches mates strictly by position, so the wrong pair would be consensed and
   nothing downstream could tell. An **unaligned** BAM cannot be coordinate-sorted, so migec skips
   collate there and the record order is preserved; ``@SQ`` in the header is what decides, and a
   header that cannot be read counts as aligned.

.. note::

   The temporary FASTQ is roughly **4x the BAM** on disk while the stage runs, and it is written
   next to your output, not in ``/tmp`` — the same siting as the range-partition buckets, for the
   same reason.

Getting the UMI into ``RX`` in the first place
----------------------------------------------

If your BAM does not have it yet, one of these does it. All three are the vendor's own tool; migec
does not reimplement any of them.

.. code-block:: bash

   # Illumina-style: R1, R2 and a separate UMI read
   fgbio FastqToBam -i R1.fq.gz R2.fq.gz UMI.fq.gz --read-structures +T +T +M -o tagged.bam

   # a BAM plus the UMI read
   fgbio AnnotateBamWithUmis -i in.bam -f UMI.fq.gz -o tagged.bam

   # Picard, from FASTQ
   picard FastqToSam F1=R1.fq.gz F2=R2.fq.gz UMI_FASTQ=UMI.fq.gz ... O=tagged.bam

Then run ``migec refine`` on ``tagged.bam``.

.. warning::

   ``umi_tools extract`` puts the UMI in the **read name**, not in a tag, and migec does not parse
   read names — the name is only ever copied through. Convert it, or use a tool that writes ``RX``.
   This is deliberate: a name is free-form and a tag is not, so guessing at a name is how one
   pipeline's ``_`` separator becomes another's silent truncation.

From a FASTQ somebody else tagged
---------------------------------

Nothing special is needed. If the comment carries ``RX:Z:`` after a tab, it works:

.. code-block:: text

   @A00123:45:HXXX:1:1101:1000:1000 RX:Z:ACGTACGT<TAB>QX:Z:IIIIIIII
   TATCAGAGTAGTGGTATTTCACAGGCGGCCAGCAGGG
   +
   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

``samtools fastq -T RX,QX`` writes exactly this, so a round trip through a BAM and back is a no-op
as far as migec is concerned.

What this does not change
-------------------------

Grouping is still on **sample + cell + UMI** and never on the alignment position. That is the
difference between migec and a map-first tool such as fgbio's ``GroupReadsByUmi`` or
``umi_tools group``, and it is measured rather than asserted — see :doc:`grouping`. Reading their
file format is not the same as adopting their model: a tool that *transports* the barcode composes
with migec, and a tool that *deduplicates* on it replaces a stage of migec
(:doc:`downstream`).
