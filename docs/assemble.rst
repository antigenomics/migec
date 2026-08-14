assemble -- one consensus per molecule
======================================

Reads carrying the same barcode are reads of one molecule, so collapse them into one consensus.

.. code-block:: bash

   migec assemble out/S1.fq.gz -o cons/
   migec assemble out/S1.fq.gz -o cons/ --contig      # random-primed reads that tile a molecule
   migec assemble out/S1.fq.gz -o cons/ --rt-error 3e-5

The input is a per-sample FASTQ written by :doc:`checkout`, with the barcode in the ``RX``/``CB``
tags. The output is ordinary FASTQ, one record per molecule, plus a per-molecule table.

The barcode is the whole key
----------------------------

A molecule is identified by **sample, cell barcode and UMI together**, never by the UMI alone. The
same UMI turning up in two cells or two samples is the normal case, not an error: a UMI is only
ever unique inside the compartment it was added to, and 4\ :sup:`12` random tags reused across ten
thousand cells is the design, not a defect.

So the sort key is ``(cell, umi, src_index)``, the sample is the file, and the partition is on the
**cell** whenever there is one — which is also what makes a per-cell scope contiguous on disk.

.. code-block:: text

   @S1.AAAACCCCGGGGTTTT.ACGTACGTACGT RX:Z:ACGTACGTACGT	BC:Z:S1	CB:Z:AAAACCCCGGGGTTTT	MI:Z:S1.AAAACCCCGGGGTTTT.ACGTACGTACGT	cD:i:12

Nothing scales with the library
-------------------------------

Grouping needs a molecule's reads together, and there is no hash map keyed by barcode anywhere in
this pipeline — at NovaSeq scale that is 19 GB. Reads are **range partitioned** on the packed key
into ``.mig`` buckets, and one bucket is sorted in RAM at a time.

Range, never hash. A hash sends a barcode and its 1-substitution neighbours to uncorrelated
buckets, which makes UMI correction impossible to apply locally and splits the molecule
permanently — and each half looks like a well-formed MIG, so nothing detects it. Range costs the
same, because a packed barcode is close to uniform, and it has the property a hash cannot: bucket
order *is* key order, so the output comes out sorted for free.

The bucket count comes from the input size, not from a flag, and the writer buffer budget is
**split across** the open buckets rather than being per-bucket — otherwise cutting the input
finer, which exists to use less memory, would use more.

.. code-block:: text

   2,423,777 reads/s = 3,741,782 partitioning + 1,713,585 groups/s consensus

    1 bucket:     203 MB
   16 buckets:    123 MB

``tests/benchmark/test_assemble_speed.py`` asserts that raising the bucket count lowers peak RSS
and does not change the output. Each configuration is measured in its own process, because
``peak_rss_bytes`` is a process high-water mark and two runs in one interpreter cannot be compared.

...or the partition arrives already built. ``migec checkout --mig`` writes its reads into exactly
these buckets, on exactly this key, so pointing ``assemble`` at them skips the pass above
entirely:

.. code-block:: bash

   migec assemble out/S1.000.mig -o asm     # one bucket names the whole partition

The layout is read out of the file headers rather than chosen — the buckets *are* the partition,
and re-deriving the bucket count from the input size would address them wrongly. Two things are
refused rather than guessed at: buckets from **two samples** (assemble is a per-sample stage, and
a UMI repeats across samples by design), and ``--limit-read``/``--limit-umi``, because a limit is
a prefix of the input and a partition has no prefix left — the first records of bucket 0 are one
corner of the barcode space. Limit at ``checkout`` instead.

Never: only buckets this stage wrote are deleted after they are consumed. In ``--mig`` mode they
are checkout's output, and a second ``assemble`` over an eaten partition finds a fraction of the
library and reports it as a smaller one.

The consensus and its quality
-----------------------------

Per column, over the reads that reach it:

.. math::

   LL[j][b] = \sum_i \left[ r_{ij} = b \right] \log(1 - e_{ij}) +
              \left[ r_{ij} \neq b \right] \log(e_{ij}/3)

``log(1-e)`` and ``log(e/3)`` depend only on the reported Phred, so they tabulate; only the
posterior needs a transcendental, three times per *output* base rather than once per read base.

Then the part that matters:

.. math::

   Q(j) = -10 \log_{10}\left( p_\text{cons}(j) + p_\text{floor} \right)

The floor is **added, not compared**. An error made during reverse transcription or the first PCR
cycle is in every read of the molecule and no consensus removes it, so the two failure modes are
independent and the emitted quality has to carry both. :doc:`X2 measured the floor <quality_floor>`
at 1.54·10\ :sup:`-4` on an HIV-1 Primer ID control, which caps every emitted quality at about
**Q40**, which is also 10x's stated figure for the V(D)J RT. A blanket 1e-6 is excluded for an
RT protocol by two orders of magnitude -- but it is the right floor for a DNA workflow with a
proofreading polymerase, which is why ``--rt-error`` names the chemistry: ``rt`` (1e-4, the
default), ``medium`` (1e-5), ``high`` (1e-6), or the rate itself. ``--pre-amp-error auto``
:ref:`fits it from the dataset <pre-amp-auto>` when the library is clonal and deep enough to carry
the measurement, and says which class it fell back to when it is not.

A tie is resolved by base order rather than by an ``N``: the posterior is then 0.5 and the emitted
quality says so at about Q3. An ``N`` would discard the information that it is one of two.

.. warning::

   The floor is the **one-molecule** floor, and every record migec emits is one molecule. 10x put
   it exactly: bases covered by a single UMI are assigned Q40, and only bases covered by at least
   two UMIs are assigned Q60. The difference is not a better consensus — it is a second molecule.
   An RT error is common-mode within a molecule and independent between them, so the floor divides
   out only when independent molecules agree, and combining molecules is
   `arda <https://github.com/antigenomics/arda>`_'s job. A per-molecule record claiming Q50 would
   be claiming two-UMI confidence on one-UMI evidence.

Coverage is capped
------------------

At most **10,000 reads** of a barcode enter the consensus. The 10,001st moves no call and no
emitted quality — the posterior saturated thousands of reads earlier and the quality is capped by
the floor regardless — while the group still costs time and resident memory in proportion to its
size. 10x cap the same way for the same reason: *"Very high coverage (greater than 10,000 reads)
of transcripts can be problematic because it degrades computational performance and adds little
information."*

.. warning::

   The cap is on the reads that are **consensed**, never on the reads that are **counted**. The
   ``cD`` tag and the ``reads`` column stay the molecule's true depth. Capping a count would
   flatten the abundance of exactly the most-amplified molecules, which is the measurement a UMI
   pipeline exists to make.

Counting mode: ``--fast``
-------------------------

.. code-block:: bash

   migec assemble ref/S1.fq.gz -o cons/ --fast

Emit the group's most frequent **exact** sequence, with every base carrying the best quality any
read *of that sequence* reported for it. There is no column model, so there is no per-base error
correction and no sub-clustering; the RT floor still caps what is claimed.

Use it when the deliverable is a molecule **count** — expression, clonotype abundance — and the
sequence only has to be right often enough to assign. Measured against the full path on 8-read
molecules at 5·10\ :sup:`-3` per base, the column posterior clears essentially every sequencing
error and the majority string keeps whatever it carried.

.. note::

   The per-base maximum is taken over the reads carrying the **winning** sequence, not over every
   read in the group. Maximising across variants would take its highest quality from exactly the
   reads that disagree at that position, asserting most confidence where the evidence conflicts.

   ``--fast`` is refused with ``--contig``: tiling reads share no exact sequence, so a majority
   vote over whole strings returns one fragment and drops the rest of the molecule.

``support`` in ``<sample>.mig.tsv`` is how many of the group's reads carried what was emitted, and
the report prints the same thing as a share. Well below 1 means the reads of a molecule disagree —
which the full path resolves per column and this one resolves by majority vote over whole strings.

Splitting a group into two molecules
------------------------------------

The discriminator is **linkage**, not a count of polymorphic sites. Independent PCR subclones
almost never co-segregate; a real second molecule does, on the same reads.

The threshold is **8.68**, a Bonferroni'd ``-log10 p`` over pairs of callable positions, and it is
a :doc:`measured false-positive point <nulls>` rather than a derivation. The nominal ``p < 0.01``
calls 30.62% of MIGs against 1.60% — a 19× over-call, because reads are not exchangeable and a
low-quality read carries a minor base at many positions at once.

The test is two-sided: at a 50/50 split which allele is "major" is a coin toss taken separately
per column, so a genuine doublet's columns come out anti-correlated as often as not.

.. note::

   The threshold implies a **minimum group size**. The strongest evidence a pair of columns can
   carry is ``log10 C(n, n/2)``, so a 50/50 split needs about 34 reads before it can clear 8.68 at
   all. Below that the data cannot separate a subclone from two bad reads at a 1% false-positive
   rate, and migec does not pretend otherwise.

Contig assembly: ``--contig``
-----------------------------

Random priming does not give co-terminal reads. Reads sharing a barcode tile the molecule at
different starts, and :doc:`X1 measured <fragmented>` this on 10x 3' GEX: the co-terminal
assumption is false for 92% of groups, and **27.3% of groups hold more than one overlap
component**.

``--contig`` places the reads against each other by exact seed matching, cuts them into overlap
components with a union-find that carries each read's offset, and emits **one consensus per
component**. A component is never extended across a gap — two reads that share a barcode but no
sequence are two contigs of the molecule, and a single consensus over them would assert sequence
that no read covers.

This is contig assembly of *one molecule's* fragments, and that is all it is. Assembling a
cell's full-length receptor, calling doublets and filtering contaminating chains are
`arda <https://github.com/antigenomics/arda>`_'s job, downstream of here.

.. warning::

   **Contig assembly needs a barcode that is not saturated.** Two fragments of two *different*
   molecules that happen to share a barcode have no sequence in common — which is exactly what two
   fragments of one molecule look like. There is nothing in the data that separates them.

   ``assemble`` runs the same birthday arithmetic checkout does, on the barcodes this run actually
   saw, and reports ``expected_molecules_per_group`` = ``E[k | k ≥ 1]`` for the Poisson-occupied
   space. When a short UMI cannot tag every input molecule distinctly *by design*, that number is
   above 1 and says how far, and contig mode warns when more than 5% of groups hold more than one
   molecule.

Merging the mates: ``--mate2``
------------------------------

Both ends of the same fragment carry the same barcode, so a pair is two fragments of one molecule
and the same rule applies:

.. code-block:: bash

   migec assemble out/S1_R1.fq.gz --mate2 out/S1_R2.fq.gz -o cons/
   migec assemble out/S1.000.mig --merge-mates -o cons/    # the pair is already in the record

Mate 2 is reverse-complemented and **placed** against mate 1. When the mates overlap the molecule
comes back as one consensus spanning the insert; when they do not it comes back as two contigs,
because the bases between the mates are covered by no read. This is ``--contig``'s rule, not a
second one.

The offset is a property of the **molecule**, not of the pair — every pair of a group is a copy of
one fragment — so it is voted on once per group over up to eight pairs and then applied to all of
them. The first version instead ran the group through ``--contig``'s all-against-all placement and
cost **11x** the single-end path (119,820 record-pairs/s against 1,356,819 records/s); voting once
brings it to 1,288,686 record-pairs/s against 2,115,912, which is the cost of the second read's
bases and nothing else. ``tests/benchmark/test_assemble_speed.py`` holds the floor.

.. note::

   Mate 2 supplies **sequence only**. The barcode, the sample and the read name come from mate 1,
   which is the mate ``checkout`` anchored its pattern on, and the two files are matched by
   position — the contract checkout's ``_R1``/``_R2`` output keeps. A file that runs out before
   the other is refused rather than paired off by position, which would attach one molecule's mate
   to another's with a consensus that looks perfectly well formed.

Output
------

======================================= =======================================================
``<sample>.consensus.fq.gz``            one record per molecule, barcodes in the name and tags
``<sample>.mig.tsv``                    cell, umi, contig, molecule, reads, length, quality,
                                        consensus error, linkage score
``assemble.coverage.tsv``               groups per power-of-two MIG size
``assemble.json``                       all of it, machine-readable
======================================= =======================================================

The report prints the MIG size histogram rather than one error number, because a consensus over
one read *is* that read and averaging it in hides the thing being measured.
