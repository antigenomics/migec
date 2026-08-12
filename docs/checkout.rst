Checkout — barcode extraction
=============================

``migec checkout`` finds the barcode pattern in each read, extracts the sample tag and the UMI,
trims the synthetic sequence away, and puts the barcode in the read header.

.. code-block:: bash

   migec checkout reads.fq.gz --barcodes barcodes.txt --out out/

Barcode patterns
----------------

The pattern grammar is MIGEC's, so published barcode tables work unchanged:

.. code-block:: text

   S1	aaACTcagtggtatcaacgcagagtNNNNtNNNNtNNNN
   S2	aaAGAcagtggtatcaacgcagagtNNNNtNNNNtNNNN

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - symbol
     - meaning
   * - ``ACGT``
     - a scored position, matched exactly
   * - IUPAC (``R``, ``Y``, …)
     - a scored position with a degenerate set: ``R`` = A|G, ``Y`` = C|T
   * - lowercase
     - a scored position at **half weight** — the adapter, where a mismatch is expected
   * - ``N`` or ``n``
     - a UMI position: captured, never scored
   * - ``.``
     - a wildcard: neither scored nor captured

``N`` always means UMI here, never IUPAC "any base" — use ``.`` for an uncaptured wildcard.

UMI runs need not be contiguous. ``NNNNtNNNNtNNNN`` captures twelve bases into one UMI, with the
lowercase ``t`` spacers scored at half weight — which is how the published MIGEC tables are
written.

``migec sheet barcodes.txt`` prints what each row will extract without running anything.

How a match is accepted
-----------------------

Not by counting mismatches. For each scored position, with IUPAC set :math:`S` of size :math:`m`,
observed base :math:`b` and error probability :math:`e`:

.. math::

   s_i = \begin{cases}
     \log_2 4\left[\frac{1-e}{m} + \frac{(m-1)e}{3m}\right] & b \in S \\
     w \cdot \log_2 \frac{4e}{3} & b \notin S
   \end{cases}

which is the log-likelihood ratio of "the tag is here" against "this is random sequence".
:math:`w` is 1.0 for uppercase and 0.5 for lowercase. At :math:`m = 1` a match is worth **+2.00
bits**, a mismatch **−9.55 bits at Q30** and **−0.60 bits at Q2**.

So a mismatch on a bad base is nearly free and a mismatch on a good base is fatal. That is what
MIGEC's "good mismatch / bad mismatch" counting was reaching for, done continuously — and without
its two defects, both of which are reproduced as tests here: v1 read the quality string from the
*start of the read* rather than from the match offset, and a dangling ``else`` meant low-quality
mismatches were never counted at all.

The acceptance threshold defaults to a Bonferroni bound over the offsets actually scanned,
:math:`\log_2(n_\text{offsets} \cdot n_\text{patterns} / \alpha)`.

.. note::

   Reads are not i.i.d. uniform ACGT — shared primers and composition bias violate the null badly.
   Treat the default threshold as a starting point and calibrate it against shuffled decoy
   patterns on your own data.

Ambiguous is not the same as unmatched
--------------------------------------

A read whose best sample tag does not beat the runner-up by ``--min-margin`` bits is reported as
**ambiguous** rather than assigned to one of them. That is a different diagnosis from *unmatched*:
ambiguous means the barcodes in the sheet are too close together, unmatched means the pattern is
wrong or absent. One counter cannot say both, so there are two.

Trimming
--------

``--trim pattern`` (the default) drops everything up to and including the matched pattern — the
adapter, the sample tag and the UMI. What is left is exactly the biological payload. This is what
you want before alignment: the tag is synthetic and will soft-clip at best, mismap at worst.

``--trim none`` keeps the read whole and puts the UMI in the header only.

Barcode transfer to the header
------------------------------

The barcode travels on in SAM-style tags::

   @r0 RX:Z:GCTAAAGACAAT	QX:Z:IIIIIIIIIIII	BC:Z:S1
   TACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTT
   +
   IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

``RX``/``QX`` are the UMI and its qualities (read by fgbio, Picard and umi_tools); ``BC`` is the
sample barcode. Tags are separated by **TAB**, with one space after the read name.

.. warning::

   The TABs are not cosmetic. ``bwa mem -C`` and ``minimap2 -y`` copy the FASTQ comment verbatim
   into the SAM record, so it has to be SAM-conformant there or the resulting BAM is malformed.

.. warning::

   ``dnaio`` — used by arda's rnaseq module — **drops the comment entirely**. Anything a
   downstream Python tool must see has to be in the read *name*.

Quality filtering of the UMI
----------------------------

``--min-umi-quality`` defaults to **0**, i.e. nothing is dropped. MIGEC used 15 and MAGERI 20, both
as hard drops. The default here is deliberately different: a low-quality UMI base is a reason to be
less certain about which molecule a read belongs to, not a reason to throw the read away. The
correction step can usually recover it, and a molecule seen three times is information.

Output
------

===================================  =========================================================
file                                 content
===================================  =========================================================
``<sample>.fq.gz``                   trimmed reads with barcodes in the header
``unmatched.fq.gz``                  reads matching no pattern (``--write-unmatched``)
``checkout.summary.tsv``             one row per sample: yields, UMI statistics, correction
``checkout.coverage.tsv``            reads and distinct UMIs per power-of-two MIG size
``checkout.umi_composition.tsv``     per-position base usage, entropy, information, collision
``checkout.json``                    everything above, machine-readable
===================================  =========================================================

See :doc:`umi_statistics` for what the last two contain and how to read them.
