Downstream: what consumes the consensus
=======================================

``migec assemble`` writes ordinary FASTQ. One record is one molecule, and the molecule's identity
is carried twice:

.. code-block:: text

    @PBMC.AAACCTGCAAAAGCAA.TTTGCCGATA RX:Z:TTTGCCGATA<TAB>BC:Z:PBMC<TAB>CB:Z:AAACCTGCAAAAGCAA<TAB>MI:Z:PBMC.AAACCTGCAAAAGCAA.TTTGCCGATA<TAB>cD:i:1
    TATCAGAGTAGTGGTATTTCACAGGCGGCCAGCAGGGCCGGCGGACCCCGCCCC
    +
    IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII

* the **name**, ``<sample>.<cell>.<umi>``, up to the first space;
* the **comment**, tab-separated SAM tags: ``RX`` UMI, ``QX`` its qualities, ``CB`` cell barcode,
  ``BC`` sample barcode, ``MI`` the molecule id, ``cD`` reads in the molecule.

The split matters, because the two halves survive different tools. A tool that keeps the comment
gets the tags; a tool that does not still gets a self-sufficient name. Nothing downstream needs
both.

Never: the tags are separated by **tabs**, not spaces. That is what makes the comment a valid SAM
record once an aligner appends it -- ``bwa mem -C`` and ``minimap2 -y`` copy the comment through
verbatim, and a space-separated comment would produce a SAM line that no parser accepts.

What was measured
-----------------

Every row below was run on this machine against a ``migec assemble`` output: 600 consensuses from a
synthetic 10x-shaped library (20 cells x 30 molecules x 3 reads), plus the 10x VDJ-T fixture from
``isalgo/umi_data`` for arda.

.. list-table::
   :header-rows: 1
   :widths: 14 30 26 30

   * - tool
     - command
     - what arrives
     - result
   * - minimap2 2.31
     - ``minimap2 -ax sr -y ref.fa cons.fq.gz``
     - name + all tags
     - 600/600 records carry ``RX``, ``CB``, ``MI``; ``samtools quickcheck`` valid
   * - bwa 0.7.19
     - ``bwa mem -C ref.fa cons.fq.gz``
     - name + all tags
     - 600/600 records carry ``RX``, ``CB``, ``MI``
   * - minibwa 0.7-r424
     - ``minibwa map -y ref.fa cons.fq.gz``
     - name + all tags
     - 600/600 records carry ``RX``, ``CB``, ``MI``, ``BC``, ``cD``; sorted BAM valid
   * - samtools 1.24
     - ``samtools sort``, ``samtools view``
     - name + all tags
     - tags round-trip through BAM unchanged
   * - arda 2.20.0
     - ``arda amplicon --r1 cons.fq.gz -p out``
     - name only
     - 385/2,651 mapped, TRA 120 / TRB 265, 7 clonotypes; ``sequence_id`` **is** the molecule id
   * - kallisto 0.52.0
     - ``kallisto quant --single``
     - sequence only
     - 171/600 pseudoaligned, ``est_counts`` are molecule counts
   * - salmon 2.5.1
     - ``salmon quant -l U -r``
     - sequence only
     - ``NumReads`` are molecule counts

The comment flag is per-tool and there is no majority convention:

.. list-table::
   :header-rows: 1
   :widths: 30 18 52

   * - invocation
     - flag
     - note
   * - ``minimap2 -ax sr``
     - ``-y``
     -
   * - ``bwa mem``
     - ``-C``
     -
   * - ``minibwa map``
     - ``-y``
     - the minimap2 spelling, not bwa's, even though minibwa is bwa-mem's successor
   * - ``minibwa mem``
     - ``-C``
     - the legacy CLI keeps bwa's spelling. lh3 marks it "not recommended"

Both wrong combinations exit non-zero with an error rather than dropping the tags quietly
(``minibwa map -C`` says ``unknown option``), so a mistake here costs a rerun and not a silent
loss of the molecule ids.

Note: **arda reads the name, not the comment.** It goes through ``dnaio``, which drops FASTQ
comments -- so the name has to be self-sufficient, and it is: ``sequence_id`` in the AIRR TSV comes
out as ``PBMC.AAACCTGCAGCCTGTT.CGTTTTTATC``, which is sample, cell and molecule without a join.
This is why the molecule id is in the name at all rather than only in ``MI:Z:``.

Not verified here
-----------------

* **STAR.** The contract is documented: STAR truncates the read name at the first whitespace and
  drops the comment, so the name survives and the tags do not -- the same position arda is in. It
  could not be confirmed on this machine: the Homebrew arm64 build of STAR 2.7.11b reports
  ``Number of input reads = 0`` for **any** FASTQ, including a one-record file with a plain name,
  so the failure is the build's and says nothing about migec's output.
* **bwa-mem2.** Same ``-C`` flag and the same comment-copying code path as ``bwa``, which was
  verified. No arm64 build to run.
* **NASC-seq2.** Not downstream of ``assemble`` at all -- it is an **alternative to** ``checkout``.
  Its demultiplexing is zUMIs, declared as ``UMI(12-19)`` / ``cDNA(23-200)`` / ``BC(201-220)`` with
  a ``find_pattern`` anchor, which is the same job migec's pattern does; its molecule tagging then
  runs on the aligned BAM. Feeding it a consensus FASTQ would collapse twice. :doc:`layouts` has
  the conversion, and the trap in it: zUMIs ranges are **1-based and inclusive**, migec slices are
  **0-based and half-open**, so ``UMI(12-19)`` is ``11:19`` and not ``12:19``.

Map first, or collapse first?
-----------------------------

There are two orders in use, and they are not interchangeable:

.. code-block:: text

   map first       raw reads --align--> BAM --group on (position, UMI)--> consensus
   collapse first  raw reads --group on (sample, cell, UMI)--> consensus --align--> BAM

migec is the second. `fgbio <https://fulcrumgenomics.github.io/fgbio/>`_, ``UMI-tools`` and
`UMIErrorCorrect <https://doi.org/10.1093/clinchem/hvac136>`_ are the first: UMIErrorCorrect aligns
with ``bwa mem`` and only then groups reads "based on target DNA region (i.e., chromosomal position
and UMI sequence)".

What the position buys
~~~~~~~~~~~~~~~~~~~~~~

**Extra key bits, for free.** Two molecules that drew the same UMI are one group under
``(sample, cell, UMI)`` and two groups under ``(position, UMI)``, so mapping first raises the
effective barcode length by however many bits the locus contributes. That is not a small
correction when the barcode is short: TSO500's 5 nt UMI is 1,024 barcodes, which does not identify
a molecule in a panel of tens of thousands, and the protocol is position-aware downstream for
exactly this reason (:doc:`layouts`).

It is worth the most when :doc:`barcode space <barcode_space>` is saturated, and migec tells you
whether it is -- ``checkout`` reports occupancy, and warns when the space is full enough that the
collision arithmetic stops being reliable.

What it costs
~~~~~~~~~~~~~

* **You align every raw read.** At 4 reads per molecule that is 4x the alignment; at 30 it is 30x.
  Collapsing first aligns one record per molecule, and it is the *clean* one.
* **A mismapping becomes a grouping error.** Position is only a key if it is right. A read placed
  in the wrong repeat copy joins the wrong molecule, and nothing downstream can see that it did.
* **It needs a reference.** Repertoire sequencing largely does not have a useful one -- every read
  of a TRB library maps to the same handful of V genes, so the locus adds close to zero bits, which
  is the case migec was built for.
* **The aligner sees uncorrected sequence**, so its own error model is doing work the consensus
  would have done better.

What migec does instead
~~~~~~~~~~~~~~~~~~~~~~~

The discriminating power of the position is really the discriminating power of *the sequence at*
that position -- and the payload is already in hand before any aligner runs. So ``assemble``
sub-clusters within a barcode group on the reads' own disagreement, splitting a group whose reads
carry two co-segregating haplotypes at the threshold :doc:`the permutation null puts it at <nulls>`
(8.68, not the nominal 2.00). Two molecules that collided on a UMI *and* came from different loci
disagree at many positions at once, which is precisely what that test detects, and it detects it
without a reference, an index or a second alignment.

Where that argument runs out is a barcode so saturated that two fragments of two molecules share no
sequence at all -- the case ``--contig`` is careful about, and the reason
``expected_molecules_per_group`` is reported (:doc:`fragmented`).

UMI-aware means two different things
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 26 54

   * - tool
     - what "UMI support" means there
     - feed it a consensus?
   * - ``minimap2``, ``bwa``, ``minibwa``
     - carries ``RX``/``CB``/``MI`` through to the SAM
     - **yes** -- this is the intended input, and it is 1/N of the alignment work
   * - ``salmon``, ``kallisto``
     - nothing; they count records
     - **yes**, plainly -- one record is one molecule, so ``NumReads`` already is a molecule count
   * - ``fgbio``, ``UMI-tools``, ``UMIErrorCorrect``, ``gencore``
     - groups an aligned BAM by ``(position, UMI)`` and consenses
     - **no** -- it is an alternative to ``assemble``, not a stage after it
   * - ``alevin``, ``bustools``, ``STARsolo``
     - reads the cell barcode and UMI out of a **raw** barcode read and deduplicates
     - **no** -- the barcode read no longer exists, and both collapsing merges molecules twice

The rule underneath the table: a tool that only *transports* the barcode composes with migec, and a
tool that *deduplicates* on it replaces a stage of migec. Running two deduplicators in series
counts each molecule once and then collapses the result again, which silently merges distinct
molecules that happen to share a sequence.

Counting after collapsing
-------------------------

**One consensus is one molecule**, so a plain quantifier's read count already is a molecule count
and the UMI-aware modes must not be run on top of it:

.. code-block:: bash

    salmon quant -i tx.idx -l A -r cons/S1.consensus.fq.gz -o quant/
    kallisto quant -i tx.idx -o quant/ --single -l 200 -s 20 cons/S1.consensus.fq.gz

Never: do not feed a consensus FASTQ to ``alevin``, ``bustools`` or ``STARsolo``. Those read the
cell barcode and UMI out of a *raw* barcode read and deduplicate themselves; migec has already
deduplicated, and the barcode read no longer exists. Running both counts each molecule once and
then collapses the result again, which silently merges distinct molecules that happen to share a
sequence. Pick one: either migec collapses and you quantify plainly, or they collapse and migec is
not in the pipeline.

For alignment, the tags are what a variant caller or a duplicate-aware tool needs:

.. code-block:: bash

    minimap2 -ax sr -y ref.fa cons/S1.consensus.fq.gz  | samtools sort -o S1.bam
    bwa mem -C ref.fa cons/S1.consensus.fq.gz          | samtools sort -o S1.bam
    minibwa map -y -t8 ref.fa cons/S1.consensus.fq.gz  | samtools sort -o S1.bam
    samtools view S1.bam | grep -o 'MI:Z:[^\t]*' | sort -u | wc -l   # molecules, from the BAM

Reproducing this table
----------------------

``tests/unit/test_downstream.py`` asserts the parts that need no external tool -- that the name is
self-sufficient, that the tags are tab-separated, and that the comment is a valid SAM tag list. The
aligner and quantifier rows are skipped unless the tool is on ``PATH``, so the same file is the
check on a machine that has them.
