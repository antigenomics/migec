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

    minimap2 -ax sr -y ref.fa cons/S1.consensus.fq.gz | samtools sort -o S1.bam
    bwa mem -C ref.fa cons/S1.consensus.fq.gz          | samtools sort -o S1.bam
    samtools view S1.bam | grep -o 'MI:Z:[^\t]*' | sort -u | wc -l   # molecules, from the BAM

Reproducing this table
----------------------

``tests/unit/test_downstream.py`` asserts the parts that need no external tool -- that the name is
self-sufficient, that the tags are tab-separated, and that the comment is a valid SAM tag list. The
aligner and quantifier rows are skipped unless the tool is on ``PATH``, so the same file is the
check on a machine that has them.
