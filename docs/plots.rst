plot -- sixteen QC figures
==========================

.. code-block:: bash

   migec plot out/                  # every panel whose table is in out/
   migec plot asm/ -o assets/ --format pdf

``migec plot`` reads no reads. Every panel is a gnuplot script over a TSV a stage already wrote,
which is the same rule the pipeline follows everywhere: a figure must be redrawable from the table
next to it, months later, by someone who no longer has the FASTQ. Nothing is computed here, so a
figure and the number in the report can never disagree.

gnuplot is not a Python package and is not a dependency. ``migec plot`` writes the ``.gp`` scripts
either way and renders them when a gnuplot is on the path; without one it says so and the scripts
still stand alone. Colours are ColorBrewer Dark2, and A/C/G/T keep the same four colours in every
panel that draws bases.

What gets drawn
---------------

.. list-table::
   :header-rows: 1
   :widths: 22 30 48

   * - panel
     - table
     - what it answers
   * - ``umi_pwm``
     - ``checkout.umi_composition.tsv``
     - Did the synthesiser mix the barcode evenly? Four lines at 1/4 is what a UMI looks like; a
       position that drifts is space you are not getting.
   * - ``umi_information``
     - ``checkout.umi_composition.tsv``
     - The same thing in bits, per position. This is what ``eff len`` is the sum of.
   * - ``umi_quality``
     - ``checkout.umi_quality.tsv``
     - The reported Phred over the barcode bases -- the input to the predicted barcode error rate.
   * - ``quality_calibration``
     - ``checkout.quality_calibration.tsv``
     - Observed error against nominal, measured on the pattern's own constant bases. On a 2-colour
       instrument the two do not agree, and this is where you see by how much.
   * - ``coverage``
     - ``checkout.coverage.tsv``
     - The MIG size distribution. If the mass is at 1-3 reads, consensus is buying counting rather
       than error correction, and this is the plot that says so.
   * - ``trimming``
     - ``checkout.trimming.tsv``
     - Payload length after trimming. One spike is a clean trim; a spike one base off its expected
       place is a pattern matched one base off, which nothing else reports.
   * - ``barcode_space``
     - ``checkout.barcode_space.tsv``
     - Nominal :math:`4^L` against the usable space the composition leaves, and the barcodes
       actually observed against both.
   * - ``cycles``
     - ``suggest.cycles.tsv``
     - The per-cycle trace ``suggest`` segments: UMI cycles near 1/4, constant cycles near 1.
   * - ``kmers``
     - ``suggest.kmers.tsv``
     - Overrepresented k-mers -- synthetic sequence still in the reads.
   * - ``cell_rank``
     - ``<sample>.rank.tsv``
     - The barcode rank curve, with the knee where cells stop and ambient begins.
   * - ``refine_coverage``
     - ``refine.coverage.tsv``
     - Molecules per MIG size after barcode correction.
   * - ``consensus_quality``
     - ``<sample>.mig.tsv``
     - Emitted quality against depth. It flattens at the RT floor, not at the instrument's Q.
   * - ``consensus_error``
     - ``<sample>.mig.tsv``
     - The posterior error the consensus itself achieved, before the floor is added.
   * - ``consensus_layout``
     - ``<sample>.mig.tsv``
     - Consensus lengths, and how many contigs each barcode produced under ``--contig``.
   * - ``assemble_coverage``
     - ``assemble.coverage.tsv``
     - Groups per MIG size, as assembled.

Checking a trim, and a consensus
--------------------------------

``suggest`` profiles any FASTQ, so it is also the "what is still in my reads" tool. Run it on the
output of a stage rather than on its input:

.. code-block:: bash

   migec suggest out/S1.fq.gz -o qc-trimmed/     # did the trim remove the primer?
   migec suggest cons/S1.consensus.fq.gz -o qc-cons/
   migec plot qc-trimmed/

The k-mer panel is the answer. An 8-mer occurs by chance about every 65 kb, so a primer that
survived trimming shows up as a run of k-mers each shifted one base from the last, hundreds of
times more often than the reads' own base composition predicts. The report stitches that run back
into the sequence it came from:

.. code-block:: text

   kmer           count   obs/exp    reads  mean pos
   GGGCCATC      20,023     701.2  100.0%      22.1
   TGGGCCAT      20,018     656.4  100.0%      21.0
   TTGGGCCA      20,018     656.4  100.0%      20.0

   overlapping into: CAGTTTAACTTTTGGGCCATCCA

Overrepresentation is measured against the reads' **own** mononucleotide composition, never
against a flat 1/4: a 70% AT library makes every AT-rich k-mer look enriched against uniform, and
the table would then be a description of the GC content rather than a finding.

.. note::

   The k-mer scan covers the whole read, not the profiled prefix, because an adapter that survived
   trimming is at the 3' end -- exactly where the cycle profile does not reach.
