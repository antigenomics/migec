plot -- twenty QC figures
=========================

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
     - ``<sample>.cell_rank.tsv``
     - **Cell Ranger's barcode rank plot**: barcodes sorted by the number of distinct UMIs they
       carry, log-log, with the called cells and the ambient background drawn in different
       colours. Never reads -- see below.
   * - ``molecule_rank``
     - ``<sample>.rank.tsv``
     - The same curve one level down: reads per molecule, largest first.
   * - ``mig_size_spectrum``
     - ``<sample>.sizes.tsv``
     - Molecules and the reads they account for, against :math:`\log(1 + \text{size})`. The two
       series peak in different places on an over-sequenced library, and that gap is the finding.
   * - ``mig_size_zipf``
     - ``<sample>.sizes.tsv``
     - Molecule size against rank, log-log. A straight line is Zipf; amplification bias bends it.
   * - ``umi_error_children``
     - ``<sample>.umi_errors.tsv``
     - How many distinct error children a molecule spawned, and how many reads were in them,
       against the depth of the parent. The dashed line is :math:`3L` -- every distinct child
       there can be -- and where the points leave it is where the barcode neighbourhood filled.
   * - ``umi_error_rate``
     - ``<sample>.umi_errors.tsv``
     - The barcode error rate the two series above imply, per depth, with what ``refine`` reports
       drawn across them. One decade is ten Phred, so this is directly comparable with the
       barcode's own reported Q.
   * - ``sample_umis``
     - ``checkout.summary.tsv``
     - Unique UMIs and reads per **sample** barcode -- the multiplexed analogue of the barcode
       rank plot.
   * - ``refine_coverage``
     - ``refine.coverage.tsv``
     - Molecules per MIG size after barcode correction.
   * - ``consensus_quality``
     - ``assemble.quality_by_depth.tsv``
     - Emitted quality against depth, as a box per depth bin. It flattens at the RT floor, not at
       the instrument's Q.
   * - ``consensus_error``
     - ``<sample>.mig.tsv``
     - The posterior error the consensus itself achieved, before the floor is added.
   * - ``consensus_layout``
     - ``<sample>.mig.tsv``
     - Consensus lengths, and how many contigs each barcode produced under ``--contig``.
   * - ``assemble_coverage``
     - ``assemble.coverage.tsv``
     - Groups per MIG size, as assembled.

The four familiar ones
----------------------

**The barcode rank plot** is deliberately on `Cell Ranger's
<https://www.10xgenomics.com/support/software/cell-ranger/latest/advanced/cr-ab-barcode-rank-plot>`_
axes, because it is the figure every user of a droplet protocol already knows how to read: barcodes
on the x axis sorted by content, that content on the y, both logarithmic, and a knee where real
cells stop and ambient RNA starts.

.. warning::

   The y axis is **unique UMIs**, never reads. One over-amplified molecule would otherwise put an
   empty droplet high on the curve -- which is the exact artefact the plot exists to make visible,
   so drawing reads there hides the thing you are looking for. ``refine`` calls the cells with
   OrdMag and reports the knee beside it; the panel draws the call on the curve rather than
   describing it in a caption.

The multiplexed analogue is ``sample_umis``: unique UMIs and reads per sample barcode, off
``checkout.summary.tsv``. Same question, one compartment up.

**The MIG size spectrum** draws both series -- how many molecules were seen *n* times, and how many
reads those molecules account for -- because they peak in different places the moment a library is
over-sequenced. Most molecules are shallow; most reads are in the deep ones. A figure with only the
first says the library is fine and a figure with only the second says it is saturated. The x axis is
:math:`\log(1 + \text{size})` rather than :math:`\log(\text{size})` so that a molecule seen once
still has a place on it.

**The rank/Zipf curve** is the cumulative count of that spectrum, read from the deep end down. This
is why ``<sample>.sizes.tsv`` is written at **exact** sizes and not in power-of-two bins: four bins
make four steps, and a straight line cannot be told from a bent one. The table costs one row per
distinct depth -- a few thousand on any real library -- not one row per molecule.

**Consensus quality against depth** is a box, and this is not a stylistic choice. Emitted quality is
discrete and capped at the RT floor, so at any real depth every molecule sits on one or two
integers. A scatter of that draws a flat line whether the bin holds ten molecules or ten million,
and thinning the scatter to keep the SVG small throws away the tails, which were the only thing the
cloud could have shown that the line does not. The quartiles come off an exact ``(depth, quality)``
count grid that ``assemble`` accumulates per bucket -- both axes are small integers, so the whole
joint distribution is 61 counters per power-of-two depth bin and there is nothing to sample.

Publication defaults
--------------------

Every figure is drawn to go straight into a paper or a dark-mode README without editing:

* **Transparent background.** The SVG carries no background rectangle, so it sits on a white page,
  a dark page, or a slide.
* **One ink colour, mid grey.** ``#808080`` reads on both, which is what lets a single file serve
  a light README, a dark README and print instead of three renders.
* **The legend is inside the plot box.** A key in the margin makes every figure wider than its
  data; ``migec plot`` puts it in the corner with the most space and the frame stays 760x520.
* Data colours stay ColorBrewer Dark2, which is qualitative and colour-blind safe.

.. code-block:: bash

   migec plot asm/ -o figs/ --format pdf     # pdfcairo, 6.4 x 4.4 in
   migec plot asm/ -o figs/ --format png     # transparent pngcairo

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
