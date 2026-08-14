Cell Ranger: cell barcodes, cell calling, and what the cell count costs
=======================================================================

migec and `Cell Ranger <https://www.10xgenomics.com/support/software/cell-ranger/latest>`_ both
read a 10x droplet library front to back: given R1 = 16 nt cell barcode + 10 nt UMI, which reads
carry a real GEM barcode, which barcodes are cells, and how much of the library is in them. Those
three questions are comparable and are compared here, on ``sc5p_v2_hs_PBMC_1k`` VDJ-T.

.. code-block:: bash

   cat sc5p_v2_hs_PBMC_1k_t_S1_L00{1,2}_R1_001.fastq.gz > R1.fq.gz
   cat sc5p_v2_hs_PBMC_1k_t_S1_L00{1,2}_R2_001.fastq.gz > R2.fq.gz
   python scripts/compare_cellranger.py --r1 R1.fq.gz --r2 R2.fq.gz \
       --cellranger-dir cellranger_published/ --whitelist 737K-august-2016.txt \
       --out /tmp/cr --threads 8 --tsv assets/cellranger.tsv

Cell Ranger is run **both ways**. 10x publish their own Cell Ranger 5.0.0 output for this exact
library, and ``cellranger vdj`` 10.1.0 was run here on the same reads
(``scripts/cellranger_vdj.sbatch``, 16 cores) so the comparison has a cost axis and a second
version to bound drift against.

.. list-table:: Cell Ranger against itself, five major versions apart
   :header-rows: 1
   :widths: 22 13 13 13 13 26

   * - version
     - cells
     - contigs
     - shared
     - Jaccard
     - source
   * - 5.0.0
     - 479
     - 943
     - 477
     - 0.9938
     - published by 10x
   * - 10.1.0
     - 478
     - 939
     - 477
     - 0.9938
     - run here, 524.5 s / 936 MB

**Version drift is negligible**, so the migec-against-Cell-Ranger gap below is not a version
artefact. That control is why 10.1.0 was run at all.

.. note::

   **5.0.0's published ``Median TRB UMIs per Cell`` of 12.0 is not reproducible from its own
   deposited tables**, under either denominator -- both give 11.0. 10.1.0 reports 11.0 for TRB and
   4.0 for TRA, and both are exactly the median over **all** cells in its own tables; 5.0.0's TRA
   figure of 5.0 is the median over **cells that have that chain**. So 10.1.0 changed the
   denominator and is self-consistent, while 5.0.0's TRB value matches neither. Derive these from
   the contig tables and say which denominator you used; do not quote ``metrics_summary`` for them.

The result
----------

.. list-table::
   :header-rows: 1
   :widths: 22 13 13 13 13 13 13

   * - tool
     - reads
     - valid barcode
     - cells
     - shared
     - molecules
     - reads in cells
   * - Cell Ranger 10.1.0
     - 6,301,573
     - **90.60%**
     - 478
     - 469
     -
     - **86.40%**
   * - migec, no whitelist
     - 6,301,573
     - 88.27%
     - 890
     - 469
     - 496,373
     - 84.26%
   * - migec, 737K whitelist
     - 6,301,573
     - 88.90%
     - 888
     - 469
     - 467,810
     - 84.88%

What it costs
-------------

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - stage
     - wall clock
     - peak RSS
     - what it produces
   * - ``cellranger vdj`` 10.1.0, 16 cores
     - 524.5 s
     - 936 MB
     - cell calls **and** per-cell contigs, annotated
   * - migec ``checkout`` + ``refine``, 8 threads
     - **35.0 s**
     - **679 MB**
     - cell calls
   * - \+ ``assemble --contig`` + arda
     - ~90 s total
     - 679 MB
     - \+ per-molecule consensus, annotated

Never: **35 s against 524 s is not like for like, and the table says so.** Cell Ranger's single
invocation also assembles and annotates a contig per cell; the stages that answer the three axes
above are 35 s of migec. The comparable end-to-end number is the third row -- **~90 s against
524.5 s, still 5.8x** -- and it produces per-*molecule* resolution that Cell Ranger's output does
not have. Both numbers are on the same 6,301,573 read pairs, but on different hardware: Cell Ranger
on a 16-core cluster node, migec on 8 threads of a laptop, so read the ratio and not the seconds.

**migec calls 1.86x the cells and loses 1.5 points of reads-in-cells.** That is the whole finding,
and it says the extra barcodes are nearly empty: 419 barcodes beyond Cell Ranger's set hold about
one and a half percent of the library between them. A cell count is not an accuracy figure on its
own.

The two cell sets are **not nested and do not measure the same thing**. migec's gate is molecules of
any sequence; Cell Ranger's is "assembled a productive V(D)J contig". A B cell or a monocyte with
plenty of molecules is correctly a migec cell and correctly not a Cell Ranger VDJ cell. So the table
carries five counts -- 478, 888, 469 shared, 419 migec-only, 9 Cell-Ranger-only -- and never a
ratio. The 9 Cell Ranger cells migec misses sit below the OrdMag threshold on raw molecule counts.

Note: **the whitelist is worth 0.63 points of validity and two cells.** Snapping off-list barcodes
onto the list lifts read validity from 88.27% to 88.90% against Cell Ranger's 90.60%, and moves the
cell count by two. The remaining 1.7-point gap is reads whose R1 is not a barcode at all: the
deepest off-list 26-mers are phase shifts of one fixed sequence (``GGTCCGTCTTGCGCCG`` and its
rotations), which ``refine``'s pass 0 remaps rather than drops.

Per-cell receptor chains, without a per-cell assembler
------------------------------------------------------

Cell Ranger assembles a contig per **cell** from every read of a barcode, then annotates it. migec
assembles a consensus per **molecule**, annotates each one with
`arda <https://github.com/antigenomics/arda>`_, and lets the cell's chain be a vote over its
molecules. Same question -- which chain is in which cell -- reached with and without an assembler.

.. code-block:: bash

   migec assemble ref/PBMC.fq.gz -o asm/ --contig --min-reads 30
   python scripts/compare_cellranger_chains.py --consensus asm/PBMC.consensus.fq.gz \
       --cellranger-dir cellranger_published/ --min-reads 30 --out /tmp/chains \
       --tsv assets/cellranger_chains.tsv

.. list-table::
   :header-rows: 1
   :widths: 12 20 16 14 19 19

   * - locus
     - Cell Ranger chains
     - migec chains
     - shared
     - chain recall
     - junction agreement
   * - TRA, against 5.0.0
     - 426
     - 451
     - 426
     - **1.0000**
     - 0.9507
   * - TRB, against 5.0.0
     - 469
     - 474
     - 468
     - **0.9979**
     - **0.9915**
   * - TRA, against 10.1.0
     - 424
     - 450
     - 424
     - **1.0000**
     - 0.9505
   * - TRB, against 10.1.0
     - 469
     - 473
     - 468
     - **0.9979**
     - **0.9915**

**Every TRA chain and all but one TRB chain Cell Ranger found is also found by migec plus arda**,
in 22 seconds over 47,584 consensuses, with no per-cell assembly step anywhere -- and the result is
the same against both Cell Ranger versions, which is the point of scoring it twice. Recall is the
metric that leads because a missed chain is the unrecoverable error; junction agreement is
secondary and is scored only over the chains both tools called, which is a self-selecting
denominator.

migec calls slightly *more* chains than Cell Ranger on both loci (451 against 426, 474 against
469). Those extra calls are not scored here -- Cell Ranger's set is the reference, not the truth,
and a second productive TRA is allelic inclusion rather than an error.

.. warning::

   **Depth does not buy junction coverage on this chemistry, and the obvious arithmetic says it
   does.** The tempting model places each read uniformly over the ~508 nt amplicon: a 90 nt read
   then spans the median 42 nt junction with probability 0.114, and 30 reads give
   :math:`1 - 0.886^{30} = 0.975`. Measured, that is wrong. Reads of one ``(CB, UMI)`` are
   co-terminal in 92% of 10x groups (:doc:`fragmented`), so a molecule is a **pile at one
   position**, not a tiling, and its consensus covers one window however deep it is. At
   ``--min-reads 30`` the mean consensus is **204 nt, not 508**, and **7,855 of 47,584 molecules
   (16.5%)** carry a cell, a locus and a junction -- close to the single-window 0.32 the geometry
   predicts, nowhere near 0.975. The depth cut is still right, because a deep pile gives a clean
   consensus; it just does not extend one.

Per-cell contigs, reference-free
--------------------------------

The window a molecule covers is one window, but *different molecules start at different
positions*, so a cell's molecules tile its transcript. That is measurable before any assembler
exists, and it is the ceiling everything below is scored against: **99.90% of the 25-mers of Cell
Ranger's 943 filtered contigs are already present in their own cell's raw migec molecules**, and
942 of 943 CDR3 nucleotide sequences appear verbatim in one of them. The contig is in the data.

``arda cells`` assembles it -- adapter trim, 25-mer seeds, verified overlaps, union-find layout,
haplotype phasing, weighted column consensus -- with **no germline reference until the finished
contigs are annotated**.

.. code-block:: bash

   migec assemble ref/PBMC.fq.gz -o asm_all/ --contig --min-reads 1
   python scripts/compare_cellranger_contigs.py --consensus asm_all/PBMC.consensus.fq.gz \
       --cellranger-dir cellranger_published/ --out /tmp/contigs \
       --tsv assets/cellranger_contigs.tsv

.. list-table::
   :header-rows: 1
   :widths: 24 13 13 16 10 12 12

   * - variant
     - k-mer coverage
     - contigs at >=90%
     - CDR3 exact
     - N50
     - chain recall
     - doublets
   * - molecules, no assembly
     - 0.9990
     - 943
     - 942 (0.9989)
     - --
     - --
     - --
   * - ``arda cells``, default
     - **0.9759**
     - **892**
     - **933 (0.9894)**
     - 536
     - **0.9777**
     - 17
   * - no phasing
     - 0.9663
     - 868
     - 926 (0.9820)
     - 460
     - 0.9714
     - 12
   * - no adapter trim
     - 0.9525
     - 879
     - 907 (0.9618)
     - 580
     - 0.9491
     - 16

Per chain at the default: **TRA 454/464 (0.9784), TRB 479/479 (1.0000)**. 479 cells and 249,635
molecules in 23 s.

.. warning::

   **``--min-reads`` throws away the tiling.** A one-read molecule is a poor consensus and one
   more window of the transcript, and the second thing is what the assembly needs. This axis runs
   at ``--min-reads 1`` on purpose; the depth cut belongs to the per-molecule route above, not to
   this one.

.. warning::

   **Contig N50 is a description, never a score.** The no-adapter-trim row has the highest N50 in
   the table, 580 against 536, and the lowest value in every other column. A contig built across
   an adapter is a longer contig and a wronger one.

Doublets: two chains of the same locus
--------------------------------------

One TRA and one TRB in a droplet is a paired T cell. Two TRB is a doublet -- and until the
phasing landed it was invisible **by construction**: the two chains share their constant region,
so the overlap layout puts them in one component and the column consensus averages their
junctions into a sequence that is neither. On a synthetic cell built from two TRB receptors that
produced one 918 nt contig with no callable junction; with the phasing on, both true junctions
come back. On this library it moves chain recall 0.9714 to 0.9777 and doublet candidates 12 to 17.

What separates a real second chain from ambient RNA is **productivity first and count second**:
of the extra chains Cell Ranger agrees with, 60/60 are productive; of those it does not, 20/128
are. An extra chain carried by exactly one molecule is contamination 96-97% of the time. The
sweep that fixes the thresholds, the per-cell table and the QC panels are documented in arda's
`single-cell page <https://arda.readthedocs.io/en/latest/singlecell.html>`_.

.. _not-comparable:

Measured, and not comparable
----------------------------

Each of these was computed and then deliberately kept out of the table, because the number would
have been read as a comparison and is not one.

* **Molecules per cell.** Cell Ranger's ``umis`` counts UMIs incorporated into a *filtered contig* --
  7,623 over 479 cells, median 15 -- while migec counts every molecule of any sequence, a median of
  178 on the same barcodes. That is a different population, not a 12x over-count, so no ratio of the
  two appears.
* **migec's "100% assigned" against "Valid Barcodes 90.6%".** Not the same measurement. On
  ``^XXXXXXXXXXXXXXXXNNNNNNNNNN`` there is nothing to score, so the assigned rate is 100% by
  construction. The comparable quantity is the read share on the whitelist, which is what the table
  reports.
* **Cell Ranger's per-contig ``reads`` as a reads-in-cells figure.** Summing it gives 63.61%, which
  is not any published metric -- it excludes reads in cells that went into no contig. The 86.80%
  comes from ``metrics_summary.csv`` and is not recoverable from the per-contig tables.
* **``Median TRB UMIs per Cell``.** Not cut but corrected -- see the note at the top of the page.
  5.0.0's published 12.0 matches neither denominator in its own tables; 10.1.0 reports 11.0 and is
  self-consistent.

Two traps worth naming
----------------------

.. warning::

   **Both lanes, or neither.** ``L001`` is 3,155,166 read pairs and ``L002`` is 3,146,407; only
   their sum, 6,301,573, is what Cell Ranger reports. ``checkout`` takes one R1 and one R2, so the
   lanes are concatenated first -- and **in the same order for both mates**. ``assemble`` matches
   mates by position and refuses only on a length mismatch, so ``L001+L002`` against ``L002+L001``
   has the right length and mis-mates every pair with nothing flagging it. The script asserts the
   read total before scoring anything.

.. warning::

   **Cell Ranger's AIRR export means the opposite of arda's by the same column name.** In
   ``_airr_rearrangement.tsv``, ``consensus_count`` is READS and ``duplicate_count`` is UMIs --
   verified on all 943 rows against the contig CSV, and the inverse of what arda writes. A join on
   the column names silently swaps reads for molecules, which on this library is a factor of 526.
   The script reads counts only from ``_filtered_contig_annotations.csv`` and asserts the identity
   at runtime rather than trusting it.

Read next
---------

* :doc:`refine` -- OrdMag cell calling, the knee, and ``--cell-whitelist``.
* :doc:`umi_statistics` -- the barcode-space arithmetic behind the validity figures.
* :doc:`downstream` -- what composes with migec's output and what replaces a stage of it.
