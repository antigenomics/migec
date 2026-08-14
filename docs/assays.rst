Assays: what a consensus is worth
=================================

A :doc:`layout <layouts>` says *where the barcode is*. It does not say what the consensus over a
barcode is worth, and that turns out to matter more: the same 12 nt inline UMI serves a repertoire
census and an MRD assay, and the right settings are opposite. This page is the second axis --
``migec sheet --assay NAME`` prints it as a paste-ready recipe, and :doc:`detection` is where the
ctDNA and MRD numbers behind it were measured.

Settings per assay
------------------

``migec sheet --assay NAME`` prints the recipe for one profile, ``--assay all`` for every one:

.. code-block:: console

   $ migec sheet --assay ctdna
   ctdna  (ultrasensitive)
       ...
       layout          tso500   ^NNNNN.....
       payload         uniform -- the barcode carries the whole burden

       migec checkout READS.fq.gz --bc-pattern '^NNNNN.....' --sample S1 -o co/
       migec refine co/S1.fq.gz -o rf/
       migec assemble rf/S1.fq.gz -o as/ --min-reads 3 --rt-error 7.37e-5

Eight profiles, and the sensitivity column is the whole point of the table -- a counting assay must
not inherit a variant-calling threshold, and vice versa:

.. list-table::
   :header-rows: 1
   :widths: 14 18 14 16 38

   * - assay
     - sensitivity
     - ``--min-reads``
     - also
     - why
   * - ``airr`` (``repseq``)
     - counting
     - 1
     -
     - a clonotype seen once is still a clonotype; the rearrangement is a second identifier
   * - ``amplicon`` (``targeted``)
     - sensitive
     - 2
     -
     - a few PCR-amplified regions; families are deep, so 3 is nearly free here
   * - ``exome`` (``capture``)
     - sensitive
     - 2
     -
     - capture duplication is a few-fold, so 3 costs more than it buys
   * - ``ctdna`` (``cfdna``)
     - ultrasensitive
     - 3
     - ``--pre-amp-error 7.37e-5``
     - artifact-limited below 1%; at 1 the dark-G artifact is *additive* to true positives
   * - ``mrd``
     - ultrasensitive
     - 3
     -
     - one known clone, tracked as low as the input DNA allows
   * - ``rnaseq``
     - counting
     - 1
     - ``--fast``
     - deduplicating, not error-correcting
   * - ``10x-gex``
     - counting
     - 1
     - ``--fast``
     - 1-3 reads per (cell, UMI) is the normal case; a threshold deletes the library
   * - ``10x-vdj``
     - counting
     - 1
     - ``--contig``
     - reads under one barcode are random-primed fragments and are not co-terminal

Never: ``--min-reads 3`` on a shallow repertoire library discards **79%** of the barcodes, and
nothing downstream can tell the difference between a molecule that was filtered and one that never
existed. The threshold is for calling variants, not for tidying a count.

Never: **``amplicon`` is not an alias for ``airr``.** A targeted panel of a few PCR-amplified
regions is an amplicon assay too, and it wants the opposite settings -- variant calling on a
uniform payload rather than counting a diverse one. They are two profiles, and ``amplicon`` is the
targeted one because that is what the word means outside immunology.

What the pre-amplification floor actually is
--------------------------------------------

The flag is ``--pre-amp-error``; ``--rt-error`` is kept as an alias because it shipped under that
name, and the class names (``rt``, ``medium``, ``high``) are historical brackets rather than
mechanisms. Never: **there is no reverse transcriptase in a DNA assay**, and four of the eight
profiles above are DNA. The floor is real in both cases; what supplies it differs:

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - library
     - what sets the floor
     - substitution signature
   * - RNA (``airr``, ``rnaseq``, ``10x-*``)
     - a reverse transcriptase miscall, then the first PCR cycle
     - unbiased; 1e-4 is 10x's own figure for their V(D)J RT
   * - DNA (``amplicon``, ``exome``, ``ctdna``, ``mrd``)
     - library-prep damage, then the first PCR cycle
     - ``C>A``/``G>T`` from guanine oxidation (8-oxoG) during acoustic shearing [Costello2013]_;
       ``C>T``/``G>A`` from cytosine deamination [Do2014]_

Both are in the molecule *before* amplification, so every read of the group carries them and no
consensus removes them -- the same argument as for RT, with a different chemistry supplying it.

Never: **that damage signature is not the artifact measured here.** Ours is ``-> G``, the 2-colour
dark-G instrument artifact, and ``--min-reads 3`` removes it because it is carried by uncorrected
small families (:doc:`measured on certified cfDNA <detection>`). A ``C>A``/``G>T`` excess instead
points at oxidative damage during library preparation, which ``--min-reads`` will **not** fix:
damage predates the barcode, so every read of the molecule agrees and the consensus reproduces it
faithfully, at high confidence. The fix is antioxidant handling and enzymatic repair before
ligation, and it belongs in the wet lab.

Two things separate them, and both are free. First, the substitution table: tabulate the types of
the calls ``--min-reads`` removes against the ones it keeps. Second, and more decisive, **read
orientation**. Costello *et al* identify the oxidation artifact by exactly this: 8-oxoG pairs with
A, so the damaged base is read as ``T`` on one strand and its partner as ``A`` on the other, and
the resulting calls pile up in reads of one orientation. A real heterozygous or subclonal variant
is orientation-symmetric. ``migec assemble`` records the strand it normalised in the ``.mig`` flags
and Picard's ``CollectSequencingArtifactMetrics`` scores the bias directly from the BAM.

Note: neither of these is implemented as a filter in migec, and neither should be -- damage is a
property of the library, not of the barcode, so the place to detect it is the aligned BAM where the
orientation still exists. This section says which question to ask, not which flag to pass.

.. [Costello2013] Costello M *et al*. Discovery and characterization of artifactual mutations in
   deep coverage targeted capture sequencing data due to oxidative DNA damage during sample
   preparation. *Nucleic Acids Research* 41(6):e67, 2013. `doi:10.1093/nar/gks1443
   <https://doi.org/10.1093/nar/gks1443>`_

.. [Do2014] Do H, Dobrovic A. Sequence artifacts in DNA from formalin-fixed tissues: causes and
   strategies for minimization. *Clinical Chemistry* 61(1):64-71, 2014.
   `doi:10.1373/clinchem.2014.223040 <https://doi.org/10.1373/clinchem.2014.223040>`_
