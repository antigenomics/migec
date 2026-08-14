The RT/PCR floor — the quality no consensus can beat
=====================================================

A consensus over *c* reads suppresses sequencing error. It suppresses nothing that was already in
the molecule when the first PCR cycle started: an RT miscall, or a polymerase error in an early
cycle, is in every read of the MIG and every consensus reproduces it faithfully, at high
confidence, because every read agrees.

So the residual error of a consensus flattens out:

.. math::

   e_\text{out}(c) \longrightarrow p_\text{floor} \quad \text{as } c \text{ grows}

and ``p_floor`` is the cap on every quality migec is allowed to emit. It was a guess spanning two
orders of magnitude — 10\ :sup:`-4`, 10\ :sup:`-5` or 10\ :sup:`-6`. Experiment **X2** measures it.

.. code-block:: bash

   uv pip install pysam      # only for X1; X2 needs nothing extra
   python scripts/quality_floor.py --reads SRR1763769_2.fastq.gz --out x2/ --window 180

The script finds the UMI by per-cycle entropy, extracts it with ``migec checkout`` rather than a
bespoke parser, builds a majority consensus per MIG, and scores it against the library's modal
sequence.

What the estimator is not
-------------------------

It is **not** a least-squares fit of ``p_floor + a/c``. That model is wrong for a majority-vote
consensus: the sequencing residual is the probability that *most* reads carry the same wrong base,
which falls roughly geometrically in *c*, not as 1/*c*. Regressing on 1/*c* lets the smallest,
noisiest bin set the intercept — on simulated data with a known floor it returned a **negative
probability**. The floor is where the curve flattens, so it is measured there, with a Poisson
interval.

Three things are excluded, each for a reason that showed up in the data:

**Ties are not calls.** Two reads that disagree have no majority. Resolving that by whichever base
came first turns a coin flip into a confident wrong base, and at even MIG sizes it dominated
everything.

**Polymorphic positions are not error.** Only positions where the molecules agree are scored. On a
viral quasispecies a real 20% variant would contribute 0.2 to the "error" rate and swamp a floor of
10\ :sup:`-4` by three orders of magnitude. The threshold must sit well above 1/(molecules), or a
position where a *single* molecule erred looks polymorphic and is excluded — dropping exactly the
positions the floor lives at. The script refuses a threshold that tight.

**A wholly divergent MIG is a different template.** Over 180 bases at a floor of 10\ :sup:`-4` the
expected disagreement is 0.02 bases, so a MIG differing at 20% of positions is not an erroneous
copy — it is another region, an off-target product, or an indel-shifted read, and we model no
indels anywhere. On the control below **0.6% of MIGs sit past 20% divergence and contribute 81% of
every mismatch in the dataset.**

Calibration
-----------

The instrument is checked against a floor injected by ``tests/synthetic/_sim.py``, which puts RT
errors into every read of a molecule:

.. list-table::
   :header-rows: 1
   :widths: 25 25 35 15

   * - injected
     - recovered
     - 95% interval
     - covers truth
   * - 1·10⁻⁴
     - 9.50·10⁻⁵
     - [7.4·10⁻⁵, 1.2·10⁻⁴]
     -
   * - 1·10⁻⁵
     - 9.62·10⁻⁶
     - [3.9·10⁻⁶, 2.0·10⁻⁵]
     -
   * - 0
     - 0
     - [0, 5.1·10⁻⁶]
     -  (bound)

It resolves the decade the project needs to settle. ``tests/synthetic/test_quality_floor.py`` keeps
that true.

The measurement
---------------

``SRR1763769`` — 2.12 M reads of an HIV-1 Primer ID library (PRJNA272736). ``checkout`` recovered
the layout with no help: a 9 nt Primer ID at 2.00 bits per cycle followed by the conserved cDNA
primer ``CAGTTTAACTTTTGGGCCAT``, and assigned **95.0%** of reads. 125,236 distinct UMIs, 52,703
MIGs of 2–200 reads, 179 of 180 positions monomorphic.

.. list-table::
   :header-rows: 1
   :widths: 20 14 18 24 12

   * - MIGs of ≥
     - mismatches
     - bases
     - p_floor [95% CI]
     - Q cap
   * - 3 reads
     - 2,864
     - 5,172,364
     - 5.54·10⁻⁴ [5.3, 5.7]
     - 32.4
   * - 5 reads
     - 948
     - 3,323,502
     - 2.85·10⁻⁴ [2.7, 3.0]
     - 35.2
   * - 10 reads
     - 578
     - 2,668,104
     - 2.17·10⁻⁴ [2.0, 2.4]
     - 36.3
   * - 20 reads
     - 437
     - 2,359,001
     - 1.85·10⁻⁴ [1.7, 2.0]
     - 36.9
   * - 50 reads
     - 329
     - 2,026,812
     - 1.62·10⁻⁴ [1.5, 1.8]
     - 37.4
   * - **80 reads**
     - **260**
     - **1,685,998**
     - **1.54·10⁻⁴ [1.4, 1.7]**
     - **37.6**

The answer, and what it settles
-------------------------------

**The floor is of order 10⁻⁴, not 10⁻⁶.** At MIGs of 80 reads or more it is
1.54·10⁻⁴ with a 95% upper bound of 1.74·10⁻⁴. The guess that it might be 10\ :sup:`-6` is
excluded by more than two orders of magnitude for a protocol with an RT step, and with it any
emitted quality above **Q40**.
For comparison, the paper this library comes from measures a residual error rate of about 1 in
10,000 (Zhou, Jones, Mieczkowski & Swanstrom, *J Virol* 89:8540–8555, 2015,
`doi:10.1128/JVI.00522-15 <https://doi.org/10.1128/JVI.00522-15>`_).

The cut is not doing the measuring — the estimate is identical at 5%, 10% and 20% divergence
thresholds, and 9% lower at 2%.

.. warning::

   **The curve is still declining at 80 reads, so 1.54·10⁻⁴ is an upper bound, not a plateau.**
   One reason is visible in our own output: the 9 nt Primer ID gives an effective space of 250,902,
   and this library occupies **49.6%** of it — about 30,900 collided pairs, with the Poisson
   correction putting the true molecule count at 171,890 against 124,436 observed.
   ``checkout`` flags the sample ``saturated``. A MIG that is really two templates has a
   consensus that is a mixture, and those mismatches are counted here as error. The true floor for
   this chemistry is therefore somewhere at or below 1.5·10⁻⁴, and a library with a longer barcode
   would measure it more sharply.

Where the other two classes come from
-------------------------------------

``--rt-error`` names a chemistry rather than asking for a number, and only the ``rt`` bracket is
measured above. The other two are published polymerase fidelities: Taq 4.3·10⁻⁵, Pfu 2.8·10⁻⁶,
Phusion 2.6·10⁻⁶ and Pwo 2.4·10⁻⁶ per base pair per template duplication (`McInerney, Adams & Hadi
2014 <https://doi.org/10.1155/2014/287430>`_), which is the decade between ``medium`` and ``high``.
The **first** cycle is the one that sets the floor, because only an error made there is copied into
every read of the molecule, and a linear-amplification error runs about 5 ± 1 times the per-cycle
PCR rate (`Shagin *et al.* 2017 <https://doi.org/10.1038/s41598-017-02727-8>`_) — which is why a
class is a bracket and not the enzyme's datasheet figure. A protocol that publishes its own rate is
passed as a rate: TSO500 v2 is 7.37·10⁻⁵. Every one of these is in ``SOURCES.md`` with its
provenance.

.. _pre-amp-auto:

Fitting it per dataset: ``--pre-amp-error auto``
------------------------------------------------

The floor is a property of the enzyme, the cycle count and the chemistry — not a constant — so
``assemble`` will fit it from the dataset in front of it:

.. code-block:: bash

   migec assemble out/S1.fq.gz -o cons/ --pre-amp-error auto

It is the measurement above, run on migec's own consensuses instead of on a bespoke one: keep the
molecules with at least 20 reads, take the library's modal sequence with one vote per molecule,
drop the positions that carry real variation and the consensuses that are a different template
altogether, and the residual is the floor. The fit and every input it rested on are written to
``assemble.pre_amp_error.tsv``, so a floor that looks wrong can be argued with.

Injecting a known floor with ``tests/synthetic/_sim.py`` and reading it back
(``tests/synthetic/test_pre_amp_auto.py``, which is what keeps this table true):

.. list-table::
   :header-rows: 1
   :widths: 14 14 24 14 14 12

   * - injected
     - fitted
     - 95% interval
     - molecules
     - bases
     - mismatches
   * - 1·10⁻⁴
     - 1.20·10⁻⁴
     - [7.94·10⁻⁵, 1.75·10⁻⁴]
     - 1868
     - 224160
     - 27
   * - 1·10⁻⁵
     - 1.58·10⁻⁵
     - [8.14·10⁻⁶, 2.76·10⁻⁵]
     - 6340
     - 760800
     - 12
   * - 0
     - 0
     - [0, 9.56·10⁻⁶]
     - 3197
     - 383640
     - 0

The point estimate sits a little above the injected RT rate because the simulator also puts
early-PCR errors into the molecule — which are pre-amplification error too. What the floor is
is *what survives a consensus*, whatever made it.

.. warning::

   **``auto`` refuses more often than it answers, and that is the feature.** It needs a clonal
   template: on a diverse library every position is polymorphic, "disagrees with the modal base"
   means "is a different molecule", and the number that comes out is the library's diversity
   rather than its chemistry. It also needs molecules deep enough that the consensus curve has
   flattened. Both refusals name the class they fell back to (``rt``, 10⁻⁴) and say why, in the
   report and in the TSV:

   * 400 clones at 40 reads per molecule — refused, 120 of 120 positions polymorphic.
   * 1.2 reads per UMI — refused, no molecule reaches 20 reads.

   With no observed mismatch at all the answer is the interval's **upper** end, not zero: a floor
   of zero is a quality of infinity.

``auto`` costs a second assembly. The consensus sequences do not depend on the floor — only the
emitted quality does — so the probe assembly runs, is read back and is deleted, and the real run
emits with the fitted cap. Name the chemistry instead when you know it.

Consequences for M1
-------------------

* The **default** is 10\ :sup:`-4`, not 10\ :sup:`-6`. A wrong default here does not degrade
  gracefully: it inflates every quality above Q40 in the output, and downstream variant callers
  believe them.
* The cap is taken from the *upper* bound of the interval. Claiming a quality the data cannot
  support is the failure that matters; being a decibel conservative is not.
* Any dataset used to fit the floor must be checked for saturation first. At 50% occupancy the
  measurement is contaminated by collisions, and the contamination is in the direction that makes
  the floor look worse than it is.
