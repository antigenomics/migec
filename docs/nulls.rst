Three permutation nulls
=======================

Three numbers in this pipeline come out of a derivation, and each derivation assumes something the
data has never been asked about. Every one of them has a permutation that measures the same
quantity while assuming nothing at all, and the gap between the two is the size of the assumption.

.. code-block:: bash

   python scripts/permutation_nulls.py --reads SRR1763769_2.fastq.gz --out x3/ \
       --cycles 32 --window 180

``SRR1763769``: 2.12 M reads of an HIV-1 Primer ID library, 94.8% assigned on the pattern
:doc:`suggest <barcode_space>` recovers unaided, **124,562 distinct 9 nt barcodes** — 47.5% of the
262,144 nominal space, which is the regime where every one of these questions bites.

.. contents::
   :local:
   :depth: 1

Are the barcode positions independent?
--------------------------------------

``P_coll = Π_j Σ_a p_j(a)²`` — the number that sets the effective barcode length, the saturation
warning and the collision-corrected molecule count — assumes the positions are independent. A
synthesiser couples each step separately, so they need not be.

The null here is not a statistic, it is a **distribution**: the product measure
``q(u) = Π_j p_j(u_j)``, built from the data's own marginals. So the test compares the whole of
``P`` against the whole of ``Q`` by Jensen-Shannon divergence, rather than comparing one functional
of each. Over the full barcode ``KL(P ‖ Q)`` is exactly the multi-information; JSD is its bounded
symmetric form and survives a cell the null gives zero weight.

JSD needs a floor, because on 4\ :sup:`k` cells and *n* barcodes the empirical distribution is
sparse and JSD is bounded away from zero even under perfect independence. A **column shuffle**
draws from ``q`` at exactly the observed *n*, so the floor is measured rather than derived.

Three things had to be right before any of this meant anything
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The first version of this null reported a 1.04× excess. All of it was artefact.

- **An ``N`` is not a fifth base.** Counting it as one let ``m_j`` fall to 0.2466, below the
  mathematical floor of 1/4, and printed an *effective length of 9.01 nt for a 9 nt barcode*. The
  bound is the check. ``N`` is now folded to ``A`` exactly as ``pack_barcode`` stores it — the
  packed key is what every stage groups on, so it is the alphabet the collision rate is about.
- **The plug-in ``Σ p̂²`` is biased upward** by ``(1 − Σp²)/n``, and that bias *grows* as the
  distribution spreads — so it grows with the k-mer width and reads as dependence accumulating
  with k. The U-statistic ``Σ nₐ(nₐ−1)/(N(N−1))`` has none.
- **Distinct barcodes, never barcodes weighted by reads.** Read counts are set by amplification, so
  weighting would let one over-amplified molecule write the composition.

What it says
^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - barcode set
     - n
     - collision ratio over 9 nt
     - JSD z, k = 2 / 3 / 4 / 5
   * - all distinct
     - 124,562
     - 1.0011×
     - −0.6 / 0.9 / 1.3 / −0.7
   * - reads ≥ 2
     - 54,680
     - 1.0103×
     - **7.9 / 11.9 / 11.0 / 10.5**

On all distinct barcodes there is no signal at all — but that is saturation, not independence. At
47.5% occupancy the observed set is close to a complete enumeration of the space, and a complete
enumeration is uniform *by construction* whatever the draw distribution was. Dropping singletons
takes occupancy to 20.9% and removes most of the barcodes that are error children rather than
molecules, and then the dependence is unambiguous.

.. note::

   Dropping singletons is affordable here and is not affordable everywhere — see
   :ref:`shallow-libraries` below. It is used for this measurement, not as a pipeline default.

And the dependence is entirely nearest-neighbour
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Pairwise JSD z-score between positions *j* and *k*, on the reads ≥ 2 set:

.. code-block:: text

           0      1      2      3      4      5      6      7      8
   0       .    7.4   -1.8   -1.3   -8.3   -3.6   -3.1   -1.6   -2.2
   1       .      .   10.2   -2.8   -4.0   -4.5   -4.0   -2.4   -6.0
   2       .      .      .   20.2   -0.9    0.6   -1.1   -0.8    1.5
   3       .      .      .      .    4.9   -3.2   -2.5   -0.8   -0.6
   4       .      .      .      .      .   11.9   -4.7   -2.4   -3.3
   5       .      .      .      .      .      .    0.9    0.6   -2.8
   6       .      .      .      .      .      .      .   17.6    1.7
   7       .      .      .      .      .      .      .      .   33.0

**Every adjacent pair is positive; every non-adjacent pair is zero.** The structure is local and
position-specific, not diffuse — and it is strongest at (7, 8), the two positions abutting the
constant anchor.

That is a synthesis failure. A coupling step that does not fire leaves the barcode one base
**short**, which frameshifts every position after it — and a frameshift is precisely a
nearest-neighbour correlation, largest next to the anchor where the shift is most visible. It is
measurable directly, from where the anchor actually starts:

.. code-block:: text

   anchor at        reads    share   interpretation
           8        5,473    0.55%   1 nt shorter barcode -- one failed coupling
           9      919,070   91.91%   as designed
          10        1,519    0.15%   1 nt longer

   0.76% of anchored reads carry a barcode shorter than designed

migec models no indels anywhere, so a substitutions-only matcher either drops these reads or
captures a barcode shifted by one. Either way the rate is worth knowing, and it is the leading
candidate for the residual dependence.

What it changes: **nothing**. The magnitude is ~1% of the collision rate — an effective length of
8.916 nt against the 8.909 nt the null implies. ``Π_j m_j`` stays, and the 1.86× collision excess
:doc:`measured from the sequences <barcode_space>` is not this. That excess is the read threshold:
a collided barcode carries two molecules' reads and is over-represented among the MIGs large enough
for a split to be visible at all.

How many distance-1 pairs are really parent and child?
------------------------------------------------------

The barcode error rate is read off the *excess* of barcode pairs at Hamming distance 1 over a
chance background, and both halves of that subtraction are derived. Two permutations, because the
two halves fail differently.

**The background: shuffle the columns.** Shuffling each position independently across the observed
barcodes keeps every marginal exactly and destroys both the dependence and the error-child
structure. What is left is chance.

.. code-block:: text

   pairs at distance 1   844,243
   expected by chance    817,358    (20 column shuffles)
   excess                 26,885    -- 3.2%

At 47.5% occupancy, **97% of all distance-1 pairs are coincidence**. This is the arithmetic behind
``err_unreliable``, and here it is measured rather than derived.

**The excess: shuffle the counts.** Permuting the read counts over the *same* distance-1 graph
keeps the barcodes, the graph and the count distribution, and destroys only which count sits on
which node. An error child is far smaller than its parent; two unrelated neighbours are not.

.. list-table::
   :header-rows: 1

   * - count ratio ≥
     - observed
     - null mean
     - null sd
     - excess
     - z
   * - 2
     - 521,137
     - 511,118
     - 413
     - 10,019
     - 24.3
   * - 5
     - 242,467
     - 224,578
     - 335
     - 17,888
     - 53.5
   * - 10
     - 205,158
     - 186,908
     - 464
     - **18,250**
     - 39.3
   * - 20
     - 182,195
     - 164,592
     - 349
     - 17,603
     - 50.5
   * - 50
     - 139,949
     - 127,537
     - 321
     - 12,412
     - 38.6

The excess plateaus at ~18,000 pairs from a ratio of 5 upward, which is what a population of
genuine children looks like: a child is typically *much* smaller than its parent, so raising the
ratio removes chance pairs without removing real ones. That is **~18,000 error children among
124,562 barcodes, counted without a model** — and it is a floor, because a child whose count
happens to resemble its parent's is invisible to this statistic.

Solving the same distance-1 excess for a per-base rate:

============================================ ==========
 barcode error, from the permutation excess   1.44e-03
 ...if every distance-1 pair were a child     8.10e-02
 ``checkout``'s analytic estimate             7.97e-04
 Phred + polymerase prediction                2.07e-03
============================================ ==========

The permutation estimate is **0.70×** the Phred + polymerase prediction. ``checkout``'s analytic
estimate is 0.39× of it — on a library it already flags unreliable, in the direction it is known to
fail. Replacing the derived background with the column-shuffle background halves the gap, which is
the change M3's error model should carry.

.. warning::

   The count-ratio statistic needs count asymmetry to work, so it has nothing to say about a
   library sequenced at 1–3 reads per UMI. See :ref:`shallow-libraries`.

When is a MIG really two molecules?
-----------------------------------

Splitting a MIG into two consensuses is the decision with the worst failure mode in the pipeline:
call it too readily and one molecule becomes two, inflating every count downstream. The threshold
was to come from a Poisson argument over polymorphic positions. That argument treats reads as
exchangeable, and they are not — **a low-quality read carries a minor base at many positions at
once**, which is indistinguishable from a linked subclone if you only look at the columns.

So the null has to preserve *both* margins of the reads × positions minor-allele matrix: the
per-position error count **and** the per-read error load. That is a curveball randomisation
(Strona *et al.* 2014), and it keeps the bad reads bad while destroying the specific co-segregation
a real subclone has. The statistic is the strongest co-segregation over any pair of callable
positions, as a Bonferroni-corrected hypergeometric ``-log10 p``.

3,312 MIGs of 10–200 reads with two or more callable positions, 82,800 randomisations:

.. list-table::
   :header-rows: 1

   * - target false positive
     - threshold
     - MIGs called
     - called
   * - 0.10
     - 3.64
     - 457
     - 13.80%
   * - 0.05
     - 4.97
     - 249
     - 7.52%
   * - 0.01
     - **8.68**
     - 53
     - 1.60%
   * - 0.001
     - 29.27
     - 4
     - 0.12%

A nominal ``p < 0.01`` — score > 2.00, which is what the derivation gives — calls **1,014 MIGs,
30.62%**. The permutation puts the 1% false-positive point at **8.68**, calling 1.60%. The derived
threshold over-calls by **19×**, and every one of those splits would have become a spurious extra
molecule.

.. note::

   **The threshold is a Monte Carlo estimate and its error is quoted.** Bootstrap 95% CI
   **[8.42, 9.14]** over 82,800 null scores. An earlier run at a tenth as many randomisations gave
   9.61, and another 11.66 — a tail quantile estimated from ~80 points is not a constant. Do not
   quote digits past what the interval supports.

The test is **two-sided**, and it has to be. At a 50/50 split which allele is the "major" one is a
coin toss taken separately in each column, so the two columns of a genuine doublet come out
*anti*-correlated as often as not, and a one-sided test scores the strongest evidence there is as
nothing at all. Minor-with-major is the same evidence in the opposite phase; both are tested and
the count is halved for it.

One consequence worth stating, because it is a floor rather than a tuning knob: the strongest
evidence a pair of columns can carry is ``log10 C(n, n/2)``, so **a 50/50 split needs about 30
reads before it can clear 8.68 at all**. Below that the data cannot distinguish a subclone from a
pair of bad reads at a 1% false-positive rate, and migec does not pretend otherwise.

.. _shallow-libraries:

Shallow libraries, where most UMIs have 1–3 reads
-------------------------------------------------

Everything above was measured on a library sequenced deeply enough that a barcode carries tens of
reads. That is not the common case. When the input molecule count is large relative to the read
budget, the MIG size histogram is Poisson-ish with its mass at 1–3 — bulk TCR repertoire profiling
routinely looks like this, and so does shallow 3' single-cell GEX (X1: only 1.5% of 10x groups hold
more than one read at all).

migec runs there and reports honestly, but three of the things on this page stop applying, and it
is better to say so than to quote a number that was calibrated somewhere else:

- **The split threshold is a no-op.** It needs ~30 reads in a group; nothing in a 1–3 read library
  reaches it. That is correct behaviour — there is no evidence to split on — but it means
  ``groups_split`` being zero says nothing about whether the library has doublets.
- **The count-ratio null has no dynamic range.** A parent with 3 reads and a child with 1 is not a
  10× asymmetry, so the ~18,000 error children counted above cannot be counted this way. The
  distance-1 *excess* still works, because it does not use counts.
- **Singleton filtering is not affordable.** It is the right move for the independence null on a
  deep library, where it costs 56% of barcodes. On a shallow one it costs 79% and takes most of the
  real molecules with it — which is why it is a measurement choice here and never a pipeline
  default. ``assemble`` keeps everything: ``--min-reads`` defaults to 1, because a molecule seen
  once is still a molecule, and the correct response to a barcode error is to correct it, not to
  threshold it away.

What migec reports on such a library is what it can support: the coverage histogram, and the fact
that the UMI is buying **counting** rather than error correction.

.. code-block:: text

   molecules   40,176
     expected  1.00 molecules per group from the birthday problem at 0.2% occupancy

   mean emitted quality  Q32.8  (capped at Q40 by the RT floor of 1.0e-04)
   mean consensus error  9.42e-04 before the floor is added

       MIG size      groups    share
              1      31,888    79.4%
            2-3       8,176    20.4%
            4-7         112     0.3%

   warning: 79.4% of molecules were seen once. A consensus over one read is that read --
     the UMI is buying counting here, not error correction

What this settles
-----------------

- **Position dependence is real but tiny and purely nearest-neighbour**, ~1% of the collision rate,
  and traceable to a 0.55% rate of failed couplings that leave the barcode one base short.
  ``Π_j m_j`` stays. The 1.86× collision excess is the read threshold, not the barcode.
- **The distance-1 background is 97% chance** at half occupancy and must be permuted, not derived.
  Doing so moves the barcode error estimate from 0.39× of the Phred + polymerase prediction to
  0.70×.
- **The split threshold is 8.68 [8.42, 9.14], not 2.00.** The derivation over-calls by 19× because
  reads are not exchangeable.
