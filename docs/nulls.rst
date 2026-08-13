Three permutation nulls
=======================

Three numbers in this pipeline come out of a derivation, and each derivation assumes something the
data has never been asked about. Every one of them has a permutation that measures the same
quantity while assuming nothing at all, and the gap between the two is the size of the assumption.

.. code-block:: bash

   python scripts/permutation_nulls.py --reads SRR1763769_2.fastq.gz --out x3/ \
       --cycles 32 --window 180

``SRR1763769`` again: 2.12 M reads of an HIV-1 Primer ID library, 94.8% assigned on the pattern
:doc:`suggest <barcode_space>` recovers unaided, **125,369 distinct 9 nt barcodes** — 47.8% of the
262,144 nominal space, which is the regime where every one of these questions bites.

.. contents::
   :local:
   :depth: 1

Are the barcode positions independent?
--------------------------------------

``P_coll = Π_j Σ_a p_j(a)²`` — the number that sets the effective barcode length, the saturation
warning and the collision-corrected molecule count — assumes the positions are independent. A
synthesiser couples each step separately, so they need not be, and the assumption has never been
checked. Nothing external is needed to check it: compare the collision probability of the
*observed joint* distribution over ``k`` adjacent positions against the product of that same data's
own marginals.

.. list-table::
   :header-rows: 1

   * - k
     - windows
     - observed / independent
   * - 1
     - 9
     - 1.000×
   * - 2
     - 8
     - 1.004×
   * - 3
     - 7
     - 1.009×
   * - 4
     - 6
     - 1.016×
   * - 5
     - 5
     - 1.024×

The excess grows linearly in ``k``, at **1.0051× per added position**, which extrapolates to
**1.04× over the full 9 nt**: an effective length of 8.98 nt against the 9.01 nt the marginal
product claims. Position independence is, on this library, correct to within 4% of the collision
rate — so it is *not* the explanation for the 1.86× collision excess
:doc:`measured from the sequences <barcode_space>`. That excess is the read threshold: a collided
barcode carries two molecules' reads and so is over-represented among the MIGs large enough for a
split to be visible at all.

.. warning::

   This null has little power on a saturated library. It runs on *distinct* barcodes, so a barcode
   drawn twice appears once — and collapsing repeats removes exactly the frequency variation being
   measured. At 47.8% occupancy the observed set is close to the whole space and is nearly uniform
   by construction. The 1.04× is a **lower bound**, and settling this properly needs a sparse
   library, not a deeper one.

How many distance-1 pairs are really parent and child?
------------------------------------------------------

The barcode error rate is read off the *excess* of barcode pairs at Hamming distance 1 over a
chance background, and both halves of that subtraction are derived. Two permutations, because the
two halves fail differently.

**The background: shuffle the columns.** Shuffling each position independently across the observed
barcodes keeps every marginal exactly and destroys both the dependence and the error-child
structure. What is left is chance.

.. code-block:: text

   pairs at distance 1   839,218
   expected by chance    773,684    (20 column shuffles)
   excess                 65,534    -- 8.5%

At 47.8% occupancy, **92% of all distance-1 pairs are coincidence**. This is the arithmetic behind
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
     - 517,387
     - 505,195
     - 527
     - 12,192
     - 23.1
   * - 5
     - 241,386
     - 222,066
     - 435
     - 19,320
     - 44.5
   * - 10
     - 204,454
     - 185,072
     - 406
     - **19,382**
     - 47.7
   * - 20
     - 181,692
     - 163,155
     - 366
     - 18,537
     - 50.6
   * - 50
     - 139,698
     - 126,510
     - 283
     - 13,188
     - 46.5

The excess plateaus at ~19,400 pairs from a ratio of 5 upward, which is what a population of
genuine children looks like: a child is typically *much* smaller than its parent, so raising the
ratio removes chance pairs without removing real ones. That is **~19,400 error children among
125,369 barcodes, counted without a model** — and it is a floor, because a child whose count
happens to resemble its parent's is invisible to this statistic.

Solving the same distance-1 excess for a per-base rate:

============================================ ==========
 barcode error, from the permutation excess   3.43e-03
 ...if every distance-1 pair were a child     7.97e-02
 ``checkout``'s analytic estimate             7.97e-04
 Phred + polymerase prediction                2.07e-03
============================================ ==========

The permutation estimate lands within **1.7×** of what the reported Phred and the polymerase
predict. ``checkout``'s analytic estimate is **2.6× below** the prediction — on a library it
already flags unreliable, in the direction it is known to fail. Replacing the derived background
with the column-shuffle background recovers most of that collapse, which is the change M3's error
model should carry.

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

1,661 MIGs of 10–200 reads had two or more callable positions:

.. list-table::
   :header-rows: 1

   * - target false positive
     - threshold
     - MIGs called
     - called
   * - 0.05
     - 4.83
     - 120
     - 7.22%
   * - 0.01
     - 9.61
     - 21
     - 1.26%
   * - 0.001
     - 32.42
     - 1
     - 0.06%

A nominal ``p < 0.01`` — score > 2.00, which is what the derivation gives — calls **455 MIGs,
27.39%**. The permutation puts the 1% false-positive point at **9.61**, calling 1.26%. The derived
threshold over-calls by **22×**, and every one of those splits would have become a spurious extra
molecule.

That is the threshold M1 uses. It is a measured false-positive curve, not a derivation, and it is
the whole reason X3 blocks the error model rather than following it.

The test is **two-sided**, and it has to be. At a 50/50 split which allele is the "major" one is a
coin toss taken separately in each column, so the two columns of a genuine doublet come out
*anti*-correlated as often as not, and a one-sided test scores the strongest evidence there is as
nothing at all. Minor-with-major is the same evidence in the opposite phase; both are tested and
the count is halved for it, which is where 9.91 became 9.61.

One consequence worth stating, because it is a floor rather than a tuning knob: the strongest
evidence a pair of columns can carry is ``log10 C(n, n/2)``, so **a 50/50 split needs about 34
reads before it can clear 9.61 at all**. Below that the data cannot distinguish a subclone from a
pair of bad reads at a 1% false-positive rate, and migec does not pretend otherwise.

.. note::

   The curveball null is what makes the difference, not the hypergeometric. Permuting each column
   independently keeps the per-position error count but hands every read an average error load —
   and real reads do not have an average error load. ``tests/synthetic/test_nulls.py`` has the
   minimal case: two reads carrying a minor base at every position score high under the
   hypergeometric and sit squarely in the middle of the both-margins null.

What this settles
-----------------

- **Position independence is fine** at 4% of the collision rate, so ``Π_j m_j`` stays. The 1.86×
  collision excess is the read threshold, not the barcode.
- **The distance-1 background is 92% chance** at half occupancy and must be permuted, not derived.
  Doing so moves the barcode error estimate from 2.6× below the Phred + polymerase prediction to
  1.7× above it.
- **The split threshold is 9.61, not 2.00.** The derivation over-calls by 22× because reads are not
  exchangeable.
