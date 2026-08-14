Barcode space, collisions and the error budget
==============================================

Every ``checkout`` run reports two things it did not measure directly but computed: how crowded the
barcode space is, and what the barcode error rate ought to be. Both are printed because both change
what the output means, and neither is visible from reads-per-UMI.

How big is the space, really
----------------------------

A 12 nt UMI has :math:`4^{12}` = 16,777,216 sequences. Fixed letters written *between* the ``N``
runs are not part of that — ``NNNNtNNNNtNNNN`` captures 12 bases, not 14, and the ``t``\ s are
scored pattern positions like any other constant.

That is the **nominal** space, and it assumes a synthesiser that delivers exactly 25% of each base.
None does. The usable space is set by the *collision* probability per position:

.. math::

   m_j = \sum_a p_j(a)^2, \qquad S_\text{eff} = \frac{1}{\prod_j m_j},
   \qquad L_\text{eff} = -\sum_j \log_4 m_j

Never: This is Rényi-2 (collision) entropy, **not** Shannon. A sequence logo draws Shannon; the
probability that two molecules coincide is :math:`\sum_a p_a^2`. Since :math:`H_2 \le H_1`, using
Shannon overstates the usable space and *understates* collisions — the direction that silently
merges distinct molecules. Both are reported; only the collision form feeds any decision.

``eff len`` in the summary is :math:`L_\text{eff}`, and it is the number to look at rather than the
nominal length. A 12 nt UMI with eight positions fixed and four mixed is worth :math:`L_\text{eff}`
= 4 and a usable space of 256: it will collide constantly, and nothing about it looks wrong until
the molecules have already been merged.

On the HIV Primer ID library in :doc:`quality_floor`, a 9 nt barcode with C at 20.5% instead of 25%
gives :math:`L_\text{eff}` = 8.97 and 250,902 usable sequences against 262,144 nominal — a 4% loss
to the oligo mix. ``bias_loss`` is that number, and ``checkout`` warns past 25%.

The birthday problem, in its useful form
----------------------------------------

Molecules land in the space independently, so occupancy per barcode is Poisson. What you observe is
the number of *occupied* barcodes, and that is what pins :math:`\lambda`:

.. math::

   \text{occupied} = S_\text{eff}\,(1 - e^{-\lambda}),
   \qquad M = S_\text{eff}\,\lambda,
   \qquad P(k>1 \mid k\ge 1) = \frac{1 - e^{-\lambda} - \lambda e^{-\lambda}}{1 - e^{-\lambda}}

``p_multi`` is the last of these and it is the one that matters: the fraction of MIGs that are
really two or more molecules pooled. Their consensus is a mixture of templates, and no amount of
over-sequencing repairs it.

The familiar :math:`\binom{M}{2}/S` "expected collided pairs" is the small-:math:`\lambda` limit of
this and is badly wrong once the space is half full — which is exactly when somebody wants the
number.

.. list-table:: Two libraries, same tool
   :header-rows: 1
   :widths: 30 22 22 26

   * -
     - 12 nt, 4k molecules
     - 9 nt HIV Primer ID
     - 6 nt, 3k molecules
   * - nominal space
     - 16,777,216
     - 262,144
     - 4,096
   * - effective space
     - 16,763,376
     - 250,902
     - 4,092
   * - observed barcodes
     - 3,999
     - 125,236
     - 1,575
   * - occupancy
     - 0.02%
     - **49.9%**
     - **38.5%**
   * - MIGs holding >1 molecule
     - 0.01%
     - **30.6%**
     - **21.4%**
   * - molecules implied
     - 4,000
     - 173,482
     - 3,003

.. warning::

   Past 90% occupancy the estimate is declined rather than reported. :math:`S_\text{eff}` is
   inferred from the observed barcodes, so as the space fills the inversion collapses onto the
   observed count and would report "no collisions" for the most collided library there can be.
   The ``saturated`` flag says so instead.

Checking it against the reads
-----------------------------

The birthday number is a model. ``scripts/collision_check.py`` tests it against something that is
not: if two molecules sharing a barcode have *different sequences*, the reads say so directly.

.. code-block:: bash

   python scripts/collision_check.py --checkout out/S1.fq.gz --json out/checkout.json

Only collisions between molecules that differ are visible, so the raw count is corrected by the
probability that two random molecules differ — measured on the same data. On the HIV library that
gives 56.9% against the 30.6% predicted, **1.86×**.

The gap is expected and it has a known cause: :math:`\prod_j m_j` assumes the positions are
independent, and it is a *lower* bound on the true collision probability. Real synthesis correlates
neighbouring positions, so the true space is smaller and collisions more frequent than the
per-position product says. The discrepancy is the size of that correlation.

.. note::

   The observed figure is also biased upward by the read threshold: a barcode holding two molecules
   has roughly twice the reads, so it is over-represented among MIGs large enough to show a split.
   Read 1.86× as "meaningfully more than predicted, and in the direction the independence
   assumption implies", not as a calibrated factor.

The error budget
----------------

``checkout`` estimates the barcode error rate from the excess of barcode pairs at Hamming
distance 1. That is a measurement, and it can be checked against two predictions:

.. math::

   \varepsilon_\text{seq} = \left\langle 10^{-Q/10} \right\rangle
   \qquad \varepsilon_\text{pol} = \epsilon_\text{pol} \times n_\text{cycles}

Note: The first is the mean of the **probabilities**, not :math:`10^{-\bar{Q}/10}`. The function is
convex, so the low-Q tail carries nearly all the error and averaging Q first hides it: half the
bases at Q40 and half at Q10 is a 5% error rate, not the 0.3% that "mean Q25" suggests.

When estimate and prediction disagree by more than 3×, ``checkout`` says so, because one of them is
wrong and the ratio says which way to look. On a 2-colour instrument suspect the nominal Phred
first — there are only ~4 distinct Q values and they are not error rates.

.. warning::

   **The distance-1 estimator has a working range.** It subtracts the coincidence expectation from
   the observed pair count, and once a barcode's :math:`3L` neighbours are themselves mostly real
   barcodes, that is a small difference of two large numbers. Measured on simulated data at a
   known injected rate of 3·10⁻³:

   .. list-table::
      :header-rows: 1
      :widths: 20 20 30 30

      * - UMI
        - occupancy
        - estimated
        - fraction of truth
      * - 12 nt
        - 0.3%
        - 2.77·10⁻³
        - 0.92
      * - 10 nt
        - 4.3%
        - 2.64·10⁻³
        - 0.88
      * - 9 nt
        - 15.9%
        - 1.95·10⁻³
        - 0.65
      * - 8 nt
        - 49.7%
        - 7.0·10⁻⁴
        - 0.23
      * - 7 nt
        - 93.2%
        - 2.6·10⁻⁶
        - 0.001

   The collapse is always downward, so a crowded library reports too *little* barcode error and
   under-corrects. ``estimate_unreliable`` is set past 5% neighbourhood occupancy.

What is written out
-------------------

``checkout.barcode_space.tsv`` — one row per sample, every field above.
``checkout.umi_quality.tsv`` — the reported Phred histogram over barcode bases, which is the input
to :math:`\varepsilon_\text{seq}`.
``checkout.coverage.tsv`` — the MIG size histogram.
``checkout.umi_composition.tsv`` — per-position base usage, entropy, information, collision.

All four are also in ``checkout.json``, and ``notebooks/barcode_space.py`` draws them.
