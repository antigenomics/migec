UMI statistics
==============

Everything on this page comes out of the same pass as ``migec checkout``, and everything is written
to a TSV before it is plotted — so any figure can be redrawn from a committed table, and no number
lives only inside a C++ object.

Coverage histogram
------------------

``checkout.coverage.tsv`` bins molecules by reads-per-UMI in powers of two — MIGEC's 17 bins, kept
so published figures stay comparable::

   sample_id	mig_size	reads	units
   S1	1	20	20
   S1	2	115	45
   S1	4	80	16

``units`` is distinct UMIs, ``reads`` is reads. Both are needed and they tell different stories:
small MIGs are the *majority of UMIs* but a *minority of reads*, which is why a histogram weighted
by UMIs and one weighted by reads look nothing alike.

Two summary numbers are reported per sample:

``mean_reads_per_umi``
   Over-sequencing. Below about 5 most molecules are seen once and consensus assembly has nothing
   to work with.

``reads_in_migs_ge5``
   The fraction of reads sitting in MIGs of at least five reads — monotone in depth and directly
   interpretable.

.. note::

   MIGEC decided over-sequencing with a "is there a peak in the smoothed histogram" test. That test
   *inverts* on deeply sequenced libraries: the error-child fraction grows monotonically with depth,
   so a massively over-sequenced library gets classified as not over-sequenced. It is not
   reproduced here.

Base composition, entropy and information
-----------------------------------------

``checkout.umi_composition.tsv`` has one row per UMI position::

   sample_id	position	A	C	G	T	entropy_bits	information_bits	collision
   S1	0	0.220000	0.306667	0.300000	0.173333	1.962861	0.037139	0.262489

``entropy_bits``
   Shannon entropy :math:`H_j = -\sum_a p_j(a)\log_2 p_j(a)`. 2.0 for a uniform position.

``information_bits``
   :math:`2 - H_j`. This is the letter height in a sequence logo: 0 for a uniform position, 2 for a
   fixed one. Summed over positions it is ``total_information``, the number of bits the UMI is
   *wasting*.

``collision``
   :math:`m_j = \sum_a p_j(a)^2`, the probability two independent draws agree at this position.
   0.25 when uniform.

Which entropy, and why it matters
---------------------------------

.. warning::

   A logo draws Shannon entropy. **Collision arithmetic must not use it.**

The probability that two independent molecules receive the same UMI is

.. math::

   P_\text{collision} = \sum_u p_u^2 = \prod_j m_j

— the Rényi entropy of order 2, not order 1. Since :math:`H_2 \le H_1`, a Shannon-derived barcode
space is always **larger** than the true one, so using it *overestimates* the usable space and
*underestimates* collisions. That is the direction that silently merges distinct molecules into
one, and the error is invisible downstream.

Two numbers follow, and both are reported per sample:

.. math::

   L_\text{eff} = -\sum_j \log_4 m_j
   \qquad
   S_\text{eff} = \frac{1}{\prod_j m_j}
   \qquad
   E[\text{collisions}] \approx \frac{M^2}{2}\prod_j m_j

``effective_length`` is what your barcode is *worth* in bases. A 12 nt UMI whose first eight
positions are fixed has an effective length of 4 and a usable space of 256 — and will collide
constantly. The nominal length tells you nothing on its own.

.. note::

   It is measured over the **whole barcode**, cell then UMI, because that is what the counters are
   keyed on: a molecule is sample + cell + UMI. So on a single-cell run compare it against
   ``barcode_length`` (26 nt for 10x), never against ``umi_length`` (10 nt) — all three lengths are
   columns of ``checkout.summary.tsv`` for exactly that reason. On a bulk library they coincide.

.. note::

   Position independence is itself an assumption; oligo synthesis produces position-correlated
   bias, so :math:`\prod_j m_j` is a *lower bound* on the true collision rate. Measure it where it
   matters.

Count correction
----------------

For each pair of UMIs one substitution apart, three hypotheses are weighed.

**A sequencing miscall.** A miscall lands on one specific alternative base, so the rate per
neighbour is :math:`\varepsilon/3`, not :math:`\varepsilon`, and the child's size follows a
zero-truncated Poisson at :math:`c_\text{parent}\varepsilon/3` — truncated because a child with
zero reads is never observed and must not carry probability mass.

**A polymerase error.** Under a branching process the child's share :math:`f` of the family has
density :math:`\propto 1/f^2` (Luria–Delbrück). This component is the one a sequencing-only model
misses, and missing it is the dominant residual error in UMI counting: a substitution introduced in
PCR cycle 1–3 is present in roughly 50/25/12 % of the descendants and carries **high** quality in
every read, so a Poisson on the sequencing rate assigns it almost no probability and it survives as
a spurious second molecule. It also explains why MIGEC merged children below 10 % of their parent
and MAGERI below 1/20 — both far above anything :math:`\varepsilon/3` predicts.

**Another real molecule.** The barcode simply belongs to a different molecule that happens to sit
one substitution away. Its probability is :math:`n \cdot P_\text{collision}`, and its read count is
drawn from the library's own MIG size distribution — which means the test adapts to sequencing
depth without another tunable.

A child is merged when the posterior exceeds ``min_posterior`` (0.95). Consequences worth knowing:

* **An isolated low-coverage UMI keeps its reads.** A molecule seen 3–5 times with no plausible
  parent is information. It is never discarded, and it is never quality-derated either — if it
  *were* an error child of some parent, all of its reads would be clean reads of that parent's
  sequence, so the consensus would be right and only the molecule count wrong.
* **A neighbour of comparable size is not merged.** No error turns 10 000 reads into 9 000.
* **Reads are always conserved.** Correction moves reads between barcodes; it never deletes them.

The error rate itself is estimated from the data, from the excess of 1-mismatch neighbours over
what independent draws would produce:

.. math::

   E[D_1](\varepsilon) = \binom{n}{2}P_\text{collision}
     + 3L\sum_i \left(1 - e^{-c_i\varepsilon}\right)
     + 3L\sum_i \left(1 - e^{-c_i\varepsilon}\right)^2

.. note::

   The third term is the one that is easy to forget. Two children of the *same* parent that differ
   at the same position by different bases are themselves at Hamming distance 1, and counting them
   as parent–child pairs inflates the estimate by up to 2×.

Molecule counts
---------------

Two molecules that draw the *same* UMI and carry the *same* sequence are invisible to any method.
The observed molecule count is therefore biased low, and inverting the Poisson occupancy recovers
the estimate:

.. math::

   \hat{M} = S_\text{eff}\cdot -\ln\left(1 - \frac{M_\text{obs}}{S_\text{eff}}\right)

.. warning::

   Above 90 % occupancy migec **declines to estimate** and sets ``saturated``. :math:`S_\text{eff}`
   is itself estimated from the observed barcodes, so at saturation it collapses onto
   :math:`M_\text{obs}` and the formula would report "no collisions" for the most collided library
   possible.

Saturation does not disable correction — MIGEC switched it off entirely in this regime, on a gate
with no statistical meaning. Here the collision prior makes correction self-limiting on its own,
and the library is flagged instead.
