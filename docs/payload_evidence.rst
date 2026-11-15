Payload draft evidence
======================

This explanation specifies the bulk compatibility safeguard. A draft is an
observed first read, not an error-free template or a consensus. Sequence and
quality are selected together, so later reads cannot replace only one half.

Observation model
-----------------

For compared base ``i``, let ``m_i`` be 1 for a mismatch and 0 for a match.
Let ``e_1i`` and ``e_2i`` be the Phred-derived error probabilities of the two
observations, capped at 0.75 for an uninformative four-base observation.
The probability that two observations of the same template disagree is
``delta_i = e_1i + e_2i - 4 e_1i e_2i / 3``. The existing ``payload_error``
parameter floors this **pair disagreement**, not each individual error.
Unknown qualities use that historical floor, an assumption rather than evidence.

For latent template mismatch fraction ``theta``, the observation mismatch
probability is ``p_i(theta) = delta_i + theta (1 - 4 delta_i / 3)``.
The same-template model uses ``theta=0``; the related-template boundary uses
``theta=payload_same_fraction`` (0.05). Each log likelihood sums
``m_i log(p_i) + (1-m_i) log(1-p_i)`` over comparable informative bases.
Zero bases, N bases and zero-information emissions are excluded.

The moment estimate is ``theta_hat = clip(sum(m_i-delta_i) /
sum(1-4 delta_i/3), 0, 1)``. A bulk incompatibility veto requires at least
eight informative bases, ``theta_hat > 0.05`` and a same/boundary likelihood
ratio below 0.1. These are the existing compatibility and refutation bounds.
Partially informative bases contribute their own emission probabilities and
contrast weights; a zero denominator supplies no template evidence.

The historical positive-evidence identical/random-template mixture and count
gates remain. This two-part rule is not a calibrated three-template posterior
or an independent estimate of false-merge rate. Reservoir clonality and residual
agreement use the same latent-fraction interpretation; cell-tagged comparisons
retain their preceding fixed-emission and observed-fraction behavior.

Storage and verification
------------------------

The public payload width remains the number of bases. Internal paired blocks
contain the bases followed by matching Phred bytes, with a separate stride.
The existing stable first-read selection copies the entire block atomically.
Resident compaction, spill reduction, rotated passes and wider-prefix
repartitioning preserve that pairing. Payload-quality spill files have their own
``payq1`` identity. Physical .mig format and read sequences are unchanged.
Memory budgets and reported table sizes include both halves of the block.

Neutral regressions cover low-quality draft errors, confidently different
templates, missing evidence, complete query candidates, cell scope, and paired
quality provenance across resident/spilled/rotated storage. Matched-data
qualification must separately report physical reads, correction assignments,
molecule groups and consensus records. More groups are not proof of more
independent molecules; supported earlier corrections and template separations
must both be traced.
