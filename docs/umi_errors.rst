Barcode error against depth
===========================

How often does a barcode base get called wrong? ``refine`` answers it twice, from two different
things in the same data, and reports both — because the two fail in opposite directions and a
single number would hide which failure you are looking at.

The first is the **distance-1 excess**: count pairs of observed barcodes one substitution apart,
subtract what independent draws would produce, and invert the remainder. That is
``estimated_error``, and it is described in :doc:`umi_statistics`.

The second is on this page. It reads the errors off the **children correction actually found**.

Two estimators, one number
--------------------------

A parent barcode carrying :math:`c` reads offered the instrument and the polymerase :math:`c L`
barcode bases to get wrong. What they got wrong is sitting in the table as that parent's error
children, so the same :math:`\varepsilon` falls out of the row two ways:

.. math::

   u(c) &= 3L \left( 1 - e^{-c\varepsilon/3} \right)
        &&\Rightarrow\quad \varepsilon = -\frac{3}{c}\ln\!\left(1 - \frac{u}{3L}\right) \\
   r(c) &= c\,L\,\varepsilon
        &&\Rightarrow\quad \varepsilon = \frac{r}{c\,L}

where :math:`u` is the **distinct child barcodes** per parent and :math:`r` is the **reads in those
children** per parent. Both are columns of ``<sample>.umi_errors.tsv``, one row per exact parent
depth.

The two differ in one respect that matters. A barcode has exactly :math:`3L` neighbours one
substitution away, so :math:`u` **saturates**: past the point where a parent has spawned all of
them, more errors add no new children and the estimate bends down. :math:`r` counts reads, of which
there is no ceiling, so it does not bend. Where the two curves separate on the figure is where this
library's barcode neighbourhood filled up — measured, not predicted.

.. note::

   ``error_from_variants`` is left **blank** past saturation rather than reported as a small
   number. Inverting a full neighbourhood returns "no errors" for the most error-ridden library
   there can be, and a blank cell is harder to misread than a zero.

What it is worth, and where it stops
------------------------------------

Both estimators are bounded by the merges ``correct_umis`` actually made. Neither is
saturation-free, and this table must not be read as though it were. Measured against a known
injected rate on simulated libraries, as a fraction of the truth:

.. list-table::
   :header-rows: 1
   :widths: 30 14 14 14 14 14

   * - occupancy
     - 0.2%
     - 2.3%
     - 9.8%
     - 33%
     - 100%
   * - distance-1 excess
     - 0.97
     - 0.96
     - 0.76
     - 0.45
     - 0.001
   * - from the children
     - **0.99**
     - **0.95**
     - **0.88**
     - **0.62**
     - 0.00

So the children estimate is the better of the two wherever either works, and at 0.2% occupancy on
a library sequenced 25 deep it lands within 1% of the injected rate. ``tests/synthetic/
test_umi_errors.py`` is what holds those numbers; the bounds there are loose on purpose, because
the claim is the trend and not a third decimal place.

.. warning::

   At 100% occupancy **both** go to zero, and for the same reason: on a full barcode space
   ``correct_umis`` refuses to merge — correctly, because a distance-1 neighbour there is more
   likely a real molecule than an error child. The ``saturated`` flag in the report is what says
   the answer is a floor. Read the flag. Do not read this table instead of it.

Read it at depth
----------------

A child whose parent was never sequenced cannot be merged into anything, so it is never counted.
At 1–3 reads per UMI that is about 80% of all barcode errors (:doc:`grouping`), which makes the
shallow end of this table a lower bound rather than a measurement.

That is why the report gives two figures:

* ``error_at_depth`` — restricted to parents seen at least ``error_depth`` (10) times, where
  correction is close to complete. This is the number to quote, and ``error_phred`` is the same
  thing as a Phred so it can be put beside the barcode's own reported quality.
* ``error_from_children`` — every depth, stated as the lower bound it is.

.. code-block:: text

   barcode error   9.73e-04 per base, estimated from the distance-1 excess
                   9.89e-04 from the reads in the children of molecules seen >= 10 times -- Q30
                   9.98e-04 over all depths, which is a LOWER bound: a child with no sequenced
                   parent cannot be found

The library there had 1e-3 injected. Three routes to the same number, agreeing to within 3%, is
what makes any of them believable — no one of them is the reference for the others.

The figures
-----------

``migec plot`` draws two panels from this table.

``umi_error_children`` puts the distinct children and the reads in them against the parent's depth,
with :math:`3L` as a dashed ceiling. Both series climb with :math:`c`; only one of them can climb
forever.

``umi_error_rate`` puts the two implied error rates against the same axis, with what ``refine``
reports drawn across them. The y axis is a log error rate, so one decade is exactly ten Phred and
1e-3 is Q30 — which is the comparison the panel exists for.

.. note::

   Both panels are drawn with **points, never lines**. One row is one exact depth, and past the
   head of the distribution most depths hold a handful of parents, so a line would render integer
   quantisation as structure and bridge gaps in the support where nothing was observed at all. The
   same correction was made to ``mig_size_spectrum``, where a line through one-molecule-per-size
   rows was drawing the :math:`y = x` diagonal as the most prominent feature of the figure.

The table
---------

``<sample>.umi_errors.tsv``, one row per distinct parent depth:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - column
     - meaning
   * - ``parent_reads``
     - reads carried by the parent, exact rather than binned
   * - ``parents``
     - surviving barcodes at that depth. Weight by this, or filter on it
   * - ``child_barcodes``
     - distinct barcodes merged into those parents
   * - ``child_reads``
     - reads carried by those children
   * - ``children_per_parent``
     - :math:`u`
   * - ``reads_per_parent``
     - :math:`r`
   * - ``neighbours``
     - :math:`3L`, the saturation ceiling. Constant down the column
   * - ``error_from_variants``
     - :math:`\varepsilon` from :math:`u`; ``.`` past saturation
   * - ``error_from_reads``
     - :math:`\varepsilon` from :math:`r`
   * - ``phred_from_reads``
     - the same, as a Phred
   * - ``estimate``
     - ``error_at_depth``, repeated so the panel can draw it. Constant down the column

``neighbours`` and ``estimate`` are constant on every row deliberately. The panels draw them as
reference lines, and a figure that needs a value its own table does not carry is a figure that will
one day disagree with the report.
