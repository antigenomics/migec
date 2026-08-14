subsample -- a smaller library that is still a library
======================================================

A smaller library that is still a library.

.. code-block:: bash

   migec subsample out/S1.fq.gz -o small.fq.gz --keep 0.5

Keeps **all the reads of a fraction of the barcodes**, selected by hashing. One streaming pass, no
sort, no memory.

What it is for
--------------

Sampling *reads* is trivial — ``seqtk sample``, ``head``, an awk one-liner. Sampling **molecules**
is not, because the reads of one molecule are scattered through the file and nothing in a raw FASTQ
says which they are. That is the whole job:

* **Equalise molecule counts across samples.** A cohort sequenced to different depths cannot be
  compared on anything that scales with molecules -- diversity, clonality, a per-target count --
  until every sample holds the same number. Down to the smallest sample's count, by molecule.
* **Saturation and rarefaction curves.** How many molecules would another lane buy? The answer is a
  curve over *molecule* fractions, and each point has to be a real library with its MIG size
  distribution intact, not a thinned read file.
* **Cheap iteration on a big run.** 1% of the barcodes with all of their reads runs in seconds and
  behaves like the library, so a pattern, a threshold or a whole pipeline can be checked before
  committing a lane to it.
* **Test fixtures**, which is the same property used for a different purpose.

The selection is deterministic and nested, so a curve computed at 1, 5, 10 and 50% is computed on
subsets of one another rather than on four independent draws.

Why not a fraction of the reads
-------------------------------

.. warning::

   **Never: Never subsample reads.** At 16 reads per molecule, keeping 0.5% of the *reads* gives
   molecules seen once each: the MIG size distribution is destroyed and every consensus collapses
   to a single read. The file still looks like a FASTQ, the pipeline still runs, and every fixture
   built from it silently tests a library nobody has.

   **Never: Nor the first N barcodes.** A UMI with 100 reads appears in the first thousand reads about a
   hundred times more often than a singleton, so first-appearance order oversamples large MIGs —
   destroying the very distribution the fixture exists to show.

``tests/synthetic/test_subsample.py`` runs the read-sampling comparison rather than asserting it:
at the same rate, barcode-sampling holds the mean MIG size to within 15% while read-sampling
collapses it below half.

Measured on the HIV Primer ID library: the shipped fixture keeps **15.77 reads per barcode against
the full library's 16.05**.

The selection
-------------

``splitmix64(packed barcode) % 10000 < keep × 100``.

- **Deterministic** — the same reads on any machine, so a fixture is reproducible from its
  definition rather than from a copy.
- **Nested** — 5% is a subset of 10%, because the test is on the same hash. A fixture shrinks
  without being re-derived.
- **Written out in** ``include/migec/subsample.hpp`` so the selection can be reproduced by anything,
  in any language. Not blake2b: this needs no dependency, and the only property asked of it is that
  the low bits are uncorrelated with barcode content.

Cells are kept whole
--------------------

When the reads carry a cell barcode the hash is taken on the **cell**, so a kept cell keeps every
molecule in it. Hashing cell+UMI together would sample molecules independently and give a fixture
of thousands of cells holding a handful of molecules each — the read-sampling mistake wearing a
different hat. ``--by-umi-only`` overrides it.

.. code-block:: text

   read  110,349
   kept  5,285 (4.8%) in 1,126 barcodes
         4.69 reads per barcode -- the same distribution as the input, which is the point
