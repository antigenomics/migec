Method
======

Why the numbers are what they are. Every page here exists because a first-pass answer was wrong and
the correction needed writing down -- these are not derivations for their own sake, they are the
measurements the defaults are set from.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - page
     - the question it answers
   * - :doc:`UMI statistics <umi_statistics>`
     - how long is my barcode really, and how much of it is the same base in every molecule
   * - :doc:`Barcode space <barcode_space>`
     - is the barcode long enough for this library, and how many molecules collided
   * - :doc:`Grouping accuracy <grouping>`
     - do the reads assigned to one molecule belong to one molecule
   * - :doc:`Fragmented libraries <fragmented>`
     - the reads of one 10x molecule are not co-terminal, so what is a consensus of them
   * - :doc:`The RT/PCR floor <quality_floor>`
     - the highest quality a consensus is allowed to claim, and why it is not Q60
   * - :doc:`Permutation nulls <nulls>`
     - what "significant" means for a split, when both margins have to be preserved
   * - :doc:`Validation <validation>`
     - the spike-in metric, and what it says the pipeline gets right

.. toctree::
   :maxdepth: 1

   umi_statistics
   barcode_space
   grouping
   fragmented
   quality_floor
   nulls
   validation
