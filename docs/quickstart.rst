Quickstart: a first molecular consensus
=======================================

This tutorial uses four synthetic reads from one molecule. It checks installation and
stage outputs; it is not biological validation. Install MIGEC as described in :doc:`installation`.

Create a tiny library
---------------------

The first eight bases are the UMI; the remaining bases are the observed payload.

.. code-block:: bash

   python - <<'PYTHON'
   from pathlib import Path
   sequence = "ACGTACGT" + "GCTAGCTAGGCTAACCGGTTAACCGGTTACGT"
   Path("reads.fq").write_text("".join(
       f"@read{i}\n{sequence}\n+\n{'I' * len(sequence)}\n" for i in range(4)))
   PYTHON

Extract, correct and assemble
-----------------------------

.. code-block:: bash

   migec checkout reads.fq --bc-pattern '^NNNNNNNN' --sample S1 -o out/
   migec refine out/S1.fq.gz -o ref/
   migec assemble ref/S1.fq.gz -o cons/

Inspect ``cons/S1.consensus.fq.gz`` and the stage JSON/TSV reports. The four input
reads share one UMI and identical payload, so expect one consensus molecule with
four supporting reads. A consensus count measures molecules, not input reads.

For your own library, replace the synthetic layout with the protocol's actual
:doc:`layout <layouts>`. Continue to :doc:`examples` for paired or cell-barcoded input,
or :doc:`formats` to interpret the output tags and tables.
