Reference
=========

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - page
     - what is in it
   * - :doc:`File formats <formats>`
     - the inter-stage contract: tags, ``.mig`` buckets, every TSV column
   * - :doc:`Speed and memory <performance>`
     - what each stage costs, what it scales with, and how the thread invariance is asserted
   * - :doc:`Pipelines <nextflow>`
     - the nf-core-style modules with their ctDNA and AIRR subworkflows, the SLURM templates, and
       the three things that go wrong in a pipeline
   * - :doc:`Roadmap <roadmap>`
     - what is done, what is not, and what the version numbers mean

.. toctree::
   :maxdepth: 1

   formats
   performance
   nextflow
   roadmap

License
-------

MIGEC is licensed under GPL-3.0-or-later. The authoritative text is
`LICENSE <https://github.com/antigenomics/migec/blob/master/LICENSE>`_.
ARDA is a separate package with its own license; pipeline interoperability does
not make their licenses interchangeable.
