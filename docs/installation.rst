Installation
============

From PyPI
---------

.. code-block:: bash

   pip install migec

Wheels are built for CPython 3.10–3.13 on Linux x86-64 and macOS arm64. There are no Windows
wheels: this is a pipeline tool and nobody has asked. A source build on Windows configures, but is
not tested in CI.

Optional extras
---------------

.. code-block:: bash

   pip install "migec[notebooks]"   # polars + marimo, for notebooks/

The pipeline requires ``typer`` and ``seqtree``. Every stage writes plain TSV with the
standard library, and :doc:`migec plot <plots>` draws those tables with **gnuplot**, which is not a
Python package — install it from your package manager (``brew install gnuplot``,
``apt install gnuplot-nox``). Without it ``migec plot`` still writes the ``.gp`` scripts, so the
figures can be drawn anywhere.

**samtools** is the second external tool, and it is needed only to read a BAM, SAM or CRAM
(:doc:`Bring your own UMI <byo_umi>`). ``conda install -c bioconda samtools`` or
``brew install samtools``. A FASTQ run never calls it.

``seqtree`` owns sequence-neighbour search, including one-mismatch searches. It is installed
automatically; the historical ``migec[seqtree]`` extra remains a compatible alias. MIGEC owns
barcode-error likelihoods, collision accounting, quality and payload evidence, and merge policy.
A candidate within an edit budget is not by itself evidence that two molecules should merge.

From source
-----------

.. code-block:: bash

   git clone https://github.com/antigenomics/migec && cd migec
   bash setup.sh

``setup.sh`` creates ``.venv`` with uv, builds the extension in editable mode, and then asserts
that ``migec._core`` imports. That last step matters: without it a failed C++ build looks like a
successful install and only fails much later, in the middle of a long run.

Native source builds fetch the immutable seqtree 1.3.4 source archive and verify its SHA256.
For an offline build, provide that source explicitly with CMake's
``-DFETCHCONTENT_SOURCE_DIR_SEQTREE=/path/to/seqtree-source``. Wheel installation does not
fetch native source. Development editable installs use ``setup.sh``'s
``--no-build-isolation`` recipe so rebuild tools remain available in the isolated environment.

Running the tests
-----------------

.. code-block:: bash

   cmake -S . -B build -DMIGEC_TESTS=ON && cmake --build build -j
   ctest --test-dir build --output-on-failure
   python -m pytest tests/unit tests/synthetic -q
