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

   pip install "migec[seqtree]"     # whitelist lookup at >=2 substitutions, N-wildcard barcodes
   pip install "migec[notebooks]"   # polars + marimo, for notebooks/

The pipeline itself has one runtime dependency, ``typer``. Every stage writes plain TSV with the
standard library, and :doc:`migec plot <plots>` draws those tables with **gnuplot**, which is not a
Python package — install it from your package manager (``brew install gnuplot``,
``apt install gnuplot-nox``). Without it ``migec plot`` still writes the ``.gp`` scripts, so the
figures can be drawn anywhere.

``seqtree`` is deliberately *not* a core dependency. Every hot path here searches a fixed-length
barcode at at most one substitution, and enumerating the ``3L`` neighbours in a hash table beats a
trie by orders of magnitude at that shape. seqtree earns its place only where enumeration cannot
express the query.

From source
-----------

.. code-block:: bash

   git clone https://github.com/antigenomics/migec && cd migec
   bash setup.sh

``setup.sh`` creates ``.venv`` with uv, builds the extension in editable mode, and then asserts
that ``migec._core`` imports. That last step matters: without it a failed C++ build looks like a
successful install and only fails much later, in the middle of a long run.

Running the tests
-----------------

.. code-block:: bash

   cmake -S . -B build -DMIGEC_TESTS=ON && cmake --build build -j
   ctest --test-dir build --output-on-failure
   python -m pytest tests/unit tests/synthetic -q
