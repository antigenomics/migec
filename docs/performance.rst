Speed and memory
================

Both are reported on every run, not behind a flag. They are the two numbers that decide whether a
pipeline can be run at all, and a tool that prints them only when asked gets run once without
asking.

.. code-block:: text

   2.2 s (917,339 reads/s on 8 threads), peak RSS 204.0 MB of which UMI counters 44.4 MB

The same fields are in ``checkout.json`` as ``wall_seconds``, ``reads_per_second``, ``threads``,
``peak_rss_bytes`` and ``umi_memory_bytes``.

Threads
-------

``--threads/-t`` defaults to one per core. **The output is byte-identical whatever it is set to.**
Reads are matched in fixed-size chunks and the chunks are written back in input order, so ``-t``
changes the wall clock and nothing else — a demultiplexer whose output depended on its thread count
would produce results that could not be compared between runs.

Measured on 2 M single-end 115 nt reads over four barcode patterns, on an M-series laptop:

.. list-table::
   :header-rows: 1
   :widths: 15 25 30 30

   * - threads
     - wall clock
     - reads/s
     - peak RSS
   * - 1
     - 13.1 s
     - 152,822
     - 139 MB
   * - 2
     - 6.7 s
     - 298,581
     - 119 MB
   * - 4
     - 3.7 s
     - 545,211
     - 151 MB
   * - 8
     - 2.2 s
     - 917,339
     - 204 MB
   * - 16
     - 1.7 s
     - 1,178,648
     - 295 MB

Two things had to be true for that to scale, and neither is obvious:

**Compression runs on the workers, not on the writer.** zlib at its default level 6 compresses
random DNA at about **7 MB/s**. Read payload is close to incompressible, so leaving compression on
the serial path caps checkout at a fraction of what the matcher can do no matter how many threads
are matching. Each worker gzips its own chunk and the writer only appends bytes — concatenated
gzip members are a valid gzip stream (:rfc:`1952` §2.2), so the result is an ordinary ``.fq.gz``
that ``zcat`` and every reader accept.

**The default compression level is 1, not 6.** On random DNA level 1 runs at 137 MB/s for 13% more
bytes. Paying twenty times the CPU for a tenth off the file is not a trade anyone would make
deliberately.

There is also no transcendental in the scoring loop. The log-likelihood only ever depends on the
reported Phred and the size of the IUPAC set, both small integers, so the whole score function
tabulates into 1.2 kB. It was 90% of runtime before it did.

.. note::

   Non-scaling parts are the gzip *read* of the input, which is one thread by construction, and the
   ``fwrite`` of already-compressed blocks. At 16 threads on this machine the reader is the wall.

Memory
------

Two allocations matter, and they scale differently.

**Per-worker buffers** are bounded by ``chunk_reads × threads``, about 5 MB per thread. They do not
grow with the input, which is why a 2 M-read run and a 2 G-read run have the same buffer footprint.

**The UMI counters** grow with the number of *distinct* UMIs, and are the reason this section
exists. They are a sorted ``(key, count)`` array with a bounded append buffer, not a hash map:

.. list-table::
   :header-rows: 1
   :widths: 45 25 30

   * - structure
     - bytes per distinct UMI
     - at 4·10⁸ UMIs
   * - ``unordered_map<uint64_t, uint32_t>``
     - ~48
     - 19 GB
   * - sorted ``(key, count)`` array
     - **~22 measured**
     - **8.8 GB**

Four hundred million distinct UMIs is an ordinary NovaSeq output at five reads per molecule, so the
difference is the difference between a run fitting and not. Sorted order is not a side effect
either: it is what the range partition and the 1-substitution neighbourhood search both want, and
it turns the structure into a flat scan instead of a pointer chase.

The append buffer grows with the data rather than sitting at a fixed ceiling per sample. A fixed
buffer costs that ceiling for every sample whatever the sample holds, which on a 96-plex sheet is
gigabytes of empty space.

.. warning::

   8.8 GB still does not fit a laptop. The counters are **not yet partitioned**, so a full run is
   held in one piece. The fix is the range partition — process one bucket of the barcode space at a
   time, so the counter only ever holds 1/2\ :sup:`bits` of the library — and it lands with ``.mig``
   bucket output in M2. Until then ``checkout`` warns when the counters pass 1 GB rather than
   letting you find out from the OOM killer.

Benchmarks
----------

``tests/benchmark/`` holds the regressions. They are off by default because CI runners vary by more
than the thing being measured:

.. code-block:: bash

   RUN_BENCHMARK=1 python -m pytest tests/benchmark -q -s

The thresholds are deliberately loose — they exist to catch a 10× regression, such as a
transcendental finding its way back into the scoring loop or compression migrating back onto the
serial path, not to police a 10% one.
