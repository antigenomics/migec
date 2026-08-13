Speed and memory
================

Both are reported on every run, not behind a flag. They are the two numbers that decide whether a
pipeline can be run at all, and a tool that prints them only when asked gets run once without
asking.

.. code-block:: text

   2.2 s (903,599 reads/s) = 1.5 s matching on 8 threads + 0.7 s UMI statistics, serial
   peak RSS 131.0 MB of which UMI counters 21.2 MB

The same fields are in ``checkout.json`` as ``wall_seconds``, ``match_seconds``,
``reads_per_second``, ``threads``, ``peak_rss_bytes`` and ``umi_memory_bytes``.

Two clocks, because they scale differently
------------------------------------------

``match_seconds`` is the demultiplexing driver: read a chunk, match it, compress it, append it.
That is the part ``--threads`` speeds up, and it scales with the number of *reads*.

``wall_seconds`` also covers the per-sample statistics — the coverage histogram, the composition,
and the count correction — which run once at the end, on one thread, and scale with the number of
*distinct UMIs* at roughly 1.5–2 µs each. On a shallow library that is most of the run.

They are reported separately because a single number would hide which one to attack, and because a
throughput figure that stopped at the driver would be measuring the matcher rather than checkout.

Threads
-------

``--threads/-t`` defaults to one per core. **The output is byte-identical whatever it is set to.**
Reads are matched in fixed-size chunks and the chunks are written back in input order, so ``-t``
changes the wall clock and nothing else — a demultiplexer whose output depended on its thread count
would produce results that could not be compared between runs.

Measured on 2 M single-end 129 nt reads over four barcode patterns at 4 reads per molecule — the
corpus ``tests/benchmark/`` builds — on an M-series laptop:

.. list-table::
   :header-rows: 1
   :widths: 10 15 18 15 20 22

   * - threads
     - wall clock
     - reads/s
     - matching
     - matching reads/s
     - peak RSS
   * - 1
     - 9.9 s
     - 202,717
     - 9.2 s
     - 217,506
     - 60 MB
   * - 2
     - 5.5 s
     - 362,486
     - 4.8 s
     - 412,855
     - 78 MB
   * - 4
     - 3.3 s
     - 603,773
     - 2.6 s
     - 758,801
     - 98 MB
   * - 8
     - 2.2 s
     - 910,500
     - 1.5 s
     - 1,335,791
     - 136 MB
   * - 16
     - 1.9 s
     - 1,056,472
     - 1.2 s
     - 1,684,654
     - 217 MB

Measured 2026-08-13 by ``python scripts/benchmark_threads.py --reads 2000000 -o assets/``, which
writes ``assets/benchmark_threads.tsv``; ``migec plot assets/`` draws the figure from that table,
so the two cannot drift apart. Earlier prose quoted 1.18 M reads/s here, from a run whose corpus
sent every read to one sample of four -- the matcher still scored all four patterns, so the
matching column was sound, but the per-sample counters were not, and the memory figure was a
single sample's.

The serial 0.7 s of statistics is why the end-to-end column flattens at 16 threads while the
matching column is still scaling. Amdahl, and it is the next thing to fix rather than the thread
count.

Two things had to be true for the matching to scale, and neither is obvious:

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

   Non-scaling parts are the gzip *read* of the input, which is one thread by construction, the
   ``fwrite`` of already-compressed blocks, and the per-sample statistics. At 16 threads on this
   machine the statistics are the larger of the two walls.

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
