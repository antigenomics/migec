"""Reading the `ci/` fixtures.

Both fixtures are `checkout` OUTPUT -- a per-sample FASTQ whose barcodes live in SAM-style header
tags and whose barcode bases have already been trimmed out of the sequence. So every test here
starts by reading those tags, and the two that want to exercise `checkout` itself have to put the
barcode back into a read first. Each module says what that costs it.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import Iterator


def records(path: str | Path) -> Iterator[tuple[str, dict[str, str], str, str]]:
    """Yield `(name, tags, sequence, quality)` from a gzipped FASTQ written by a migec stage.

    Note: a tag value is split off at the SECOND colon, never with a plain `split(":")`. `QX` and
    `CY` are quality strings and Phred 25 is `:`, so a naive split tears every second record in
    half -- and it tears the ones with bad barcode quality, which are exactly the reads the
    matcher's acceptance bar is being tested on.
    """
    with gzip.open(path, "rt") as fh:
        while True:
            header = fh.readline()
            if not header:
                return
            sequence = fh.readline().rstrip("\n")
            fh.readline()
            quality = fh.readline().rstrip("\n")
            fields = header.rstrip("\n").lstrip("@").split()
            tags = {}
            for field in fields[1:]:
                key, _, value = field.split(":", 2)
                tags[key] = value
            yield fields[0], tags, sequence, quality
