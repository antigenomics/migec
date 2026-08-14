"""Which `.mig` buckets a path names.

`checkout --mig` and `refine` both write `<sample>.<bbb>.mig`, one file per range-partition
bucket, and `refine` and `assemble` both take them. One bucket file names the whole partition:
a single bucket is a corner of the barcode space, and running a stage on it alone would silently
process a fraction of the sample rather than fail.
"""

from __future__ import annotations

from pathlib import Path


def mig_buckets(path: str | Path) -> list[str]:
    """The buckets `path` names, or an empty list when it is a FASTQ.

    A `.mig` file brings its siblings; a directory holding exactly one sample's buckets is that
    sample. Never: a directory holding two samples is refused by name rather than assembled
    together -- a UMI repeats across samples by design, so grouping them invents molecules that
    never existed, and nothing downstream can tell.
    """
    p = Path(path)
    if p.suffix == ".mig":
        sample = p.name.split(".")[0]
        return sorted(str(f) for f in p.parent.glob(f"{sample}.*.mig"))
    if not p.is_dir():
        return []
    by_sample: dict[str, list[str]] = {}
    for f in sorted(p.glob("*.mig")):
        by_sample.setdefault(f.name.split(".")[0], []).append(str(f))
    if len(by_sample) > 1:
        names = ", ".join(sorted(by_sample))
        raise ValueError(
            f"{p} holds buckets for {len(by_sample)} samples ({names}). This is a per-sample "
            f"stage: point it at one sample's buckets, e.g. {p}/{sorted(by_sample)[0]}.000.mig, "
            f"which brings the rest of that sample with it"
        )
    return next(iter(by_sample.values()), [])
