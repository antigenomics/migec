"""Which `.mig` buckets a path names.

`checkout --mig` and `refine` both write `<sample>.<bbb>.mig`, one file per range-partition
bucket, and `refine` and `assemble` both take them. One bucket file names the whole partition:
a single bucket is a corner of the barcode space, and running a stage on it alone would silently
process a fraction of the sample rather than fail.
"""

from __future__ import annotations

from pathlib import Path


def _sample_of(name: str) -> str:
    """The sample id in `<sample>.<bbb>.mig`, where the id itself may contain periods.

    Never: not `name.split(".")[0]`. A sample id with a period in it -- `S1.rep2` -- would then
    read as sample `S1`, and two samples whose ids share a prefix would be collected into one run.
    The suffix is fixed: three digits and `.mig`, so the id is everything before them.
    """
    stem = name[: -len(".mig")] if name.endswith(".mig") else name
    head, _, bucket = stem.rpartition(".")
    return head if head and bucket.isdigit() else stem


def mig_buckets(path: str | Path) -> list[str]:
    """The buckets `path` names, or an empty list when it is a FASTQ.

    A `.mig` file brings its siblings; a directory holding exactly one sample's buckets is that
    sample. Never: a directory holding two samples is refused by name rather than assembled
    together -- a UMI repeats across samples by design, so grouping them invents molecules that
    never existed, and nothing downstream can tell.
    """
    p = Path(path)
    # Never: `glob` on a name built from a sample id treats `[`, `*` and `?` as pattern syntax, so
    # an id carrying one selects the wrong files -- or none, silently. The bucket suffix is always
    # three digits and the extension is fixed, so the sibling test is a plain string comparison.
    def siblings(directory: Path, sample: str) -> list[str]:
        out = [
            f
            for f in directory.iterdir()
            if f.suffix == ".mig" and _sample_of(f.name) == sample
        ]
        return sorted(str(f) for f in out)

    if p.suffix == ".mig":
        return siblings(p.parent, _sample_of(p.name))
    if not p.is_dir():
        return []
    by_sample: dict[str, list[str]] = {}
    for f in sorted(p.iterdir()):
        if f.suffix == ".mig":
            by_sample.setdefault(_sample_of(f.name), []).append(str(f))
    if len(by_sample) > 1:
        names = ", ".join(sorted(by_sample))
        raise ValueError(
            f"{p} holds buckets for {len(by_sample)} samples ({names}). This is a per-sample "
            f"stage: point it at one sample's buckets, e.g. {p}/{sorted(by_sample)[0]}.000.mig, "
            f"which brings the rest of that sample with it"
        )
    return next(iter(by_sample.values()), [])
