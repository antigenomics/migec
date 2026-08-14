"""BAM, SAM and CRAM input: convert once with samtools, run the stage on the result.

For a capture, exome, ctDNA or MRD library the UMI is in the index read, so it never reaches the
user as a FASTQ -- it arrives as an fgbio/Picard/xGen BAM carrying `RX` (and usually `QX`). That
file is a migec input in every respect except that migec could not open it, which is how a real
dataset ends up with "has UMIs, we cannot use them" written next to it (`SOURCES.md`, SRP578416).

Never: no htslib. `project/design-io-interop.md` ratified that out -- its own build system, its own
wheel problem -- and named `samtools` as the route in both directions. This is the reverse of a
conversion the repo already documents, so it uses the same tool. samtools is an external program
like gnuplot, not a Python dependency.

Never: not a FIFO, a real file. `refine` opens its input THREE times in one call (whitelist pass,
table pass, rewrite) and `assemble` sizes its range partition from `std::filesystem::file_size`,
which on a pipe falls back to the minimum bucket count and blows pass 2's memory budget. A named
pipe would break both, silently in the second case.
"""

from __future__ import annotations

import contextlib
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

# The tags the stages actually read out of the comment: `RX`/`CB` are the barcode, `QX`/`CY` its
# quality, `BC` the sample id. Carrying more would be noise in every record.
TAGS = "RX,QX,CB,CY,BC"

_INSTALL = "install it with `conda install -c bioconda samtools` or `brew install samtools`"


def is_alignment(path: str | Path) -> bool:
    """True for BAM, CRAM or SAM. Content first; the name only decides where there is no magic."""
    p = Path(path)
    if p.is_dir():
        return False
    try:
        with open(p, "rb") as fh:
            head = fh.read(4)
    except OSError:
        return False
    if head == b"CRAM" or head == b"BAM\x01":
        return True
    if head[:2] == b"\x1f\x8b":
        # A BGZF-wrapped BAM and a gzipped FASTQ have the same first two bytes, so the gzip check
        # alone cannot tell them apart -- inflate far enough to read the block's own magic.
        try:
            with gzip.open(p, "rb") as fh:
                return fh.read(4) == b"BAM\x01"
        except OSError:
            return False
    if head.startswith(b"@HD\t") or head.startswith(b"@SQ\t"):
        return True
    return p.suffix.lower() == ".sam"


@contextlib.contextmanager
def as_fastq(source: str | Path, workdir: str | Path, need_rx: bool = True):
    """Convert `source` to FASTQ in a temporary directory, yielding `(mate1, mate2 or None)`.

    The temporary directory is created inside `workdir` -- the output directory the user chose --
    rather than in the system temp, for the same reason `.assemble_buckets` and `.refine_spill`
    live there: it is the filesystem they have room on. It is removed when the stage returns.

    Cost: the FASTQ is roughly 4x the BAM on disk while the stage runs.
    """
    src = Path(source)
    if shutil.which("samtools") is None:
        raise ValueError(f"reading {src.name} needs samtools on PATH -- {_INSTALL}")
    work = Path(workdir)
    work.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(work), prefix=".bam2fq.") as tmp:
        d = Path(tmp)
        r1, r2, r0, singles = d / "R1.fq", d / "R2.fq", d / "R0.fq", d / "S.fq"
        # Never: on an ALIGNED file, collate is not optional. `samtools fastq -1/-2` pairs by
        # ADJACENCY, so on a coordinate-sorted BAM it writes one molecule's mate 1 beside
        # another's mate 2 -- and `assemble` matches mates strictly by position, so the wrong pair
        # is consensed and nothing downstream can tell. Sorting by name first is what makes
        # position-matching true.
        #
        # An UNALIGNED file cannot be coordinate-sorted -- there are no coordinates -- so collate
        # buys nothing there and costs the record order, which is the order refine rewrites its
        # reads in. `@SQ` is what tells the two apart, and it is the conservative way round: a
        # header we cannot read is treated as aligned.
        fastq = ["samtools", "fastq", "-n", "-T", TAGS,
                 "-1", str(r1), "-2", str(r2), "-0", str(r0), "-s", str(singles)]
        # collate's stderr goes to a file rather than a pipe: reading a pipe only after the second
        # process exits deadlocks if the first fills it, and there is nothing here to gain by
        # streaming it.
        log = d / "samtools.log"
        failed = 0
        with open(log, "w") as errfh:
            if _has_references(src):
                collate = subprocess.Popen(
                    ["samtools", "collate", "-u", "-O", "-T", str(d / "collate"), str(src)],
                    stdout=subprocess.PIPE, stderr=errfh,
                )
                convert = subprocess.run(
                    fastq + ["-"], stdin=collate.stdout, capture_output=True, text=True,
                    check=False,
                )
                if collate.stdout is not None:
                    collate.stdout.close()
                collate.wait()
                failed = collate.returncode
            else:
                convert = subprocess.run(
                    fastq + [str(src)], capture_output=True, text=True, check=False,
                )
        if failed or convert.returncode:
            raise ValueError(_why(src, log.read_text(errors="replace") or convert.stderr))

        paired = r1.stat().st_size > 0 and r2.stat().st_size > 0
        loose = [p for p in (r0, singles) if p.stat().st_size > 0]
        if paired and loose:
            raise ValueError(
                f"{src.name} holds both paired and unpaired reads "
                f"({', '.join(p.name for p in loose)} are not empty). Split it before running "
                f"migec: which of the two the run should describe is not ours to guess"
            )
        if paired:
            first, second = r1, r2
        elif loose:
            first, second = loose[0], None
        elif r1.stat().st_size > 0:
            first, second = r1, None          # mate 1 only, e.g. a filtered BAM
        else:
            raise ValueError(f"{src.name} converted to no reads at all")

        if need_rx:
            _require_rx(first, src)
        yield first, second


def _has_references(source: Path) -> bool:
    """True if the header declares any `@SQ`, i.e. the file is aligned. Unreadable counts as yes."""
    try:
        header = subprocess.run(
            ["samtools", "view", "-H", str(source)], capture_output=True, text=True, check=False,
        )
    except OSError:
        return True
    if header.returncode:
        return True
    return any(line.startswith("@SQ\t") for line in header.stdout.splitlines())


def _require_rx(fastq: Path, source: Path) -> None:
    """Fail now if the barcode is not there, rather than succeeding with zero molecules."""
    with open(fastq) as fh:
        head = fh.readline()
    if "RX:Z:" not in head:
        raise ValueError(
            f"{source.name} carries no RX tag: its first read converted to '{head.strip()}'. "
            f"refine and assemble group on RX, so without it the run would report zero molecules "
            f"and no error. If the UMI is in a separate index read, put it in RX first "
            f"(`fgbio AnnotateBamWithUmis`, `picard FastqToSam UMI_FASTQ=`); if the barcode is "
            f"still inline in the read, run `migec checkout` on this file instead"
        )


def _why(source: Path, stderr: str) -> str:
    """samtools' own last line, plus the hint for the one failure whose message is unhelpful."""
    last = next((ln for ln in reversed(stderr.strip().splitlines()) if ln.strip()), "")
    msg = f"samtools could not convert {source.name}: {last or 'no message'}"
    if "reference" in stderr.lower() or source.suffix.lower() == ".cram":
        msg += (
            " -- a CRAM needs its reference: set REF_PATH or REF_CACHE, or convert it yourself "
            "with `samtools fastq --reference <ref.fa>`"
        )
    return msg
