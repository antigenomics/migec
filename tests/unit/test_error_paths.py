"""What every stage does with input that is not a library.

A job killed mid-write, a download that stopped, a sheet pasted out of a spreadsheet: each of
these reaches a stage as bytes that are nearly right. The crash is not the failure that matters --
the failure that matters is the run that reads half a file, drops the rest and reports success,
because a consensus assembled from part of a molecule's reads is still a consensus and nothing
downstream can tell. Every path here was untested.

The stages are exercised through their Python entry points rather than the CLI, because that is
where the C++ exception crosses into Python: a message that never becomes a `MigecError` is a
`std::terminate` with no output at all.
"""

from __future__ import annotations

import gzip
import random
import re
import shutil
import struct
from pathlib import Path

import pytest

from tests.conftest import requires_core

pytestmark = requires_core

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def write_sheet(path: Path, sample_id: str = "S1", umi_len: int = 12) -> Path:
    path.write_text(f"{sample_id}\t{'N' * umi_len}{ADAPTER.lower()}\n")
    return path


def write_corpus(path: Path, n_umis: int = 300, reads_per_umi: int = 4, seed: int = 0) -> Path:
    """A whole library: every UMI carries all of its reads, never a sample of them."""
    rng = random.Random(seed)
    with gzip.open(path, "wt") as fh:
        for i in range(n_umis):
            umi = "".join(rng.choice("ACGT") for _ in range(12))
            payload = "".join(rng.choice("ACGT") for _ in range(60))
            for r in range(reads_per_umi):
                s = umi + ADAPTER + payload
                fh.write(f"@r{i}_{r}\n{s}\n+\n{'I' * len(s)}\n")
    return path


def empty_gzip(path: Path) -> Path:
    with gzip.open(path, "wt"):
        pass
    return path


def gzip_text(path: Path) -> str:
    """The file's contents, which also asserts it is a complete gzip member."""
    with gzip.open(path, "rt") as fh:
        return fh.read()


def chop(src: Path, dst: Path, tail_bytes: int = 40) -> Path:
    dst.write_bytes(src.read_bytes()[:-tail_bytes])
    return dst


@pytest.fixture(scope="module")
def library(tmp_path_factory):
    """One small library, checked out both ways: FASTQ, and `.mig` buckets.

    Module-scoped because three of the tests only need something well-formed to damage, and
    checking out 1,200 reads twice per test is the whole cost of this file.
    """
    from migec.checkout import run as checkout

    d = tmp_path_factory.mktemp("errpaths")
    reads = write_corpus(d / "reads.fq.gz")
    sheet = write_sheet(d / "bc.txt")
    checkout(reads, sheet, d / "co")
    mig = checkout(reads, sheet, d / "co_mig", mig=True)
    return {
        "dir": d,
        "reads": reads,
        "sheet": sheet,
        "tagged": d / "co" / "S1.fq.gz",
        "mig_dir": d / "co_mig",
        "mig_paths": [Path(p) for p in mig["mig_paths"]],
    }


# ------------------------------------------------------------------ an empty FASTQ


def test_an_empty_fastq_checks_out_into_a_file_the_next_stage_can_read(tmp_path):
    """Delete this and an empty input is free to write a zero-byte `S1.fq.gz`, which is not a
    gzip member at all -- the next stage then reports a truncated file, blaming corruption for a
    run that simply had no reads."""
    from migec.checkout import run as checkout

    empty = empty_gzip(tmp_path / "empty.fq.gz")
    summary = checkout(empty, write_sheet(tmp_path / "bc.txt"), tmp_path / "co")

    assert summary["total"] == 0
    assert summary["assigned"] == 0
    out = tmp_path / "co" / "S1.fq.gz"
    assert out.exists(), "the sample's file must exist even with nothing in it"
    assert gzip_text(out) == ""
    # And the sample still has a row in the summary table: an absent row and a zero row mean
    # different things to anything that joins on sample id.
    rows = (tmp_path / "co" / "checkout.summary.tsv").read_text().splitlines()
    assert len(rows) == 2 and rows[1].split("\t")[:3] == ["S1", "0", "0"]


def test_refine_refuses_an_empty_fastq_and_names_what_is_missing(tmp_path):
    """The one stage that raises here. Delete this and refine may start reporting zero molecules
    for a file with no barcodes in it, which is the same output an untagged input gives -- the
    user then looks for the missing molecules instead of for the missing `checkout` step."""
    from migec.refine import run as refine

    empty = empty_gzip(tmp_path / "empty.fq.gz")
    with pytest.raises(RuntimeError, match="RX"):
        refine(empty, tmp_path / "ref")


def test_assemble_and_subsample_report_an_empty_input_as_empty(tmp_path):
    """Both accept it. What must hold is that the outputs are complete gzip members holding zero
    records: an empty consensus file that cannot be opened turns "no reads" into "corrupt file"
    one stage later, and that is the diagnosis a user acts on."""
    from migec.assemble import run as assemble
    from migec.subsample import run as subsample

    empty = empty_gzip(tmp_path / "empty.fq.gz")

    a = assemble(empty, tmp_path / "asm", sample_id="S1")
    assert (a["reads"], a["groups"], a["molecules"]) == (0, 0, 0)
    assert gzip_text(tmp_path / "asm" / "S1.consensus.fq.gz") == ""
    # The molecule table is a header and nothing else, not an absent file.
    assert (tmp_path / "asm" / "S1.mig.tsv").read_text().startswith("cell\tumi\t")

    s = subsample(empty, tmp_path / "sub.fq.gz", keep_percent=50.0)
    assert (s["reads"], s["reads_kept"], s["barcodes"]) == (0, 0, 0)
    assert gzip_text(tmp_path / "sub.fq.gz") == ""


# ------------------------------------------------------------------ a truncated gzip


def test_a_truncated_gzip_raises_and_names_the_file(library, tmp_path):
    """The failure this exists for is a SHORT READ: zlib hitting the end of a chopped member and
    the stage treating it as end of input. Every read past the cut then vanishes with a summary
    that says the run succeeded, and the molecule counts are quietly wrong. Delete this and
    nothing distinguishes "the file ended" from "the file ends here".

    The message must carry the path because a pipeline runs four stages over dozens of files and
    "truncated record 1197" alone does not say which one to re-fetch.
    """
    from migec.assemble import run as assemble
    from migec.checkout import run as checkout
    from migec.refine import run as refine
    from migec.subsample import run as subsample

    raw = chop(library["reads"], tmp_path / "raw.fq.gz")
    tagged = chop(library["tagged"], tmp_path / "tagged.fq.gz")
    cases = {
        "checkout": lambda: checkout(raw, library["sheet"], tmp_path / "co"),
        "refine": lambda: refine(tagged, tmp_path / "ref"),
        "assemble": lambda: assemble(tagged, tmp_path / "asm"),
        "subsample": lambda: subsample(tagged, tmp_path / "sub.fq.gz", keep_percent=100.0),
    }
    for stage, call in cases.items():
        with pytest.raises(RuntimeError) as exc:
            call()
        message = str(exc.value)
        named = raw if stage == "checkout" else tagged
        assert str(named) in message, f"{stage} did not name the file: {message}"
        assert "truncat" in message, f"{stage} did not say what was wrong: {message}"


# ------------------------------------------------------------------ a malformed FASTQ


def test_a_record_missing_its_quality_line_is_refused_by_record_number(library, tmp_path):
    """A three-line record is what a `head -n` or a killed writer leaves behind. Accepting it
    shifts every following record by one line, so sequences are read as headers and the run
    produces barcodes that were never sequenced. Delete this and that becomes silent."""
    from migec.checkout import run as checkout
    from migec.subsample import run as subsample

    bad = tmp_path / "three_lines.fq"
    bad.write_text(f"@r0\n{'ACGT' * 20}\n+\n")

    for call in (lambda: checkout(bad, library["sheet"], tmp_path / "co"),
                 lambda: subsample(bad, tmp_path / "sub.fq.gz", keep_percent=100.0)):
        with pytest.raises(RuntimeError) as exc:
            call()
        message = str(exc.value)
        assert str(bad) in message
        assert re.search(r"record 1\b", message), message


def test_a_quality_string_shorter_than_the_sequence_is_refused(library, tmp_path):
    """The one malformation that parses. Every stage indexes quality by sequence position -- the
    consensus posterior, the barcode-error estimator, the UMI quality filter -- so a short quality
    string is an out-of-bounds read on the hot path. Delete this and the check that stops it can
    go with it, and what replaces the failure is whatever was next in memory."""
    from migec.assemble import run as assemble
    from migec.checkout import run as checkout
    from migec.refine import run as refine

    raw = tmp_path / "short_qual.fq"
    raw.write_text(f"@r0\n{'ACGT' * 20}\n+\n{'I' * 4}\n")
    tagged = tmp_path / "short_qual_tagged.fq"
    tagged.write_text(
        f"@r0 RX:Z:ACGTACGTACGT\tQX:Z:{'I' * 12}\tBC:Z:S1\n{'ACGT' * 20}\n+\n{'I' * 4}\n"
    )

    for path, call in ((raw, lambda: checkout(raw, library["sheet"], tmp_path / "co")),
                       (tagged, lambda: refine(tagged, tmp_path / "ref")),
                       (tagged, lambda: assemble(tagged, tmp_path / "asm"))):
        with pytest.raises(RuntimeError) as exc:
            call()
        message = str(exc.value)
        assert str(path) in message
        assert "quality" in message, message


# ------------------------------------------------------------------ a damaged .mig bucket


def _biggest_bucket(paths: list[Path]) -> Path:
    return max(paths, key=lambda p: p.stat().st_size)


def _block_header_offset(blob: bytes) -> tuple[int, int]:
    """(offset of the first block header, its `stored_bytes`).

    Found by validating candidates rather than by a fixed offset: the plaintext file header
    carries a variable-length sample id and a provenance blob in front of the first block. The
    layout is `[Header][BlockHeader][deflate]...[terminator BlockHeader][u64][magic]`, so the last
    block's payload ends 32 bytes before EOF, and that is what identifies the real header among
    byte sequences that merely look like one.
    """
    for h in range(len(blob) - 20):
        n_records, raw_bytes, stored_bytes, _crc, codec = struct.unpack_from("<IIIIB", blob, h)
        if codec == 1 and n_records and raw_bytes and h + 20 + stored_bytes == len(blob) - 32:
            return h, stored_bytes
    raise AssertionError("no .mig block header found -- the on-disk layout changed")


def _damaged_copy(library, tmp_path, name, mutate) -> Path:
    """A whole bucket set with one bucket damaged. Returns the path to point a stage at."""
    dst = tmp_path / name
    shutil.copytree(library["mig_dir"], dst)
    target = dst / _biggest_bucket(library["mig_paths"]).name
    target.write_bytes(mutate(bytearray(target.read_bytes())))
    return dst / sorted(p.name for p in library["mig_paths"])[0]


def test_a_mig_bucket_truncated_mid_block_is_an_error_never_a_short_read(library, tmp_path):
    """A bucket is one range of the barcode space, so stopping early at a torn block silently
    deletes every molecule whose barcode sorts after the cut -- and the remaining ones assemble
    perfectly, which is why nothing else detects it. Delete this and a killed `checkout --mig`
    produces a smaller library that looks entirely healthy."""
    from migec.assemble import run as assemble
    from migec.refine import run as refine

    entry = _damaged_copy(library, tmp_path, "cut", lambda b: bytes(b[:-20]))
    for call in (lambda: assemble(entry, tmp_path / "asm"),
                 lambda: refine(entry, tmp_path / "ref")):
        with pytest.raises(RuntimeError, match="mig_reader"):
            call()


def test_a_flipped_byte_in_a_mig_block_never_reaches_a_consensus(library, tmp_path):
    """Bit rot and a half-written page both look like this. The record layout is column-major and
    fixed-width, so a flipped length turns the whole block's sequences and qualities into a
    consistent-looking shift -- consensuses of bases that were never called. Delete this and there
    is nothing between a damaged intermediate and a plausible FASTQ."""
    from migec.assemble import run as assemble

    def flip_payload(b: bytearray) -> bytes:
        h, stored = _block_header_offset(bytes(b))
        b[h + 20 + stored // 2] ^= 0xFF
        return bytes(b)

    entry = _damaged_copy(library, tmp_path, "flip", flip_payload)
    with pytest.raises(RuntimeError, match="mig_reader"):
        assemble(entry, tmp_path / "asm")


def test_the_mig_block_crc_is_checked_and_not_merely_written(library, tmp_path):
    """zlib's own check covers a payload that fails to inflate; the CRC32 in the block header is
    what covers one that inflates into the wrong bytes, and a checksum nothing compares is a
    checksum that will be wrong silently. Corrupting the RECORDED value is the only way to reach
    that branch without hand-building a deflate stream -- it fails the moment the comparison is
    dropped, which is the regression worth catching."""
    from migec.assemble import run as assemble

    def flip_recorded_crc(b: bytearray) -> bytes:
        h, _ = _block_header_offset(bytes(b))
        b[h + 12] ^= 0xFF  # the crc32 field: n_records, raw_bytes, stored_bytes, then it
        return bytes(b)

    entry = _damaged_copy(library, tmp_path, "crc", flip_recorded_crc)
    with pytest.raises(RuntimeError, match="CRC"):
        assemble(entry, tmp_path / "asm")


# ------------------------------------------------------------------ what a sheet may declare


def test_a_barcode_longer_than_the_packed_key_is_refused_with_both_numbers(tmp_path):
    """A barcode is 2-bit packed into a u64, so 32 bases is the ceiling and a 40 nt UMI would be
    truncated to its first 32 -- silently merging every molecule that shares a prefix. Delete this
    and the failure is a molecule count that is too low with no message anywhere."""
    from migec.checkout import run as checkout

    reads = write_corpus(tmp_path / "reads.fq.gz", n_umis=10)
    sheet = write_sheet(tmp_path / "bc.txt", umi_len=40)
    with pytest.raises(RuntimeError) as exc:
        checkout(reads, sheet, tmp_path / "co")
    message = str(exc.value)
    assert "40" in message and "32" in message, message


@pytest.mark.parametrize("sample_id", ["../evil", "sub/dir", "/tmp/absolute"])
def test_a_sample_id_that_is_a_path_is_refused_before_anything_is_written(tmp_path, sample_id):
    """A sample id becomes a file name, so a sheet is a filesystem write primitive: `../evil`
    put `evil.fq.gz` in the output directory's PARENT and reported every read assigned. Delete
    this and a pasted sheet can overwrite whatever it names, with a summary that says success."""
    from migec.checkout import run as checkout

    reads = write_corpus(tmp_path / "reads.fq.gz", n_umis=10)
    sheet = write_sheet(tmp_path / "bc.txt", sample_id=sample_id)
    out = tmp_path / "co"

    with pytest.raises(RuntimeError) as exc:
        checkout(reads, sheet, out)
    assert sample_id in str(exc.value), str(exc.value)
    # Nothing escaped: no reads file anywhere under the output directory or beside it.
    assert not list(out.glob("*.fq.gz"))
    assert sorted(p.name for p in tmp_path.iterdir()) == ["bc.txt", "co", "reads.fq.gz"]
