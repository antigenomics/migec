#!/bin/bash -l
#SBATCH --job-name=migec-umivarcal
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# UMI-VarCal on the same certified arms: the third pipeline, and the one whose input convention is
# "UMI in the read NAME".
#
# Never: like UMIErrorCorrect, it REPLACES assemble rather than following it -- it builds its own
# UMI-aware pileup from an alignment of raw reads and calls from that. So it runs from the raw
# FASTQ with its own `extract`, and its row says "pipeline" rather than "caller".
#
# Its extraction is also the convention migec deliberately does not read: `extract` moves the first
# `-u` bases of each read into the read name. migec never parses read names (docs/byo_umi.rst) --
# a name is free-form and a tag is not -- so this is the comparison that shows what that costs and
# what it buys, from the other side.

set -euo pipefail

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
FASTQ=$ROOT/fastq
PANEL=$ROOT/panel/panel_named.bed
MMI=$ROOT/ref/GRCh38.sr.mmi
# Never: give UMI-VarCal the SYMLINK in our own ref/, not the shared reference path. `ParseFASTA`
# caches a msgpack index of the FASTA beside it as `<ref>.dic` -- so pointing it at
# /projects/cdr3_common writes a multi-GB index into a shared, read-mostly reference tree, and a
# killed run leaves a 0-byte `.dic` there that every later run of every user then "finds".
LOCALREF=$ROOT/ref/GRCh38.fa
WORK=$ROOT/umivarcal
UVDIR=$ROOT/env-uv
SRC=$ROOT/src-umivarcal
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$WORK" "$ROOT/logs"

declare -A ARMS=(
  [SRR17220921]=WT [SRR17220920]=WT [SRR17220919]=WT
  [SRR17220931]=0.125 [SRR17220930]=0.125 [SRR17220929]=0.125
  [SRR17220928]=0.25 [SRR17220926]=0.25 [SRR17220925]=0.25
  [SRR17220924]=1 [SRR17220923]=1 [SRR17220922]=1
)

fastq_for() {
  for c in "$FASTQ/$1_1.fastq.gz" "$FASTQ/$1.fastq.gz"; do [ -f "$c" ] && { echo "$c"; return 0; }; done
  return 1
}

echo "host=$(hostname) job=${SLURM_JOB_ID:-local} threads=$THREADS"

# Never: `if not os.path.isfile(dictPath)` treats a 0-BYTE index as a valid one, so one killed run
# poisons every run after it with "Unpack failed: incomplete input" from a line that says
# "Index found". Drop an empty one, here and in the shared tree an earlier attempt polluted.
for d in "${LOCALREF%.fa}.dic" \
         /projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.dic; do
    [ -f "$d" ] && [ ! -s "$d" ] && { rm -f "$d" && echo "removed empty index $d"; }
done || true

export MAMBA_ROOT_PREFIX=$ROOT/mamba
MICROMAMBA=$HOME/bin/micromamba
# Never: test for `conda-meta/history`, not for `bin/python`. A half-deleted prefix still has a
# python and micromamba no longer calls it an environment, so a `bin/python` guard skips the
# rebuild and every later `install -p` fails with "Expected environment not found at prefix".
# Never: the check has to require msgpack >= 1.0, not merely that it imports. The patch below
# DELETES `encoding="utf-8"`, and that is only equivalent on 1.x, where `raw=False` is the default;
# on 0.6.2 `raw` still defaults to True and the tool would silently get BYTES keys where it expects
# strings. An env left over from an earlier attempt satisfies a bare `import msgpack` and turns the
# patch from a no-op into a behaviour change.
usable_env() {
    [ -s "$UVDIR/conda-meta/history" ] || return 1
    [ -x "$UVDIR/bin/bwa" ] || return 1
    "$UVDIR/bin/python" -c "import psutil, tqdm, msgpack, sys; sys.exit(0 if msgpack.version >= (1, 0) else 1)" 2>/dev/null
}

if ! usable_env; then
    echo "=== creating the UMI-VarCal environment ==="
    rm -rf "$UVDIR"
    # Never: `psutil` is NOT in its README's requirements list and the tool does not start without
    # it -- the import is inside the compiled `functions.pyx`, so it surfaces as a
    # ModuleNotFoundError from a .pyx line and not as a missing dependency. Same for `tqdm`.
    "$MICROMAMBA" create -y -p "$UVDIR" -c conda-forge -c bioconda \
        python=3.10 cython msgpack-python scipy statsmodels pysam numpy pip samtools \
        minimap2 bwa psutil tqdm
    "$UVDIR/bin/pip" install --quiet pyfasta || echo "pyfasta from pip failed"
fi
export PATH="$UVDIR/bin:$PATH"
"$UVDIR/bin/python" -c "import msgpack; print('msgpack', msgpack.version)"

# It is a git checkout with a Cython extension to build, not a package.
if [ ! -d "$SRC/.git" ]; then
    git clone --recursive https://gitlab.com/vincent-sater/umi-varcal.git "$SRC"
fi
# Never: THE ONE PATCH, and it is a serialization CALL, not the algorithm.
#
# UMI-VarCal calls `msgpack.pack(..., encoding="utf-8")` and
# `msgpack.unpack(..., encoding="utf-8", max_buffer_size=<huge>)` in eleven places. `encoding` was
# removed in msgpack 1.0 (2019) and `max_buffer_size` is not an `unpackb` argument in 0.6 either,
# so the tool runs on NO msgpack a solver will give you: 1.x rejects `encoding`, 0.6.2 rejects
# `max_buffer_size`, and the last release accepting both is 0.5.6, which conda-forge builds only
# for Python <= 3.6.
#
# Both kwargs are simply DELETED, which is msgpack's own documented migration: `raw=False` (what
# `encoding="utf-8"` did) is the 1.x default for unpack, `use_bin_type=True` is the 1.x default for
# pack, and `max_buffer_size` was set to 10^12 to disable a limit that does not apply to a bytes
# input. Nothing about the pileup, the grouping or the calling changes -- but it IS a patched
# comparator and every row it produces has to say so.
if [ -d "$SRC/functions" ] && [ ! -f "$SRC/.migec-msgpack-patched" ]; then
    echo "=== patching the msgpack calls (see the note above) ==="
    grep -rn 'encoding="utf-8"' "$SRC/functions"/*.pyx "$SRC/functions"/*.py \
        > "$WORK/msgpack_calls_before.txt" 2>/dev/null || true
    sed -i -E 's/, *encoding="utf-8" *, *max_buffer_size=[^)]*\)/)/g; s/, *encoding="utf-8" *\)/)/g' \
        "$SRC/functions"/*.pyx "$SRC/functions"/*.py
    touch "$SRC/.migec-msgpack-patched"
    echo "patched $(wc -l < "$WORK/msgpack_calls_before.txt") call sites"
fi
if [ -d "$SRC/functions" ]; then
    (cd "$SRC/functions" && "$UVDIR/bin/python" setup.py build_ext --inplace) \
        > "$WORK/build.log" 2>&1 || { echo "BUILD FAILED -- see $WORK/build.log"; tail -20 "$WORK/build.log"; }
fi
ls "$SRC" | head

# Never: `extract` takes a BAM OR SAM, NOT a FASTQ; its flags are `-i`, `-f` and
# `--umi_length | -l` -- there is no `-u` and no `-o` -- and it writes into `EXTRACTED/<name>/`
# relative to the CURRENT DIRECTORY. Reading `functions/Config.py` and `functions/Extract.py` once
# would have replaced five failed guesses at an interface that is written down.
#
# Never: `--umi_length 0`, which the SOURCE says reads the UMI from the RX tag instead of from the
# read start, CANNOT BE PASSED: the validator tests the value for truthiness, so 0 is
# indistinguishable from absent and it exits with "The parameter --umi_length | -l is required but
# missing". So the RX path is unreachable from the CLI and the tool gets the reads with the UMI
# still at the 5' end, which is what it was written for.
#
# The alignment in the input BAM is not used -- `extract` reads sequences and reverts to FASTQ,
# then realigns with bwa itself -- so an UNALIGNED bam built by `samtools import` is enough, and it
# costs seconds instead of an alignment we would throw away.
process() {
    local acc=$1 out=$WORK/$acc
    local fq
    fq=$(fastq_for "$acc") || { echo "  no fastq $acc" >&2; return 1; }
    mkdir -p "$out"
    local ubam="$out/$acc.raw.bam"
    # Never: the final BAM is `<name>_extracted.bam`, not `<name>_extracted_sorted.bam`. The
    # source sorts INTO the `_sorted` name and then renames it back over the unsorted one
    # (`Extract.py` 187-191), so the `_sorted` file exists only between two statements.
    local extracted="$out/EXTRACTED/$acc.raw/$acc.raw_extracted.bam"

    [ -s "$ubam" ] || samtools import -0 "$fq" -o "$ubam"

    if [ ! -s "$extracted" ]; then
        # 12 nt UMI at the read start -- the same layout `migec checkout --bc-pattern '0:12'`
        # reads. The 16 nt SiMSen-Seq spacer behind it is left in place; bwa soft-clips it, which
        # is what this tool does with it. UMIErrorCorrect has `-sl` for the spacer and this one has
        # no concept of it.
        ( cd "$out" && "$UVDIR/bin/python" "$SRC/umi-varcal.py" extract \
            -i "$ubam" -f "$LOCALREF" -l 12 ) > "$out/extract.log" 2>&1 \
            || { echo "  extract failed for $acc"; tail -5 "$out/extract.log"; return 1; }
    fi
    # Never: EMPTY OUTPUT HERE MEANS THE LIBRARY IS SINGLE-END, not that something misfired.
    # `Extract.py` pairs reads BY READ ID -- it holds each read in a `seen` dict and writes to R1/R2
    # only when the same id turns up a SECOND time (lines 138-154). On a single-end library every
    # id appears once, so both FASTQs stay empty, `bwa mem` aligns nothing, and the extracted BAM
    # has 0 records after reporting "Working 100 %". The tool is for paired-end libraries, as its
    # own description says, and `PRJNA788522` is single-end 151 nt.
    if [ ! -s "$extracted" ] || [ "$(samtools view -c "$extracted" 2>/dev/null || echo 0)" -eq 0 ]; then
        echo "  extract produced an EMPTY bam for $acc -- UMI-VarCal pairs reads by id and this"
        echo "  library is single-end, so there is nothing it can do with it. Not a config error."
        return 1
    fi

    if ! compgen -G "$out"/*.vcf > /dev/null; then
        /usr/bin/time -f "umivarcal\t$acc\t%e\t%M" -o "$out/time.txt" \
            "$UVDIR/bin/python" "$SRC/umi-varcal.py" call -i "$extracted" -b "$PANEL" \
            -f "$LOCALREF" -o "$out" -c "$THREADS" > "$out/call.log" 2>&1 \
            || { echo "  call failed for $acc"; tail -5 "$out/call.log"; return 1; }
    fi
    echo "  done $acc (${ARMS[$acc]})"
}

echo "=== UMI-VarCal: ${#ARMS[@]} runs ==="
for acc in "${!ARMS[@]}"; do
    process "$acc" || true
done

{
  printf "run\tarm\tfiltered\tchrom\tpos\tref\talt\tdp\taf\n"
  for acc in "${!ARMS[@]}"; do
    for v in "$WORK/$acc"/*.vcf; do
      [ -s "$v" ] || continue
      "$UVDIR/bin/python" - "$v" "$acc" "${ARMS[$acc]}" <<'PY'
import sys
path, acc, arm = sys.argv[1:4]
for line in open(path):
    if line.startswith("#") or not line.strip():
        continue
    f = line.rstrip("\n").split("\t")
    if len(f) < 8:
        continue
    info = dict(kv.split("=", 1) for kv in f[7].split(";") if "=" in kv)
    print("\t".join([acc, arm, f[6], f[0], f[1], f[3], f[4],
                     info.get("DP", ""), info.get("AF", info.get("VAF", ""))]))
PY
    done
  done
} > "$ROOT/umivarcal_variants.tsv"
echo "wrote $ROOT/umivarcal_variants.tsv ($(( $(wc -l < "$ROOT/umivarcal_variants.tsv") - 1 )) calls)"
echo done
