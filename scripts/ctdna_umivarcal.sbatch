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
        minimap2 psutil tqdm
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

process() {
    local acc=$1 out=$WORK/$acc
    mkdir -p "$out"
    local fq
    fq=$(fastq_for "$acc") || { echo "  no fastq $acc" >&2; return 1; }

    # 12 nt UMI at the read start, the same layout `migec checkout --bc-pattern '0:12'` reads.
    if [ ! -s "$out/extracted.fastq" ]; then
        "$UVDIR/bin/python" "$SRC/umi-varcal.py" extract -i "$fq" -u 12 -o "$out/extracted.fastq" \
            > "$out/extract.log" 2>&1 || { echo "  extract failed for $acc"; tail -5 "$out/extract.log"; return 1; }
    fi
    if [ ! -s "$out/$acc.bam" ]; then
        minimap2 -ax sr -t "$THREADS" "$MMI" "$out/extracted.fastq" 2>/dev/null \
            | samtools sort -@ "$THREADS" -o "$out/$acc.bam" -
        samtools index "$out/$acc.bam"
    fi
    if [ ! -s "$out/$acc.vcf" ]; then
        # Never: a failed `call` leaves a truncated `.msgpack` pileup behind, and the next attempt
        # reads it and dies with "Unpack failed: incomplete input" -- which reads as a new bug and
        # is the previous one's litter. Keep the expensive artefacts, drop everything else.
        find "$out" -maxdepth 1 -type f ! -name 'extracted.fastq' ! -name "$acc.bam" \
            ! -name "$acc.bam.bai" ! -name '*.log' -delete
        /usr/bin/time -f "umivarcal\t$acc\t%e\t%M" -o "$out/time.txt" \
            "$UVDIR/bin/python" "$SRC/umi-varcal.py" call -i "$out/$acc.bam" -b "$PANEL" \
            -f "$REF" -o "$out" -c "$THREADS" > "$out/call.log" 2>&1 \
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
