#!/bin/bash -l
#SBATCH --job-name=migec
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#
# One sample through all three migec stages.
#
#   sbatch --export=ALL,R1=s1_R1.fq.gz,R2=s1_R2.fq.gz,SAMPLE=s1,PRESET=10x-v2 migec_sample.sbatch
#
# Runs without SLURM too -- every SLURM variable below has a fallback, so
#   R1=s1_R1.fq.gz SAMPLE=s1 PRESET=umi bash migec_sample.sbatch
# does the same work on a laptop. That is how this file is tested.
#
# Never: `-t` changes the wall clock and NOTHING else. Every stage is byte-identical at any thread
# count, which is what makes an escalating `--requeue` or a retry on a different node safe: two
# runs of the same input cannot disagree, so a resumed pipeline never has to explain a diff.

set -euo pipefail

# --- inputs -------------------------------------------------------------------------------------
: "${R1:?set R1=<reads.fq.gz>}"
R2="${R2:-}"
SAMPLE="${SAMPLE:-$(basename "${R1%%.*}")}"
OUT="${OUT:-migec_out/$SAMPLE}"

# Exactly one of these. PRESET is the friendliest; see `migec sheet --presets`.
PRESET="${PRESET:-}"
BC_PATTERN="${BC_PATTERN:-}"
READ_STRUCTURE="${READ_STRUCTURE:-}"

# Which mate still carries sequence after the barcode is trimmed. 2 for 10x, where R1 is barcode
# and nothing else -- trimming the pattern there leaves R1 empty and every payload test vacuous.
PAYLOAD_MATE="${PAYLOAD_MATE:-1}"
RT_ERROR="${RT_ERROR:-rt}"

THREADS="${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN)}"
mkdir -p logs "$OUT"

echo "host=$(hostname)  job=${SLURM_JOB_ID:-local}  sample=$SAMPLE  threads=$THREADS"

layout=()
if   [ -n "$PRESET" ];         then layout=(--preset "$PRESET")
elif [ -n "$BC_PATTERN" ];     then layout=(--bc-pattern "$BC_PATTERN")
elif [ -n "$READ_STRUCTURE" ]; then layout=(--read-structure "$READ_STRUCTURE")
else echo "set one of PRESET, BC_PATTERN or READ_STRUCTURE" >&2; exit 2
fi

# --- checkout: scales with READS ------------------------------------------------------------------
# ~1.5 M reads/s at 16 threads end to end. This is the stage to give cores to.
migec checkout "$R1" ${R2:+"$R2"} "${layout[@]}" \
    --sample "$SAMPLE" -o "$OUT/checkout" -t "$THREADS"

# After checkout the barcode has moved into the read name and been trimmed off the sequence, so
# the payload mate is the one the next two stages read.
if [ -n "$R2" ] && [ "$PAYLOAD_MATE" = "2" ]; then
    TAGGED="$OUT/checkout/${SAMPLE}_R2.fq.gz"
else
    TAGGED="$OUT/checkout/${SAMPLE}.fq.gz"
    [ -f "$TAGGED" ] || TAGGED="$OUT/checkout/${SAMPLE}_R1.fq.gz"
fi

# --- refine: bounded by DISTINCT BARCODES, not reads ----------------------------------------------
# It holds the barcode table in memory and streams the reads three times, so this is the stage that
# wants memory rather than cores. `<sample>.json` reports `table_bytes` -- size the next run from it.
migec refine "$TAGGED" -o "$OUT/refine" -t "$THREADS"

# --- assemble: bounded by ONE BUCKET at a time ----------------------------------------------------
# The range partition means peak memory is set by the bucket count, not the library size.
migec assemble "$OUT/refine/${SAMPLE}.fq.gz" -o "$OUT/assemble" \
    -t "$THREADS" --rt-error "$RT_ERROR"

# --- figures (optional; needs gnuplot) ------------------------------------------------------------
if command -v gnuplot >/dev/null 2>&1; then
    for d in checkout refine assemble; do migec plot "$OUT/$d" -o "$OUT/plots" || true; done
fi

echo "consensus: $OUT/assemble/${SAMPLE}.consensus.fq.gz"
echo done
