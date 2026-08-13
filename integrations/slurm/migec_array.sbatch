#!/bin/bash -l
#SBATCH --job-name=migec-array
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#
# One array task per sample. This is the shape a real cohort wants: samples are independent, they
# finish at different times, and a failure re-runs one task rather than the whole cohort.
#
#   sbatch --array=1-$(($(wc -l < samples.tsv) - 1)) migec_array.sbatch samples.tsv
#
# samples.tsv is a header plus one row per sample:
#
#   sample  r1                r2                layout      payload_mate
#   PBMC1   PBMC1_R1.fq.gz    PBMC1_R2.fq.gz    10x-v2      2
#   HIV     SRR1763769_2.fq.gz  -               primerid    1
#
# `layout` is a preset name, a pattern (`^NNNNNNNN`), or a slice list (`0:12`); the script tells
# them apart. Use `-` for an absent R2.
#
# Runs without SLURM: with no SLURM_ARRAY_TASK_ID it does row 1, so
#   bash migec_array.sbatch samples.tsv
# is a single-sample smoke test of the same code path. That is how this file is tested.

set -euo pipefail

SHEET="${1:-samples.tsv}"
[ -f "$SHEET" ] || { echo "no sample sheet at $SHEET" >&2; exit 2; }

# Row 1 is the first DATA row -- the header is skipped, not counted. Getting this off by one is
# the classic array-job bug: it silently processes the header as a sample and the run "succeeds".
TASK="${SLURM_ARRAY_TASK_ID:-1}"
ROW=$(tail -n +2 "$SHEET" | sed -n "${TASK}p")
[ -n "$ROW" ] || { echo "no row $TASK in $SHEET (it has $(($(wc -l < "$SHEET") - 1)) samples)" >&2; exit 2; }

IFS=$'\t ' read -r SAMPLE R1 R2 LAYOUT PAYLOAD_MATE <<<"$ROW"
[ "${R2:-}" = "-" ] && R2=""
PAYLOAD_MATE="${PAYLOAD_MATE:-1}"
OUT="${OUT:-migec_out}/$SAMPLE"
THREADS="${SLURM_CPUS_PER_TASK:-$(getconf _NPROCESSORS_ONLN)}"
RT_ERROR="${RT_ERROR:-rt}"

mkdir -p logs "$OUT"
echo "host=$(hostname) task=$TASK sample=$SAMPLE threads=$THREADS layout=$LAYOUT"

# A preset is a bare name; a pattern has N/X/^ or a colon. Deciding here rather than asking the
# sheet to say which keeps the sheet readable.
if [[ "$LAYOUT" == *:* || "$LAYOUT" == *N* || "$LAYOUT" == ^* ]]; then
    layout=(--bc-pattern "$LAYOUT")
else
    layout=(--preset "$LAYOUT")
fi

migec checkout "$R1" ${R2:+"$R2"} "${layout[@]}" --sample "$SAMPLE" -o "$OUT/checkout" -t "$THREADS"

if [ -n "$R2" ] && [ "$PAYLOAD_MATE" = "2" ]; then
    TAGGED="$OUT/checkout/${SAMPLE}_R2.fq.gz"
else
    TAGGED="$OUT/checkout/${SAMPLE}.fq.gz"
    [ -f "$TAGGED" ] || TAGGED="$OUT/checkout/${SAMPLE}_R1.fq.gz"
fi

migec refine "$TAGGED" -o "$OUT/refine" -t "$THREADS"
migec assemble "$OUT/refine/${SAMPLE}.fq.gz" -o "$OUT/assemble" -t "$THREADS" --rt-error "$RT_ERROR"

echo "consensus: $OUT/assemble/${SAMPLE}.consensus.fq.gz"
echo done
