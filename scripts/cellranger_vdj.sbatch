#!/bin/bash
#SBATCH --job-name=cellranger-vdj
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=%x-%j.out
# 2026-08-14
#
# `cellranger vdj` on sc5p_v2_hs_PBMC_1k VDJ-T, so the migec comparison has a Cell Ranger side that
# was actually RUN rather than only read off 10x's published tables. Two things come out of it that
# the published output cannot give:
#   1. wall clock and peak RSS, which is the cost axis every other compare_* script reports;
#   2. a SECOND Cell Ranger version. The published calls are 5.0.0; this is 10.1.0, so the gap
#      between them bounds how much of any migec-vs-Cell-Ranger difference is just version drift.
#
# Never: `/projects/cdr3*` is read-only for us. This lives in /projects/immunestatus, which had
# 4.6T free at 52% used -- checked before anything was downloaded.
#
# Never: compute never runs on an aldan3 frontend. The download was I/O and ran there deliberately;
# this is 16 cores for hours and belongs in SLURM.
#
# Submit:
#   aldan3 push scripts/cellranger_vdj.sbatch /projects/immunestatus/migec-cellranger/
#   aldan3 slurm submit /projects/immunestatus/migec-cellranger/cellranger_vdj.sbatch
set -euo pipefail

D=/projects/immunestatus/migec-cellranger
CR="$D/cellranger-10.1.0/cellranger"
REF="$D/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0"
FASTQS="$D/sc5p_v2_hs_PBMC_1k_t_fastqs"
RUN="$D/run"

mkdir -p "$RUN" && cd "$RUN"
# Cell Ranger refuses to overwrite an existing run id, and a half-finished one is worse than none.
rm -rf PBMC_vdj

echo "=== cellranger $($CR --version) on $(hostname), $SLURM_CPUS_PER_TASK cores ==="
echo "=== started $(date) ==="

/usr/bin/time -v "$CR" vdj \
    --id=PBMC_vdj \
    --reference="$REF" \
    --fastqs="$FASTQS" \
    --sample=sc5p_v2_hs_PBMC_1k_t \
    --localcores="$SLURM_CPUS_PER_TASK" \
    --localmem=60 \
    --disable-ui

echo "=== done $(date) ==="
echo "=== metrics ==="
cat PBMC_vdj/outs/metrics_summary.csv
echo "=== cells ==="
tail -n +2 PBMC_vdj/outs/filtered_contig_annotations.csv | cut -d, -f1 | sort -u | wc -l
