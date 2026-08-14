#!/bin/bash -l
#SBATCH --job-name=migec-uec
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
# Never: `medium` tops out at 16 hours and a job asking for more sits at PENDING
# (PartitionTimeLimit) for ever rather than being rejected. The bwa index is ~1 h of this.
#SBATCH --time=15:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# UMIErrorCorrect on the same certified arms: PIPELINE against pipeline, not caller against caller.
#
# Never: UMIErrorCorrect does NOT follow `assemble`, it replaces it. It aligns raw reads first and
# then groups on (chromosomal position, UMI) at edit distance <= 1, which is the opposite order to
# migec; running it on a migec consensus would collapse an already-collapsed library. So the honest
# comparison is both pipelines from the same FASTQ to a call set, and the rows say "pipeline"
# rather than "caller" (docs/variants.rst).
#
# It is also the RIGHT comparator for this dataset in particular: Osterlund et al. wrote both the
# tool and `PRJNA788522`, so its defaults were chosen on exactly these runs.
#
# Never: MATCH THE GROUP SIZE. UMIErrorCorrect's consensus cutoff is >= 3 reads per UMI by default
# and migec's `--min-reads` is 1. `ctdna_callers.sbatch` already sweeps 1/3/5 on our side, so 3 is
# the row that compares like with like, and 1 is reported next to it rather than instead of it.

set -euo pipefail

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
FASTQ=$ROOT/fastq
PANEL=$ROOT/panel/panel_named.bed
WORK=$ROOT/uec
UECDIR=$ROOT/env-uec
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$WORK" "$ROOT/logs" "$ROOT/ref"

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
if [ ! -x "$UECDIR/bin/bwa" ]; then
    echo "=== creating the UMIErrorCorrect environment ==="
    "$MICROMAMBA" create -y -p "$UECDIR" -c conda-forge -c bioconda \
        python=3.10 umierrorcorrect bwa samtools
fi
export PATH="$UECDIR/bin:$PATH"
# Logged before anything is run: if a flag has moved between versions, this is the line that says
# so, rather than twelve identical failures.
run_umierrorcorrect.py --help > "$WORK/usage.txt" 2>&1 || true
head -40 "$WORK/usage.txt"

# --- bwa index, once ---------------------------------------------------------------------------
# UMIErrorCorrect maps with `bwa mem`, and there is no plain bwa index of GRCh38 on this cluster --
# only bwameth and bwa-mem2. Building it here keeps the comparison on GRCh38 coordinates, which is
# what every other table in this directory is on; extracting the panel into a small reference would
# have needed a coordinate conversion nobody could check.
LOCALREF=$ROOT/ref/GRCh38.fa
[ -e "$LOCALREF" ] || ln -s "$REF" "$LOCALREF"
if [ ! -s "$LOCALREF.bwt" ]; then
    echo "=== bwa index (once, ~1 h) ==="
    time bwa index "$LOCALREF"
fi
ls -la "$LOCALREF".* | awk '{print $9, $5}'

# --- one run -------------------------------------------------------------------------------------
process() {
    local acc=$1 out=$WORK/$acc
    [ -s "$out/$acc.vcf" ] && { echo "  have $acc"; return 0; }
    local fq
    fq=$(fastq_for "$acc") || { echo "  no fastq $acc" >&2; return 1; }
    mkdir -p "$out"
    # `-ul 12 -sl 16` is SiMSen-Seq's layout as the tool's own paper states it, and it is the same
    # 12 nt UMI `migec checkout --bc-pattern '0:12'` reads.
    /usr/bin/time -f "umierrorcorrect\t$acc\t%e\t%M" -o "$out/time.txt" \
        run_umierrorcorrect.py -r1 "$fq" -r "$LOCALREF" -o "$out" \
        -ul 12 -sl 16 -t "$THREADS" -bed "$PANEL" > "$out/run.log" 2>&1 || true
    ls "$out" > "$out/listing.txt"
    echo "  done $acc (${ARMS[$acc]})"
}

echo "=== UMIErrorCorrect: ${#ARMS[@]} runs ==="
for acc in "${!ARMS[@]}"; do
    process "$acc" || true
done

# --- collect whatever it wrote ---------------------------------------------------------------------
# Its output naming is version-dependent, so the glob is wide and the file that matched is named in
# the table rather than assumed.
{
  printf "run\tarm\tsource\tchrom\tpos\tref\talt\tdp\taf\n"
  for acc in "${!ARMS[@]}"; do
    for v in "$WORK/$acc"/*.vcf "$WORK/$acc"/*/*.vcf; do
      [ -s "$v" ] || continue
      awk -v acc="$acc" -v arm="${ARMS[$acc]}" -v src="$(basename "$v")" \
          'BEGIN{FS=OFS="\t"} !/^#/ {
             dp=""; af="";
             n=split($8, a, ";");
             for (i=1;i<=n;i++) { if (a[i] ~ /^DP=/) sub("DP=","",a[i]), dp=a[i];
                                  if (a[i] ~ /^AF=/) sub("AF=","",a[i]), af=a[i] }
             print acc, arm, src, $1, $2, $4, $5, dp, af }' "$v"
    done
  done
} > "$ROOT/uec_variants.tsv"
echo "wrote $ROOT/uec_variants.tsv ($(( $(wc -l < "$ROOT/uec_variants.tsv") - 1 )) calls)"
cat "$WORK"/*/time.txt > "$ROOT/uec_runtime.tsv" 2>/dev/null || true
echo done
