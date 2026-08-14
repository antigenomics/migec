#!/bin/bash -l
#SBATCH --job-name=migec-minreads
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# Does --min-reads explain the sensitivity gap?
#
# Published SiMSen-Seq claims <0.1% VAF; TSO500 ctDNA v1 is specified at 0.5% and v2 at 0.2%.
# Scoring our calls against certified truth gave 8% sensitivity at 0.125% and a 25% false-positive
# rate on the 0% arm. That gap says the pipeline is wrong before it says the literature is.
#
# The prime suspect: migec defaults to `--min-reads 1`. A consensus over ONE read IS that read --
# no error correction at all, just counting. UMIErrorCorrect and the SiMSen-Seq pipeline both
# require a family of >= 3 before calling a consensus. Including singletons feeds raw read quality
# straight to the caller, which is exactly where a -> G dark-G artifact would survive.
#
# So: rerun assemble at --min-reads 1, 3 and 5 on a certified subset, and score each.

set -euo pipefail

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
FASTQ=$ROOT/fastq
PANEL=$ROOT/panel/panel_named.bed
MMI=$ROOT/ref/GRCh38.sr.mmi
WORK=$ROOT/minreads
ENVDIR=$ROOT/env
THREADS=${SLURM_CPUS_PER_TASK:-16}
export PATH="$ENVDIR/bin:$PATH"
mkdir -p "$WORK" "$ROOT/logs"

# One replicate per certified arm at the input/depth with the most molecules, so the comparison is
# about the filter and not about running out of evidence. 20 ng / 10x throughout.
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

echo "host=$(hostname) job=${SLURM_JOB_ID:-local}"
echo "arms: ${#ARMS[@]} runs, min-reads 1 / 3 / 5"

for acc in "${!ARMS[@]}"; do
  fq=$(fastq_for "$acc") || { echo "no fastq $acc"; continue; }
  base=$WORK/$acc
  mkdir -p "$base"
  # checkout+refine once: only assemble changes with --min-reads.
  if [ ! -f "$base/rf/$acc.fq.gz" ]; then
      migec checkout "$fq" --bc-pattern '0:12' --sample "$acc" -o "$base/co" -t "$THREADS" \
          > "$base/checkout.log" 2>&1
      migec refine "$base/co/$acc.fq.gz" -o "$base/rf" -t "$THREADS" > "$base/refine.log" 2>&1
      rm -rf "$base/co"
  fi
  for mr in 1 3 5; do
      out=$base/mr$mr
      [ -f "$out/$acc.vcf.gz" ] && continue
      mkdir -p "$out"
      migec assemble "$base/rf/$acc.fq.gz" -o "$out/as" -t "$THREADS" --min-reads "$mr" \
          > "$out/assemble.log" 2>&1
      minimap2 -ax sr -y -t "$THREADS" "$MMI" "$out/as/$acc.consensus.fq.gz" 2>/dev/null \
          | samtools sort -@ "$THREADS" -o "$out/$acc.bam" -
      samtools index "$out/$acc.bam"
      lofreq call -f "$REF" -l "$PANEL" -o "$out/$acc.vcf" "$out/$acc.bam" 2>/dev/null || true
      [ -f "$out/$acc.vcf" ] && { bgzip -f "$out/$acc.vcf"; tabix -f -p vcf "$out/$acc.vcf.gz"; }
      # molecules surviving the filter, at the hotspot's target
      n=$(samtools view -F 0x900 "$out/$acc.bam" "3:179234248-179234329" 2>/dev/null \
          | grep -o 'MI:Z:[^[:space:]]*' | sort -u | wc -l | tr -d ' ')
      echo -e "$acc\t${ARMS[$acc]}\t$mr\t$n" >> "$WORK/molecules.tsv"
  done
  echo "  done $acc (${ARMS[$acc]})"
done

echo "=== results ==="
{ printf "run\tarm\tmin_reads\tchrom\tpos\tref\talt\tdp\taf\n"
  for acc in "${!ARMS[@]}"; do
    for mr in 1 3 5; do
      v=$WORK/$acc/mr$mr/$acc.vcf.gz
      [ -f "$v" ] || continue
      bcftools query -f "${acc}\t${ARMS[$acc]}\t${mr}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n" \
        "$v" 2>/dev/null
    done
  done; } > "$ROOT/minreads_variants.tsv"
echo "wrote $ROOT/minreads_variants.tsv ($(( $(wc -l < "$ROOT/minreads_variants.tsv") - 1 )) calls)"
sort -u "$WORK/molecules.tsv" > "$ROOT/minreads_molecules.tsv"
echo done
