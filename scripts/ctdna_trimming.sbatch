#!/bin/bash -l
#SBATCH --job-name=migec-trimming
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=10:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# Does adapter trimming actually remove the off-target load, and does it change the calls?
#
# `docs/detection.rst` diagnoses the off-target share -- 5-7% at 80 ng rising to 47-58% at 5 ng --
# as TruSeq adapter read-through: 85-95 bases soft-clipped at MAPQ 4-16, carrying
# GATCGGAAGAGCACACGTCTGAACTCCAGTCAC, mismapping into a 97%-G locus. The diagnosis has been in the
# docs since the run that found it and has never been tested by applying the fix.
#
# `ctdna_callers.sbatch` has since measured the OTHER half of that table and found it wrong: a
# MAPQ 20 floor removes the mismapped alignments and RAISES the false-positive burden 2.6x at
# --min-reads 1, because it moves the caller's threshold rather than removing the artifact. So the
# remaining claim -- that trimming removes the cause -- has to be measured, not inherited.
#
# One thing changes: cutadapt before checkout. Same barcode pattern, same stages, same panel, same
# caller, same thresholds.

set -euo pipefail

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
FASTQ=$ROOT/fastq
PANEL=$ROOT/panel/panel_named.bed
MMI=$ROOT/ref/GRCh38.sr.mmi
WORK=$ROOT/trimming
ENVDIR=$ROOT/env
TRIMDIR=$ROOT/env-trim
THREADS=${SLURM_CPUS_PER_TASK:-16}
# Illumina TruSeq. The full 33-mer is what was found in the soft-clipped tails; cutadapt matches a
# prefix of it, so a read that ran only a few bases into the adapter is trimmed too.
ADAPTER=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

export PATH="$ENVDIR/bin:$PATH"
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
if [ ! -x "$TRIMDIR/bin/cutadapt" ]; then
    echo "=== creating the trimming environment ==="
    "$MICROMAMBA" create -y -p "$TRIMDIR" -c conda-forge -c bioconda cutadapt
fi
"$TRIMDIR/bin/cutadapt" --version

# Molecules inside the panel against molecules anywhere, by MI tag -- one consensus record is one
# molecule, so this is a count and not a coverage average.
off_target() {
    local bam=$1
    local total on
    total=$(samtools view -F 0x900 "$bam" | grep -o 'MI:Z:[^[:space:]]*' | sort -u | wc -l)
    on=$(samtools view -F 0x900 -L "$PANEL" "$bam" | grep -o 'MI:Z:[^[:space:]]*' | sort -u | wc -l)
    echo "$total $on"
}

process() {
    local acc=$1 mr=$2 out=$WORK/$acc
    mkdir -p "$out"
    local fq
    fq=$(fastq_for "$acc") || { echo "  no fastq $acc" >&2; return 1; }

    # Never: the UMI is at the 5' end and the adapter read-through is at the 3' end, so `-a` only.
    # A `-g` here would eat the barcode `checkout` is about to read.
    if [ ! -s "$out/trimmed.fq.gz" ]; then
        "$TRIMDIR/bin/cutadapt" -a "$ADAPTER" -m 40 -j "$THREADS" \
            -o "$out/trimmed.fq.gz" "$fq" > "$out/cutadapt.log" 2>&1
    fi
    if [ ! -s "$out/rf/$acc.fq.gz" ]; then
        migec checkout "$out/trimmed.fq.gz" --bc-pattern '0:12' --sample "$acc" \
              -o "$out/co" -t "$THREADS" > "$out/checkout.log" 2>&1
        migec refine "$out/co/$acc.fq.gz" -o "$out/rf" -t "$THREADS" > "$out/refine.log" 2>&1
        rm -rf "$out/co"
    fi
    local d=$out/mr$mr
    if [ ! -s "$d/$acc.bam" ]; then
        mkdir -p "$d"
        migec assemble "$out/rf/$acc.fq.gz" -o "$d/as" -t "$THREADS" --min-reads "$mr" \
              > "$d/assemble.log" 2>&1
        minimap2 -ax sr -y -t "$THREADS" "$MMI" "$d/as/$acc.consensus.fq.gz" 2>/dev/null \
            | samtools sort -@ "$THREADS" -o "$d/$acc.bam" -
        samtools index "$d/$acc.bam"
    fi
    [ -s "$d/$acc.vcf.gz" ] || {
        lofreq call -f "$REF" -l "$PANEL" -o "$d/$acc.vcf" "$d/$acc.bam" 2>/dev/null || true
        [ -f "$d/$acc.vcf" ] && bgzip -f "$d/$acc.vcf"
    }
    echo "  done $acc mr$mr (${ARMS[$acc]})"
}

echo "=== trimmed pipeline: ${#ARMS[@]} runs x min-reads 1/3 ==="
for acc in "${!ARMS[@]}"; do
    for mr in 1 3; do
        process "$acc" "$mr" || true
    done
done

# --- the two tables ---------------------------------------------------------------------------------
echo "=== off-target share, trimmed against untrimmed, on the SAME runs ==="
{
  printf "run\tarm\tmin_reads\ttrimmed\tmolecules_total\tmolecules_on_target\toff_target_share\n"
  for acc in "${!ARMS[@]}"; do
    for mr in 1 3; do
      for kind in untrimmed trimmed; do
        if [ "$kind" = untrimmed ]; then bam=$ROOT/minreads/$acc/mr$mr/$acc.bam
        else bam=$WORK/$acc/mr$mr/$acc.bam; fi
        [ -s "$bam" ] || continue
        read -r total on < <(off_target "$bam")
        [ "$total" -gt 0 ] || continue
        share=$(python3 -c "print(f'{1 - $on / $total:.4f}')")
        flag=0; [ "$kind" = trimmed ] && flag=1
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$acc" "${ARMS[$acc]}" "$mr" "$flag" "$total" "$on" "$share"
      done
    done
  done
} > "$ROOT/trimming_offtarget.tsv"
echo "wrote $ROOT/trimming_offtarget.tsv"

{
  printf "run\tarm\tmin_reads\ttrimmed\tchrom\tpos\tref\talt\tdp\taf\n"
  for acc in "${!ARMS[@]}"; do
    for mr in 1 3; do
      v=$WORK/$acc/mr$mr/$acc.vcf.gz
      [ -s "$v" ] && bcftools query -f \
        "${acc}\t${ARMS[$acc]}\t${mr}\t1\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n" \
        "$v" 2>/dev/null || true
    done
  done
} > "$ROOT/trimming_variants.tsv"
echo "wrote $ROOT/trimming_variants.tsv ($(( $(wc -l < "$ROOT/trimming_variants.tsv") - 1 )) calls)"
echo done
