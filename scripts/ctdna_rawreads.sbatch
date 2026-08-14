#!/bin/bash -l
#SBATCH --job-name=migec-rawreads
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# The baseline every other row is measured against: the same reads, the same aligner, the same
# callers -- and no consensus.
#
# Without this row the caller comparison says which caller is better on migec output and nothing
# about whether collapsing was worth doing. With it, one thing differs between the two rows and it
# is the only thing migec does.
#
# Never: it is NOT "the raw FASTQ". The reads here are checkout's and refine's output -- barcode
# found, trimmed and corrected -- so the alignment and the trimming are identical to the consensus
# arm and the ONLY difference is that a read is a read rather than a molecule. Aligning the
# untrimmed FASTQ instead would soft-clip 28 bases of barcode and spacer off every read and the
# comparison would be measuring that.

set -euo pipefail

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
PANEL=$ROOT/panel/panel_named.bed
MMI=$ROOT/ref/GRCh38.sr.mmi
WORK=$ROOT/rawreads
ENVDIR=$ROOT/env
GATKDIR=$ROOT/env-gatk
LOCALREF=$ROOT/callers/ref/GRCh38.fa
THREADS=${SLURM_CPUS_PER_TASK:-16}
export PATH="$ENVDIR/bin:$PATH"
GATK=$GATKDIR/bin/gatk
mkdir -p "$WORK" "$ROOT/logs"

declare -A ARMS=(
  [SRR17220921]=WT [SRR17220920]=WT [SRR17220919]=WT
  [SRR17220931]=0.125 [SRR17220930]=0.125 [SRR17220929]=0.125
  [SRR17220928]=0.25 [SRR17220926]=0.25 [SRR17220925]=0.25
  [SRR17220924]=1 [SRR17220923]=1 [SRR17220922]=1
)

echo "host=$(hostname) job=${SLURM_JOB_ID:-local} threads=$THREADS"

for acc in "${!ARMS[@]}"; do
    src=$ROOT/minreads/$acc/rf/$acc.fq.gz
    [ -s "$src" ] || { echo "  no refined reads for $acc" >&2; continue; }
    out=$WORK/$acc
    mkdir -p "$out"
    if [ ! -s "$out/$acc.bam" ]; then
        minimap2 -ax sr -y -t "$THREADS" "$MMI" "$src" 2>/dev/null \
            | samtools sort -@ "$THREADS" -o "$out/reads.bam" -
        samtools addreplacerg -r "$(printf '@RG\tID:1\tSM:%s\tPL:ILLUMINA\tLB:%s' "$acc" "$acc")" \
            -o "$out/$acc.bam" "$out/reads.bam"
        samtools index "$out/$acc.bam"
        rm -f "$out/reads.bam"
    fi
    if [ ! -s "$out/lofreq.vcf.gz" ]; then
        lofreq call -f "$LOCALREF" -l "$PANEL" -o "$out/lofreq.vcf" "$out/$acc.bam" 2>/dev/null || true
        [ -f "$out/lofreq.vcf" ] && bgzip -f "$out/lofreq.vcf"
    fi
    if [ ! -s "$out/mutect2.vcf.gz" ]; then
        # Same `--max-reads-per-alignment-start 0` as the consensus arm. Here it matters even more:
        # an amplicon's reads all start at the same coordinate, so the default 50 would leave a
        # thousandth of them.
        "$GATK" Mutect2 -R "$LOCALREF" -L "$PANEL" -I "$out/$acc.bam" \
            --max-reads-per-alignment-start 0 -O "$out/mutect2.raw.vcf.gz" \
            > "$out/mutect2.log" 2>&1 || true
        [ -s "$out/mutect2.raw.vcf.gz" ] && "$GATK" FilterMutectCalls -R "$LOCALREF" \
            -V "$out/mutect2.raw.vcf.gz" -O "$out/mutect2.vcf.gz" >> "$out/mutect2.log" 2>&1 || true
    fi
    echo "  done $acc (${ARMS[$acc]})"
done

{
  printf "run\tarm\tmin_reads\tcaller\tmapq_floor\tfiltered\tchrom\tpos\tref\talt\tdp\taf\n"
  for acc in "${!ARMS[@]}"; do
    v=$WORK/$acc/lofreq.vcf.gz
    [ -s "$v" ] && bcftools query -f \
      "${acc}\t${ARMS[$acc]}\t0\tlofreq\t0\tPASS\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n" \
      "$v" 2>/dev/null || true
    v=$WORK/$acc/mutect2.vcf.gz
    [ -s "$v" ] && bcftools query -f \
      "${acc}\t${ARMS[$acc]}\t0\tmutect2\t0\t%FILTER\t%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AF]\n" \
      "$v" 2>/dev/null || true
  done
} > "$ROOT/rawreads_variants.tsv"
echo "wrote $ROOT/rawreads_variants.tsv ($(( $(wc -l < "$ROOT/rawreads_variants.tsv") - 1 )) calls)"
echo done
