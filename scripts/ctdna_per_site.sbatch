#!/bin/bash -l
#SBATCH --job-name=migec-persite
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# Phase 2: per-TARGET molecule counts and variant calls, using the panel phase 1 inferred.
#
# This replaces `molecules / n_amplicons` with a count. The old figure assumed even coverage and
# the panel is not evenly covered (20.4 / 17.4 / 16.9 / 15.9 / 7.6%), so it overstated the weakest
# target by 2.6x -- the difference between P(>=3 copies) = 0.66 and 0.14 at 5 ng and 0.125% VAF.
#
# A molecule is counted by its MI tag, not by its reads: one consensus record is one molecule, so
# `samtools view <region> | cut MI | sort -u | wc -l` is the molecule count at that target.

set -euo pipefail

# Never: THE FASTQ NAME IS NOT ONE PATTERN. `sra_fetch.py get` writes `<acc>_1.fastq.gz` when it
# converted an `.sra` with fasterq-dump, but `<acc>.fastq.gz` when the archive had no `.sra` and it
# took the ORIGINAL SUBMITTED file instead. 13 of these 100 runs are the second kind, and a fixed
# `${acc}_1.fastq.gz` skips them silently -- which is how the first pass lost 2 of its 5 panel runs
# and still reported success.
fastq_for() {
    local acc=$1
    for c in "$FASTQ/${acc}_1.fastq.gz" "$FASTQ/${acc}.fastq.gz"; do
        [ -f "$c" ] && { echo "$c"; return 0; }
    done
    return 1
}

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
FASTQ=$ROOT/fastq
PANEL=${PANEL:-$ROOT/panel/panel_named.bed}
WORK=$ROOT/persite
ENVDIR=$ROOT/env
THREADS=${SLURM_CPUS_PER_TASK:-16}
LIMIT=${LIMIT:-0}          # 0 = every run; set e.g. 12 for a quick pass

export PATH="$ENVDIR/bin:$PATH"
mkdir -p "$WORK" "$ROOT/logs"

[ -s "$PANEL" ] || { echo "no panel at $PANEL -- run infer_panel.sbatch first" >&2; exit 2; }
echo "host=$(hostname) job=${SLURM_JOB_ID:-local} threads=$THREADS"
echo "panel: $(wc -l < "$PANEL") targets"

# Never: BUILD THE INDEX ONCE. `minimap2 ... ref.fa` re-indexes the 3 GB genome on EVERY call, and
# for 100 one-sample alignments that is ~100 index builds against 100 alignments that each take
# seconds -- the first attempt managed one sample in the time this whole job should take. A
# prebuilt `.mmi` is read straight in. Same trap as any per-item setup cost inside a loop.
MMI=$ROOT/ref/GRCh38.sr.mmi
if [ ! -s "$MMI" ]; then
    echo "=== building minimap2 index once (-x sr) ==="
    time minimap2 -x sr -t "$THREADS" -d "$MMI" "$REF"
fi
ls -la "$MMI" | awk '{print "index:", $5, "bytes"}'

# --- one sample: consensus -> BAM -> per-target molecules -> variants ------------------------------
process() {
    local acc=$1 out=$WORK/$acc
    mkdir -p "$out"

    if [ ! -f "$out/$acc.consensus.fq.gz" ]; then
        # Reuse phase 1's consensus when it exists, otherwise make it.
        if [ -f "$ROOT/panel/$acc/as/$acc.consensus.fq.gz" ]; then
            cp "$ROOT/panel/$acc/as/$acc.consensus.fq.gz" "$out/$acc.consensus.fq.gz"
        else
            local fq
            fq=$(fastq_for "$acc") || { echo "  no fastq for $acc" >&2; return 1; }
            migec checkout "$fq" --bc-pattern '0:12' --sample "$acc" \
                  -o "$out/co" -t "$THREADS" > "$out/checkout.log" 2>&1
            migec refine   "$out/co/$acc.fq.gz" -o "$out/rf" -t "$THREADS" > "$out/refine.log" 2>&1
            migec assemble "$out/rf/$acc.fq.gz" -o "$out/as" -t "$THREADS" > "$out/assemble.log" 2>&1
            mv "$out/as/$acc.consensus.fq.gz" "$out/$acc.consensus.fq.gz"
            rm -rf "$out/co" "$out/rf" "$out/as"      # the consensus is the only keeper
        fi
    fi

    if [ ! -f "$out/$acc.bam.bai" ]; then
        # The .mmi, not the .fa -- see the note above.
        minimap2 -ax sr -y -t "$THREADS" "$MMI" "$out/$acc.consensus.fq.gz" 2>/dev/null \
            | samtools sort -@ "$THREADS" -o "$out/$acc.bam" -
        samtools index "$out/$acc.bam"
    fi

    # Molecules per target. Never count READS here: after assemble one record is one molecule, so
    # counting rows would count nothing new, and counting distinct MI is robust to a molecule
    # producing a supplementary alignment.
    : > "$out/molecules_per_target.tsv"
    while read -r chrom start end name; do
        n=$(samtools view -F 0x900 "$out/$acc.bam" "${chrom}:$((start+1))-${end}" 2>/dev/null \
            | grep -o 'MI:Z:[^[:space:]]*' | sort -u | wc -l | tr -d ' ')
        printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$acc" "$chrom" "$start" "$end" "${name:-.}" "$n" \
            >> "$out/molecules_per_target.tsv"
    done < "$PANEL"

    # LoFreq: a STANDARD caller, because after assemble its depth already IS a molecule count.
    # Never a UMI-aware caller here -- it would group and consense a second time.
    if [ ! -f "$out/$acc.vcf.gz" ]; then
        lofreq call -f "$REF" -l "$PANEL" -o "$out/$acc.vcf" "$out/$acc.bam" 2>/dev/null || true
        [ -f "$out/$acc.vcf" ] && { bgzip -f "$out/$acc.vcf"; tabix -f -p vcf "$out/$acc.vcf.gz"; }
    fi
}

# --- run over the cohort ---------------------------------------------------------------------------
mapfile -t ACCS < <(ls "$FASTQ"/*.fastq.gz 2>/dev/null | xargs -n1 basename \
                    | sed 's/_1\.fastq\.gz$//; s/\.fastq\.gz$//' | sort -u)
[ "$LIMIT" -gt 0 ] && ACCS=("${ACCS[@]:0:$LIMIT}")
echo "processing ${#ACCS[@]} runs"

i=0
for acc in "${ACCS[@]}"; do
    i=$((i+1))
    echo "[$i/${#ACCS[@]}] $acc"
    process "$acc" || echo "  FAILED $acc"
done

# --- collect ----------------------------------------------------------------------------------------
{ printf "run\tchrom\tstart\tend\tgene\tmolecules\n"
  cat "$WORK"/*/molecules_per_target.tsv 2>/dev/null; } > "$ROOT/molecules_per_target.tsv"
echo "wrote $ROOT/molecules_per_target.tsv ($(( $(wc -l < "$ROOT/molecules_per_target.tsv") - 1 )) rows)"

{ printf "run\tchrom\tpos\tref\talt\tdp\taf\n"
  for v in "$WORK"/*/*.vcf.gz; do
      [ -f "$v" ] || continue
      acc=$(basename "$v" .vcf.gz)
      bcftools query -f "${acc}\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n" "$v" 2>/dev/null
  done; } > "$ROOT/variants.tsv"
echo "wrote $ROOT/variants.tsv ($(( $(wc -l < "$ROOT/variants.tsv") - 1 )) calls)"

echo "--- molecules per target, summed over runs ---"
awk -F'\t' 'NR>1{k=$2":"$3"-"$4" "$5; s[k]+=$6} END{for (x in s) printf "%-40s %10d\n", x, s[x]}' \
    "$ROOT/molecules_per_target.tsv" | sort -k2 -nr

df -h /projects/tcr_bcr_rnaseq | tail -1
echo done
