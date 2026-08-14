#!/bin/bash -l
#SBATCH --job-name=migec-callers
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# Is the artifact the CALLER's, or the consensus input's?
#
# `ctdna_minreads.sbatch` scored LoFreq against the certified arms and found that the 0% arm
# returns 9-11 calls per sample at 0.4-1.4% VAF, 94% of them -> G, eight positions reproducing in
# 3 of 3 replicates. That is a 2-colour dark-G artifact, and it is a property of the reads that
# reached the consensus -- but ONE caller cannot prove that. A second caller on the IDENTICAL BAM
# can: if Mutect2 calls the same positions, the artifact is in the input; if it does not, it was
# LoFreq's threshold.
#
# So this job changes exactly one thing at a time:
#   - same consensus BAM, LoFreq vs Mutect2                     -> caller vs caller
#   - same caller, with and without a recurrence filter         -> what a background model buys
#   - same caller, with and without a MAPQ 20 floor             -> what the mismapping costs
# and it reuses the BAMs `ctdna_minreads.sbatch` already wrote, so nothing is re-assembled and the
# comparison cannot drift by a pipeline change between arms.
#
# Never: UMI-VarCal and UMIErrorCorrect are NOT here. They do their own UMI grouping and their own
# consensus, so they REPLACE assemble rather than following it -- running them on this BAM would
# collapse an already-collapsed library and report a molecule count that means nothing. They are a
# pipeline-level comparison from raw reads, which is a different job.

set -euo pipefail

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
PANEL=$ROOT/panel/panel_named.bed
ENVDIR=$ROOT/env
GATKDIR=$ROOT/env-gatk
WORK=$ROOT/callers
THREADS=${SLURM_CPUS_PER_TASK:-16}
export PATH="$ENVDIR/bin:$PATH"
mkdir -p "$WORK" "$ROOT/logs"

# The same twelve runs `ctdna_minreads.sbatch` used: one replicate per certified arm at the input
# and depth with the most molecules, so a difference between callers is not a difference in how
# much evidence they were given.
declare -A ARMS=(
  [SRR17220921]=WT [SRR17220920]=WT [SRR17220919]=WT
  [SRR17220931]=0.125 [SRR17220930]=0.125 [SRR17220929]=0.125
  [SRR17220928]=0.25 [SRR17220926]=0.25 [SRR17220925]=0.25
  [SRR17220924]=1 [SRR17220923]=1 [SRR17220922]=1
)

echo "host=$(hostname) job=${SLURM_JOB_ID:-local} threads=$THREADS"

# GNU time, for the wall clock and peak RSS of each caller. Absent, the calls still run and the
# runtime table is simply empty -- a missing measurement, never a made-up one.
TIME=/usr/bin/time
[ -x "$TIME" ] || TIME=""

# A wrapper rather than `${TIME:+...}`: the conditional expansion is not re-quoted, so the quotes
# around a format string end up IN the format string, and dropping them loses the `\t` to bash's
# own quote removal. Both failures are silent -- a runtime table full of the letter t.
run_timed() {
    local label=$1 acc=$2 mr=$3 mapq=$4 tf=$5
    shift 5
    if [ -n "$TIME" ]; then
        "$TIME" -f "$label\t$acc\t$mr\t$mapq\t%e\t%M" -o "$tf" "$@"
    else
        "$@"
    fi
}

# --- gatk in its OWN prefix ------------------------------------------------------------------------
# Never: not into $ROOT/env. gatk4 drags a JDK and a second python along, and re-solving the
# environment that produced every number so far to add a comparator is how the comparator changes
# the baseline. A separate prefix cannot.
export MAMBA_ROOT_PREFIX=$ROOT/mamba
MICROMAMBA=$HOME/bin/micromamba
if [ ! -x "$GATKDIR/bin/gatk" ]; then
    echo "=== creating the gatk environment ==="
    "$MICROMAMBA" create -y -p "$GATKDIR" -c conda-forge -c bioconda gatk4 openjdk=17
fi
GATK=$GATKDIR/bin/gatk
"$GATK" --version | head -3

# --- the reference GATK wants ----------------------------------------------------------------------
# Mutect2 needs a .fai and a .dict beside the fasta, and the reference tree is not ours to write in,
# so the fasta is linked here and the sidecars are built next to the link.
mkdir -p "$WORK/ref"
LOCALREF=$WORK/ref/GRCh38.fa
[ -e "$LOCALREF" ] || ln -s "$REF" "$LOCALREF"
[ -s "$LOCALREF.fai" ] || { [ -s "$REF.fai" ] && ln -sf "$REF.fai" "$LOCALREF.fai" || samtools faidx "$LOCALREF"; }
[ -s "$WORK/ref/GRCh38.dict" ] || "$GATK" CreateSequenceDictionary -R "$LOCALREF" -O "$WORK/ref/GRCh38.dict"

# --- one (run, min_reads, arm): every caller and every filter on ONE bam ----------------------------
call_all() {
    local acc=$1 mr=$2 out=$WORK/$acc/mr$mr
    local src=$ROOT/minreads/$acc/mr$mr/$acc.bam
    [ -s "$src" ] || { echo "  no bam for $acc mr$mr" >&2; return 1; }
    mkdir -p "$out"

    # Mutect2 refuses a BAM with no read group, and minimap2 writes none. Adding one rewrites the
    # header and the RG tag only -- no realignment, so the alignments the two callers see are the
    # same bytes.
    if [ ! -s "$out/rg.bam" ]; then
        # printf, not "\t" in double quotes: bash leaves the backslash alone there and the header
        # line arrives as one field, which samtools accepts and Mutect2 then reads as no sample.
        local rg
        rg=$(printf '@RG\tID:1\tSM:%s\tPL:ILLUMINA\tLB:%s' "$acc" "$acc")
        samtools addreplacerg -r "$rg" -o "$out/rg.bam" "$src"
        samtools index "$out/rg.bam"
    fi
    # The MAPQ 20 arm. Adapter read-through mismaps into a 97%-G locus at MAPQ 4-16, so this is the
    # cheap half of the fix: no re-run, just refuse the alignments that were never trustworthy.
    if [ ! -s "$out/q20.bam" ]; then
        samtools view -b -q 20 -o "$out/q20.bam" "$out/rg.bam"
        samtools index "$out/q20.bam"
    fi

    for bam in rg q20; do
        local mapq=0; [ "$bam" = q20 ] && mapq=20
        if [ ! -s "$out/lofreq.$bam.vcf.gz" ]; then
            run_timed lofreq "$acc" "$mr" "$mapq" "$out/lofreq.$bam.time" \
                lofreq call -f "$LOCALREF" -l "$PANEL" -o "$out/lofreq.$bam.vcf" \
                "$out/$bam.bam" 2>/dev/null || true
            [ -f "$out/lofreq.$bam.vcf" ] && bgzip -f "$out/lofreq.$bam.vcf"
        fi
        if [ ! -s "$out/mutect2.$bam.vcf.gz" ]; then
            # Never: `--max-reads-per-alignment-start 0`. Mutect2 DOWNSAMPLES to 50 reads per
            # alignment start by default, and on an amplicon consensus every molecule starts at
            # the same coordinate -- so it kept 78-207 of the 4,282-15,895 molecules LoFreq saw on
            # the identical BAM, and a 1% variant was then supported by one or two of them. It
            # reported ~1 call per sample and never the certified hotspot, which reads as "Mutect2
            # cannot do this" and is really "Mutect2 was shown 1.5% of the evidence". Same class
            # of trap as a family-size filter after `assemble` (docs/variants.rst).
            run_timed mutect2 "$acc" "$mr" "$mapq" "$out/mutect2.$bam.time" \
                "$GATK" Mutect2 -R "$LOCALREF" -L "$PANEL" -I "$out/$bam.bam" \
                --max-reads-per-alignment-start 0 \
                -O "$out/mutect2.$bam.raw.vcf.gz" > "$out/mutect2.$bam.log" 2>&1 || true
            if [ -s "$out/mutect2.$bam.raw.vcf.gz" ]; then
                # Never: FilterMutectCalls is not optional. Mutect2's own authors require it, and
                # Maruzani's benchmark skipped it on purpose to keep six callers even -- so its
                # Mutect2 column is that tool without the filter it ships with.
                "$GATK" FilterMutectCalls -R "$LOCALREF" -V "$out/mutect2.$bam.raw.vcf.gz" \
                    -O "$out/mutect2.$bam.vcf.gz" >> "$out/mutect2.$bam.log" 2>&1 || true
            fi
        fi
    done
    echo "  done $acc mr$mr (${ARMS[$acc]})"
}

echo "=== calling: ${#ARMS[@]} runs x min-reads 1/3/5 x {lofreq, mutect2} x {all, MAPQ>=20} ==="
for acc in "${!ARMS[@]}"; do
    for mr in 1 3 5; do
        call_all "$acc" "$mr" || true
    done
done

# --- one table, one row per call --------------------------------------------------------------------
# Flat on purpose: run, arm, min_reads, caller, mapq_floor are five columns, never one label. The
# scoring in `scripts/ctdna_persite.py` groups on them.
echo "=== collecting ==="
{
  printf "run\tarm\tmin_reads\tcaller\tmapq_floor\tfiltered\tchrom\tpos\tref\talt\tdp\taf\n"
  for acc in "${!ARMS[@]}"; do
    for mr in 1 3 5; do
      for bam in rg q20; do
        mapq=0; [ "$bam" = q20 ] && mapq=20
        v=$WORK/$acc/mr$mr/lofreq.$bam.vcf.gz
        [ -s "$v" ] && bcftools query -f \
          "${acc}\t${ARMS[$acc]}\t${mr}\tlofreq\t${mapq}\tPASS\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n" \
          "$v" 2>/dev/null || true
        v=$WORK/$acc/mr$mr/mutect2.$bam.vcf.gz
        # Mutect2's PASS is FilterMutectCalls' verdict; the rest are kept with their reason so the
        # filter's own contribution is visible instead of assumed.
        [ -s "$v" ] && bcftools query -f \
          "${acc}\t${ARMS[$acc]}\t${mr}\tmutect2\t${mapq}\t%FILTER\t%CHROM\t%POS\t%REF\t%ALT\t[%DP]\t[%AF]\n" \
          "$v" 2>/dev/null || true
      done
    done
  done
} > "$ROOT/callers_variants.tsv"
echo "wrote $ROOT/callers_variants.tsv ($(( $(wc -l < "$ROOT/callers_variants.tsv") - 1 )) calls)"

{
  printf "caller\trun\tmin_reads\tmapq_floor\tseconds\tpeak_rss_kb\n"
  cat "$WORK"/*/mr*/*.time 2>/dev/null | grep -E "^(lofreq|mutect2)" || true
} > "$ROOT/callers_runtime.tsv"
echo "wrote $ROOT/callers_runtime.tsv"
echo done
