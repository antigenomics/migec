#!/bin/bash -l
#SBATCH --job-name=migec-panel
#SBATCH --partition=medium
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# Phase 1: infer the panel from the data, rather than trusting a vendor BED.
#
# The amplicon count in `assets/ctdna_titration.tsv` came from a consensus-prefix tally, which
# assumes the panel is EVENLY COVERED. It is not -- observed shares are 20.4/17.4/16.9/15.9/7.6%,
# so `molecules / n_amplicons` overstates the weakest target by 2.6x. Fixing that needs real
# coordinates: align the consensus to GRCh38, merge the covered intervals, and use the resulting
# BED as the panel. Then a per-target molecule count is a count, not an average.
#
# Never: the reference here is ENSEMBL-STYLE (`>1`, not `>chr1`). Every BED, region string and
# downstream tool has to use the same convention or it silently matches nothing.

set -euo pipefail

REF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa
GTF=/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.110.chr.gtf
ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
FASTQ=$ROOT/fastq
WORK=$ROOT/panel
ENVDIR=$ROOT/env
THREADS=${SLURM_CPUS_PER_TASK:-16}

mkdir -p "$WORK" "$ROOT/logs"
echo "host=$(hostname) job=${SLURM_JOB_ID:-local} threads=$THREADS"
df -h /projects/tcr_bcr_rnaseq | tail -1

# --- environment ---------------------------------------------------------------------------------
export MAMBA_ROOT_PREFIX=$ROOT/mamba
MICROMAMBA=$HOME/bin/micromamba
if [ ! -x "$ENVDIR/bin/minimap2" ]; then
    echo "creating env"
    "$MICROMAMBA" create -y -p "$ENVDIR" -c conda-forge -c bioconda \
        python=3.12 pip minimap2 samtools bedtools bcftools lofreq
fi
export PATH="$ENVDIR/bin:$PATH"
python -m pip install --quiet "migec==2.2.1"
echo "minimap2 $(minimap2 --version)  samtools $(samtools --version | head -1 | awk '{print $2}')"
python -c "import migec; print('migec', migec.__version__)"

# --- pick the deepest runs to infer the panel from -------------------------------------------------
# A panel is a property of the assay, so it only needs the runs with the most molecules; the
# shallow ones would just add noise to the interval edges.
PANEL_RUNS="SRR17220895 SRR17220892 SRR17220901 SRR8255842 SRR8255850"

consensus_for() {
    local acc=$1 out=$WORK/$acc
    if [ -f "$out/as/$acc.consensus.fq.gz" ]; then echo "$out/as/$acc.consensus.fq.gz"; return; fi
    mkdir -p "$out"
    migec checkout "$FASTQ/${acc}_1.fastq.gz" --bc-pattern '0:12' --sample "$acc" \
          -o "$out/co" -t "$THREADS" > "$out/checkout.log" 2>&1
    migec refine   "$out/co/$acc.fq.gz" -o "$out/rf" -t "$THREADS" > "$out/refine.log" 2>&1
    migec assemble "$out/rf/$acc.fq.gz" -o "$out/as" -t "$THREADS" > "$out/assemble.log" 2>&1
    echo "$out/as/$acc.consensus.fq.gz"
}

echo "=== consensus for the panel runs ==="
CONS=""
for acc in $PANEL_RUNS; do
    [ -f "$FASTQ/${acc}_1.fastq.gz" ] || { echo "  skip $acc: no fastq"; continue; }
    c=$(consensus_for "$acc")
    echo "  $acc -> $(zcat "$c" | wc -l | awk '{print $1/4}') molecules"
    CONS="$CONS $c"
done
[ -n "$CONS" ] || { echo "no consensus produced" >&2; exit 1; }

# --- align and infer the covered intervals ---------------------------------------------------------
# `-y` carries the molecule tags (RX/MI/CB) into the BAM; migec's own downstream contract.
echo "=== aligning to GRCh38 ==="
cat $CONS > "$WORK/panel_pool.fq.gz"
minimap2 -ax sr -y -t "$THREADS" "$REF" "$WORK/panel_pool.fq.gz" 2> "$WORK/minimap2.log" \
    | samtools sort -@ "$THREADS" -o "$WORK/panel_pool.bam" -
samtools index "$WORK/panel_pool.bam"
samtools flagstat "$WORK/panel_pool.bam" | head -6

# An amplicon panel piles molecules on a handful of loci and leaves the rest of the genome empty,
# so a depth floor separates target from off-target cleanly. The floor is a fraction of the pooled
# molecule count rather than a constant: it has to scale with how much was pooled.
TOTAL=$(samtools view -c -F 0x900 "$WORK/panel_pool.bam")
FLOOR=$(python - <<PY
print(max(50, int($TOTAL * 0.002)))
PY
)
echo "=== inferring panel: $TOTAL primary alignments, depth floor $FLOOR ==="

bedtools genomecov -ibam "$WORK/panel_pool.bam" -bga \
    | awk -v f="$FLOOR" '$4 >= f' \
    | bedtools merge -d 50 -i - \
    | awk '$3 - $2 >= 40' > "$WORK/panel_inferred.bed"

echo "intervals found: $(wc -l < "$WORK/panel_inferred.bed")"
awk '{printf "%s:%d-%d\t%d bp\n", $1, $2, $3, $3-$2}' "$WORK/panel_inferred.bed"

# --- name the intervals from the GTF ---------------------------------------------------------------
# Which genes these are is the thing a prefix tally could never tell us.
echo "=== annotating against Ensembl 110 ==="
awk -F'\t' '$3=="gene"' "$GTF" \
    | sed -n 's/^\([^\t]*\)\t[^\t]*\tgene\t\([0-9]*\)\t\([0-9]*\).*gene_name "\([^"]*\)".*/\1\t\2\t\3\t\4/p' \
    | sort -k1,1 -k2,2n > "$WORK/genes.bed"

bedtools intersect -a "$WORK/panel_inferred.bed" -b "$WORK/genes.bed" -wa -wb -loj \
    | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"($7=="."?"intergenic":$7)}' \
    | sort -u -k1,1 -k2,2n > "$WORK/panel_named.bed"

echo "--- the inferred panel ---"
awk '{printf "%-6s %12s-%-12s %5d bp  %s\n", $1, $2, $3, $3-$2, $4}' "$WORK/panel_named.bed"

echo
echo "wrote $WORK/panel_named.bed"
df -h /projects/tcr_bcr_rnaseq | tail -1
echo done
