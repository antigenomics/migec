#!/bin/bash -l
#SBATCH --job-name=migec-mageri
#SBATCH --partition=medium
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.out
#SBATCH --error=/projects/tcr_bcr_rnaseq/migec_ctdna/logs/%x_%j.err
#
# MAGERI head to head, at the VARIANT level as well as the consensus level.
#
# MAGERI assembles, maps and calls in one run, so stopping the comparison at the consensus stops it
# one step short of what either tool is for. The simulator's `variant_af` turns the clone set into
# one reference plus point variants of it at a known allele fraction, which is the ctDNA shape and
# the only shape in which a call set has an answer.
#
# It runs here rather than on a laptop because LoFreq is here: migec has no caller of its own, so
# the comparison needs the one migec would hand its consensus to.
#
# Bundle expected at $HOME/migec-mageri: scripts/, tests/, mageri.jar (see `aldan3 push`).

set -euo pipefail

ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna
BUNDLE=$HOME/migec-mageri
VENV=$ROOT/venv-migec
OUT=$ROOT/mageri
THREADS=${SLURM_CPUS_PER_TASK:-8}
MIGEC_VERSION=${MIGEC_VERSION:-2.4.0}

# The tool binaries come from the environment the ctDNA work already used, so LoFreq, minimap2 and
# samtools are the same builds every other number here was produced with. Only migec is fresh:
# `$ROOT/env` holds 2.2.1, and a head to head has to be against the version being released.
export PATH="$ROOT/env/bin:$ROOT/env-gatk/bin:$PATH"
mkdir -p "$OUT" "$ROOT/logs"

if [ ! -x "$VENV/bin/python" ]; then
    "$ROOT/env/bin/python" -m venv "$VENV"
fi
"$VENV/bin/pip" install --quiet --upgrade pip
"$VENV/bin/pip" install --quiet "migec==$MIGEC_VERSION"
"$VENV/bin/python" -c "import migec; print('migec', migec.__version__)"
java -version 2>&1 | head -1
lofreq version 2>&1 | head -1
minimap2 --version

cd "$BUNDLE"
# 20,000 molecules over one reference plus five point variants, crossed with the depth.
#
# Never: 8 reads per UMI does not discriminate. Both tools call 5 of 5 with no false positive at
# 1% AND at 0.1% there, because a consensus over 8 reads removes essentially all sequencing error
# and the answer is then set by the molecule count, which is the same for both. The regime where a
# caller can be wrong is the SHALLOW one -- 1.5 reads per UMI, where most consensuses are one read
# and carry raw read quality straight into the pileup. That is the artifact-limited regime of
# `docs/detection.rst`, and it is the common one in practice.
for cov in 8 1.5; do
    for af in 0.01 0.001; do
        tag="${af/0./}_cov${cov/./}"
        "$VENV/bin/python" scripts/compare_mageri.py \
            --out "$OUT/af$tag" --jar "$BUNDLE/mageri.jar" \
            --molecules 20000 --clones 6 --coverage "$cov" --min-count 1 \
            --variant-af "$af" --caller lofreq --threads "$THREADS" \
            --tsv "$OUT/consensus_af$tag.tsv" --tsv-variants "$OUT/variants_af$tag.tsv"
    done
done

echo "=== consensus ==="
cat "$OUT"/consensus_af*.tsv
echo "=== variants ==="
cat "$OUT"/variants_af*.tsv
echo done
