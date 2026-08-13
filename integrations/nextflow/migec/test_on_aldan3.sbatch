#!/bin/bash
#SBATCH --job-name=migec-nf-test
#SBATCH --partition=short
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=migec-nf-test-%j.out
#SBATCH --error=migec-nf-test-%j.err
#
# Validate the MIGEC Nextflow module on Aldan-3.
#
# Never run this on a frontend. Frontends are shared and have two cores; the aldan3 client refuses
# to start nextflow there. Submit it:
#
#     aldan3 push integrations/nextflow/migec /home/<you>/migec-nf
#     aldan3 slurm submit /home/<you>/migec-nf/test_on_aldan3.sbatch
#     aldan3 slurm log <jobid>
#
# Four steps, in increasing order of what they need:
#   1. what the node actually has
#   2. a tiny simulated library, so the test needs no input data
#   3. `-stub-run` -- validates the module's syntax and wiring, and needs only nextflow
#   4. the real commands, which need migec installed
#
# Steps 1-3 run whether or not migec is available, because they are the ones a cluster can always
# check and the module's wiring is the part most likely to be wrong.

set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
# Not $BASH_SOURCE: SLURM runs a COPY of the script out of its spool directory, so the script's
# own path says nothing about where the module was pushed. Default to the documented push
# destination and let an env var override it.
MODULE_DIR="${MIGEC_MODULE_DIR:-$HOME/migec-nf}"
WORK="${TMPDIR:-/tmp}/migec-nf-${SLURM_JOB_ID:-local}"
mkdir -p "$WORK"
trap 'rm -rf "$WORK"' EXIT

echo "== node"
hostname; nproc; free -g | head -2
echo "module dir: $MODULE_DIR"
test -f "$MODULE_DIR/main.nf" || echo "  WARNING: no main.nf there -- set MIGEC_MODULE_DIR"
for tool in nextflow java podman python3; do
    printf '%-9s %s\n' "$tool" "$(command -v "$tool" || echo MISSING)"
done

echo
echo "== a tiny library, simulated here so the test needs no input data"
python3 - "$WORK" <<'PY'
import gzip, pathlib, random, sys
out = pathlib.Path(sys.argv[1]); rng = random.Random(0)
adapter = "CAGTGGTATCAACGCAGAGT"
with gzip.open(out / "r1.fq.gz", "wt") as fh:
    for m in range(2000):
        umi = "".join(rng.choice("ACGT") for _ in range(12))
        payload = "".join(rng.choice("ACGT") for _ in range(70))
        for r in range(4):
            s = umi + adapter + payload
            fh.write(f"@m{m}_{r}\n{s}\n+\n{'I' * len(s)}\n")
(out / "bc.txt").write_text(f"S1\t{'N' * 12}{adapter.lower()}\n")
print("wrote", out / "r1.fq.gz")
PY

echo
echo "== stub run: the module's syntax and wiring, which needs only nextflow"
if command -v nextflow >/dev/null; then
    mkdir -p "$WORK/modules"
    cp -r "$MODULE_DIR" "$WORK/modules/migec"
    # Aldan-3 ships Nextflow 21.10.6, where DSL2 is opt-in. It is the default from 22.x, but
    # declaring it costs nothing and the module is DSL2 either way.
    cat > "$WORK/main.nf" <<'NF'
nextflow.enable.dsl = 2

include { MIGEC } from './modules/migec/main'
workflow {
    Channel.of([[id: 'S1'], [file(params.r1)]]) | MIGEC
    MIGEC.out.consensus.view { meta, f -> "consensus: ${meta.id} -> ${f.name}" }
}
NF
    ( cd "$WORK" && NXF_ANSI_LOG=false timeout 600 nextflow run main.nf -stub-run \
        --r1 "$WORK/r1.fq.gz" \
        --migec_bc_pattern 'NNNNNNNNNNNNcagtggtatcaacgcagagt' \
        --outdir "$WORK/nf_out" -work-dir "$WORK/work" ) \
        && echo "STUB OK" || echo "STUB FAILED"
else
    echo "nextflow MISSING on this node -- cannot validate the module here"
fi

echo
echo "== migec"
if ! command -v migec >/dev/null; then
    # Never a bare `pip install`: a compute node may have no route to PyPI, and a job that stalls
    # for its whole time limit reports nothing at all.
    timeout 120 python3 -m pip install --quiet --user --timeout 10 --retries 1 'migec>=2.0.0' \
        || echo "pip could not install migec (no route to PyPI, or not yet published)"
fi
if ! command -v migec >/dev/null; then
    echo "SKIP: migec is not on PATH and could not be installed."
    echo "      Install it into a cluster conda env and re-submit; everything above still stands."
    exit 0
fi
migec info

echo
echo "== the commands the module wraps"
set -e
migec checkout "$WORK/r1.fq.gz" -b "$WORK/bc.txt" -o "$WORK/co" -t "${SLURM_CPUS_PER_TASK:-4}"
migec refine "$WORK/co/S1.fq.gz" -o "$WORK/ref"
migec assemble "$WORK/ref/S1.fq.gz" -o "$WORK/asm"
test -s "$WORK/asm/S1.consensus.fq.gz" && echo "PIPELINE OK" || { echo "PIPELINE FAILED"; exit 1; }

echo "== done"
