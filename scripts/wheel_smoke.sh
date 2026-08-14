#!/bin/bash
# What a freshly built wheel has to do before it is published. 2026-08-14
#
# Run by cibuildwheel (CIBW_TEST_COMMAND) inside the target environment, where the ONLY thing
# installed is the wheel and its declared dependencies -- so this also checks that the dependency
# list is complete, which importing the package does not.
#
# Never: `import migec` is not a smoke test. It leaves the console script, typer, and every code
# path that touches the extension at runtime unexercised, and all three have shipped broken behind
# a green import.
set -euo pipefail

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT
cd "$work"

python -c "import migec; print('version', migec.__version__)"
migec info
migec sheet --presets > /dev/null

# One molecule, four reads, through all three stages: an 8 nt UMI, a constant adapter, a payload.
adapter=CAGTGGTATCAACGCAGAGT
{
  for i in 1 2 3 4; do
    printf '@r%s\n%s\n+\n%s\n' "$i" "ACGTACGT${adapter}TTTTCCCCGGGGAAAA" "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
  done
} > reads.fq
printf 'S1\t%s\n' "NNNNNNNN${adapter}" | tr 'A-Z' 'a-z' | sed 's/nnnnnnnn/NNNNNNNN/' > bc.txt

migec checkout reads.fq -b bc.txt -o co
migec refine co/S1.fq.gz -o re
migec assemble re/S1.fq.gz -o as
# One molecule in, one consensus out. A wheel that gets this wrong is not shippable.
n=$(python -c "import gzip,sys; print(sum(1 for i,_ in enumerate(gzip.open('as/S1.consensus.fq.gz','rt')) ) // 4)")
test "$n" = "1" || { echo "expected 1 consensus record, got $n" >&2; exit 1; }
echo "wheel smoke: ok"
