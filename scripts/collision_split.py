#!/usr/bin/env python3
# 2026-08-14
#
# Does `assemble`'s linkage sub-clustering actually separate two molecules that collided on a
# barcode, and from what depth?
#
# This is the other half of the map-first comparison in `docs/grouping.rst`. UMI-tools and fgbio
# separate collided molecules by sending them to different alignment positions; we have no aligner
# at that point, so the claim is that the payload does it instead. That claim needs a number, and
# it has one -- the linkage test's own ceiling. The strongest evidence a pair of columns can carry
# for a 50/50 split is `log10 C(n, n/2)`, which reaches X3's threshold of 8.68 at n ~ 32:
#
#     n         6    10    20    30    34    40    80
#     ceiling  1.3   2.4   5.3   8.2   9.4  11.1  23.0
#
# so below ~32 reads on the barcode nothing can clear it, whatever the two molecules look like.
# This measures where that bites, against the simulator's own record of which barcodes truly held
# more than one molecule.
#
# Never: "collided" is defined on the TRUE barcode, never the observed one. Two molecules whose
# OBSERVED barcodes coincide because one of them picked up a sequencing error are not a collision,
# they are what `refine` corrects -- and counting them makes the collision rate grow with the read
# count, which it cannot do. Measured before the definition was fixed: 16 "collisions" at 5 reads
# per molecule rising to 169 at 80, on a library whose molecule count never moved.
#
# Usage:
#     python scripts/collision_split.py --out /tmp/split --coverage 5 20 40 80 160

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

ADAPTER = "CAGTGGTATCAACGCAGAGT"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--molecules", type=int, default=20_000)
    ap.add_argument("--clones", type=int, default=200,
                    help="distinct sequences. Never: on ONE clone there is nothing to separate -- "
                         "two molecules that collide there hold the same sequence, so no amount of "
                         "sub-clustering and no aligner can tell them apart")
    ap.add_argument("--coverage", type=float, nargs="+", default=[5.0, 20.0, 40.0, 80.0, 160.0])
    ap.add_argument("--umi-len", type=int, default=12)
    ap.add_argument("--umi-error", type=float, default=3e-3)
    ap.add_argument("--seed", type=int, default=17)
    ap.add_argument("--tsv", type=Path)
    args = ap.parse_args(argv)

    from migec.assemble import run as assemble_run
    from migec.checkout import run as checkout_run
    from migec.refine import run as refine_run
    from tests.synthetic._sim import SimConfig, simulate

    cols = ["reads_per_molecule", "clones", "true_collisions", "split_by_assemble", "fraction",
            "groups_split", "mean_reads_on_collided_barcode"]
    lines = ["\t".join(cols)]
    print(lines[0])
    for cov in args.coverage:
        d = args.out / f"c{args.clones}_x{cov}"
        cfg = SimConfig(n_molecules=args.molecules, n_clones=args.clones, coverage=cov,
                        coverage_cv=0.4, umi_len=args.umi_len, umi_error=args.umi_error,
                        adapter=ADAPTER, seed=args.seed)
        sim = simulate(cfg, d / "sim")
        (d / "bc.txt").write_text(f"S1\t{sim['pattern']}\n")
        checkout_run(sim["reads"], d / "bc.txt", d / "co")
        refine_run(d / "co" / "S1.fq.gz", d / "ref")
        st = assemble_run(d / "ref" / "S1.fq.gz", d / "asm")

        by_umi: dict[str, set[str]] = collections.defaultdict(set)
        reads_per: collections.Counter = collections.Counter()
        with open(sim["truth_reads"]) as fh:
            h = fh.readline().rstrip("\n").split("\t")
            i_umi, i_mol = h.index("umi_true"), h.index("molecule_id")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                by_umi[f[i_umi]].add(f[i_mol])
                reads_per[f[i_umi]] += 1
        collided = {u for u, m in by_umi.items() if len(m) > 1}

        emitted: collections.Counter = collections.Counter()
        tsv = next(iter(sorted((d / "asm").glob("*.mig.tsv"))), None)
        if tsv is None:
            raise SystemExit(f"assemble wrote no per-molecule table into {d / 'asm'}")
        with open(tsv) as fh:
            i_umi = fh.readline().rstrip("\n").split("\t").index("umi")
            for line in fh:
                emitted[line.split("\t")[i_umi]] += 1

        split = sum(1 for u in collided if emitted.get(u, 0) > 1)
        mean_reads = sum(reads_per[u] for u in collided) / len(collided) if collided else 0.0
        row = (f"{cov}\t{args.clones}\t{len(collided)}\t{split}\t"
               f"{split / len(collided) if collided else 0:.3f}\t{st['groups_split']}\t"
               f"{mean_reads:.1f}")
        lines.append(row)
        print(row, flush=True)

    if args.tsv:
        args.tsv.write_text("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
