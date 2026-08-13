#!/usr/bin/env python3
# 2026-08-14
# What is the lowest allele frequency this library can actually detect?
#
# The question every rare-variant application arrives with -- exome, ctDNA, MRD -- and it is
# answered by two numbers, neither of which is the variant caller:
#
#   N   molecules covering the site        (what `migec assemble` counts)
#   p   per-MOLECULE error floor at that base  (the RT/first-cycle floor, `docs/quality_floor.rst`)
#
# Sequencing deeper raises reads per molecule, not N. Only more input DNA, or more sites, raises
# the evidence. And p is a floor no consensus can go below: an error made before the barcode was
# attached is in every read of that molecule.
#
# Usage:
#   python scripts/detection_limit.py --molecules 12000                   # one site
#   python scripts/detection_limit.py --input-ng 20 --sites 5            # from DNA mass
#   python scripts/detection_limit.py --molecules 5000 --sites 30 --rt-error high   # MRD, duplex
#   python scripts/detection_limit.py --from-json asm/assemble.json --sites 5

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

PG_PER_HAPLOID_GENOME = 3.3

# The pre-amplification floor by protocol class, the same names `migec assemble --rt-error` takes.
# These are per base per molecule; `docs/quality_floor.rst` and `SOURCES.md` carry the citations.
FLOORS = {
    "rt": 1e-4,        # anything with a reverse transcription step. 10x state 1e-4 for V(D)J RT
    "medium": 1e-5,    # an ordinary polymerase, no RT
    "high": 1e-6,      # a proofreading polymerase, no RT
    "duplex": 1e-9,    # both strands must agree; the product of two independent floors
}


def poisson_at_least(lam: float, k: int) -> float:
    """P(X >= k) for X ~ Poisson(lam). k is small, so sum from below."""
    if lam <= 0:
        return 0.0 if k > 0 else 1.0
    return 1 - sum(math.exp(-lam) * lam**i / math.factorial(i) for i in range(k))


def limit_of_detection(molecules: float, min_support: int, confidence: float) -> float:
    """Lowest VAF detected with `confidence` probability, given `molecules` covering the site.

    Solves P(Poisson(N f) >= k) = confidence for f. Monotone in f, so bisect.
    """
    if molecules <= 0:
        return float("inf")
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = (lo + hi) / 2
        if poisson_at_least(molecules * mid, min_support) < confidence:
            lo = mid
        else:
            hi = mid
    return hi


def background_molecules(molecules: float, floor: float) -> float:
    """Expected FALSE variant molecules at one site, from the pre-amplification floor alone.

    Never: divide the floor by 3. `p` is the chance the base is wrong; a caller asks about ONE
    alternative allele, so only a third of those errors look like the variant being tracked.
    Using the whole floor overstates the background threefold and makes every LOD look worse
    than it is.
    """
    return molecules * floor / 3


def describe(molecules: float, sites: int, floor: float, min_support: int,
             confidence: float) -> dict:
    """The arithmetic for one configuration, per site and pooled over sites.

    Note: pooling over sites is what makes MRD work. Tracking `sites` independent patient-specific
    variants multiplies the molecules that can carry evidence, which is why an MRD panel follows
    tens of mutations rather than one -- see the note in the report.
    """
    pooled = molecules * sites
    return {
        "molecules_per_site": molecules,
        "sites": sites,
        "pooled_molecules": pooled,
        "lod_one_site": limit_of_detection(molecules, min_support, confidence),
        "lod_pooled": limit_of_detection(pooled, min_support, confidence),
        "background_per_site": background_molecules(molecules, floor),
        "background_pooled": background_molecules(pooled, floor),
        # Below this frequency the true signal is smaller than the chemistry's own false signal,
        # whatever the depth and whatever the caller. It is the floor's floor.
        "vaf_equals_background": floor / 3,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--molecules", type=float, help="molecules covering ONE site")
    src.add_argument("--input-ng", type=float, help="DNA input in ng; molecules are derived")
    src.add_argument("--from-json", type=Path, help="an assemble.json; uses its molecule count")

    ap.add_argument("--sites", type=int, default=1,
                    help="independent sites tracked. MRD pools evidence across them (default 1)")
    ap.add_argument("--efficiency", type=float, default=1.0,
                    help="fraction of input molecules that reach the consensus, with --input-ng")
    ap.add_argument("--strands", type=float, default=2.0,
                    help="strands recovered per input duplex fragment, with --input-ng (default 2)")
    ap.add_argument("--rt-error", default="rt",
                    help="floor: rt | medium | high | duplex, or a rate (default rt = 1e-4)")
    ap.add_argument("--min-support", type=int, default=3,
                    help="variant molecules a caller needs before it will call (default 3)")
    ap.add_argument("--confidence", type=float, default=0.95,
                    help="detection probability the LOD is quoted at (default 0.95)")
    a = ap.parse_args()

    floor = FLOORS.get(a.rt_error)
    if floor is None:
        try:
            floor = float(a.rt_error)
        except ValueError:
            raise SystemExit(f"--rt-error must be one of {list(FLOORS)} or a rate") from None

    if a.molecules is not None:
        molecules = a.molecules
        provenance = f"{molecules:,.0f} molecules/site, given"
    elif a.input_ng is not None:
        molecules = a.input_ng * 1000 / PG_PER_HAPLOID_GENOME * a.strands * a.efficiency
        provenance = (f"{a.input_ng:g} ng / 3.3 pg x {a.strands:g} strands x "
                      f"{a.efficiency:g} efficiency = {molecules:,.0f} molecules/site")
    else:
        summary = json.loads(a.from_json.read_text())
        total = summary.get("molecules") or summary.get("consensuses") or 0
        molecules = total / max(a.sites, 1)
        provenance = (f"{total:,.0f} molecules in {a.from_json.name} / {a.sites} sites = "
                      f"{molecules:,.0f} per site")
        # Never: a library total divided by the site count is an AVERAGE, and real panels are not
        # evenly covered -- measured 0.31-0.61x of the mean at the weakest target, plus off-target
        # product that a total cannot see. Align and count per target when a reference exists.
        print("note: --from-json divides a library total by --sites, which assumes even coverage.\n"
              "      Measured on a real panel the weakest target held 0.31-0.61x of the mean, and\n"
              "      off-target product was 8-58% of the library. Align and count per target if\n"
              "      you have a reference; otherwise read this as an optimistic average.\n",
              file=sys.stderr)

    r = describe(molecules, a.sites, floor, a.min_support, a.confidence)

    print(f"molecules      {provenance}")
    print(f"error floor    {floor:.1e} per base per molecule ({a.rt_error})")
    print(f"calling        >= {a.min_support} variant molecules, quoted at {a.confidence:.0%} "
          f"detection\n")

    print(f"{'':<26}{'one site':>16}{'pooled over ' + str(a.sites):>18}")
    print(f"{'molecules':<26}{r['molecules_per_site']:>16,.0f}{r['pooled_molecules']:>18,.0f}")
    print(f"{'limit of detection (VAF)':<26}{r['lod_one_site']:>16.2e}{r['lod_pooled']:>18.2e}")
    print(f"{'false molecules from floor':<26}{r['background_per_site']:>16.2f}"
          f"{r['background_pooled']:>18.2f}")

    print()
    lod = r["lod_pooled"] if a.sites > 1 else r["lod_one_site"]
    bg = r["vaf_equals_background"]
    if lod < bg:
        print(f"Never: the molecule count reaches {lod:.1e}, but the CHEMISTRY bottoms out at "
              f"{bg:.1e}.")
        print(f"       Below that a true variant is rarer than the floor's own false calls, so "
              f"more input\n       cannot help. Use a lower floor: --rt-error "
              f"{'high' if floor > 1e-6 else 'duplex'} (duplex needs both strands sequenced).")
    else:
        need = a.min_support / max(bg, 1e-12)
        print(f"Detectable to {lod:.1e}; the chemistry's own floor sits at {bg:.1e}, so the "
              f"molecule count is\nthe binding constraint. To go lower, add input DNA or track "
              f"more sites -- {need:,.0f} pooled\nmolecules would reach the floor.")

    if a.sites > 1:
        print(f"\nNote: pooling {a.sites} sites is what MRD does. Tracking a patient's own variant "
              f"set multiplies\nthe molecules that can carry evidence; one site at this input "
              f"reaches only {r['lod_one_site']:.1e}.")


if __name__ == "__main__":
    main()
