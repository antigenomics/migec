"""A UMI read simulator with ground truth.

This is the single most valuable test asset in the repo: every accuracy claim in M1-M3 is measured
against it. It is deliberately NOT a shipped CLI command -- it lives in tests until the assertions
that use it have stabilised.

The model, in the order the real chemistry applies it:

  1. `n_molecules` template molecules are drawn from `n_clones` distinct sequences (so that
     collisions between molecules of the SAME sequence exist -- those are undetectable by any
     consensus method and are the reason molecule counts are biased low).
  2. Each molecule gets a UMI drawn from a possibly non-uniform base composition. With
     `n_molecules` comparable to 4^umi_len, genuine collisions happen at the birthday rate, which
     is exactly what the correction step must not mistake for errors.
  3. RT introduces errors at `rt_error` per base. These happen BEFORE any amplification, so they
     are present in every read of the molecule and no consensus can remove them. This is the floor
     the emitted quality must never claim to beat.
  4. PCR runs `pcr_cycles` cycles at efficiency `pcr_efficiency` with `pcr_error` per base per
     cycle. An error in an early cycle reaches a large fraction of the descendants; if it exceeds
     half, it becomes the consensus base. Errors in the UMI itself create error-child MIGs whose
     size ratio to the parent follows the branching process, not a Poisson -- which is why the
     correction model needs a polymerase component.
  5. Sequencing draws reads and adds errors at a per-base rate derived from the assigned quality.

Truth files written by `simulate()`:

  truth_reads.tsv       read_id, molecule_id, clone_id, umi_true, umi_observed, n_errors
  truth_molecules.tsv   molecule_id, clone_id, umi_true, n_reads, rt_errors, early_pcr_errors
  truth_consensus.fa    molecule_id -> the sequence a perfect assembler would report
  clones.fa             the clone sequences, as a reference for map-first tools to align to
                        (template + RT errors + any PCR error that reached >50% of the reads)
"""

from __future__ import annotations

import gzip
import random
from dataclasses import dataclass, field
from pathlib import Path

BASES = "ACGT"


@dataclass
class SimConfig:
    n_molecules: int = 1000
    n_clones: int = 20
    seq_len: int = 120
    umi_len: int = 12
    coverage: float = 10.0  # mean reads per molecule
    coverage_cv: float = 0.8  # lognormal spread; real MIG size distributions are very skewed
    seq_error: float = 1e-3
    rt_error: float = 1e-5
    pcr_error: float = 1e-5
    pcr_cycles: int = 25
    pcr_efficiency: float = 0.9
    umi_error: float = 3e-4  # per base, sequencing error in the UMI
    umi_base_freqs: tuple[float, float, float, float] = (0.25, 0.25, 0.25, 0.25)
    mean_qual: int = 32
    # Constant region between the UMI and the payload. Empty gives a bare `<umi><payload>` read,
    # which no barcode pattern can anchor on -- set it when the reads have to survive checkout.
    adapter: str = ""
    seed: int = 20260813
    # Write a second FASTQ holding the other mate: mate 1 is the first `read_len` payload bases
    # (behind the UMI and the adapter), mate 2 is the last `read_len` bases as SEQUENCED, which is
    # their reverse complement. The two are sequenced independently, so an error in the overlap is
    # in one mate and not the other -- which is the whole reason merging them is worth anything.
    # `read_len` 0 means the whole payload, so the mates overlap completely.
    paired: bool = False
    read_len: int = 0
    # Turn the clone set into a REFERENCE plus point variants of it, each at this allele fraction.
    # 0 leaves the default: `n_clones` unrelated random sequences, which is a repertoire. Anything
    # above 0 is the ctDNA shape -- one template, `n_clones - 1` single substitutions of it, each
    # in `variant_af` of the molecules -- and it is what makes a VARIANT CALLER comparable, since
    # a caller needs one reference and a known answer at known positions.
    variant_af: float = 0.0
    truth: dict = field(default_factory=dict)


def _rand_seq(rng: random.Random, n: int) -> str:
    return "".join(rng.choice(BASES) for _ in range(n))


def _mutate(seq: str, rate: float, rng: random.Random) -> tuple[str, int]:
    """Substitutions only. Illumina indel rates are ~1e-6/base and we model none anywhere."""
    if rate <= 0:
        return seq, 0
    out = list(seq)
    n = 0
    for i, b in enumerate(out):
        if rng.random() < rate:
            out[i] = rng.choice([x for x in BASES if x != b])
            n += 1
    return "".join(out), n


def _revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _draw_umi(rng: random.Random, length: int, freqs) -> str:
    return "".join(rng.choices(BASES, weights=freqs, k=length))


def _early_pcr_variants(cfg: SimConfig, rng: random.Random, seq_len: int) -> dict[int, str]:
    """Positions where a PCR error reached more than half the molecule's descendants.

    An error first appearing at cycle k ends up in a fraction f of the final population; under a
    Galton-Watson model with efficiency e, E[f] ~ 1/(1+e)^k. Only f > 0.5 flips the consensus, and
    that essentially requires k <= 1. Rather than simulating the whole tree per molecule (which is
    where a naive simulator spends all its time), we sample the cycle-1 and cycle-2 events
    directly -- the tail beyond that cannot reach 50% at any realistic efficiency.
    """
    variants: dict[int, str] = {}
    e = cfg.pcr_efficiency
    for cycle, share in ((1, 1.0 / (1.0 + e)), (2, 1.0 / (1.0 + e) ** 2)):
        expected = seq_len * cfg.pcr_error * (1.0 + e) ** (cycle - 1)
        n_events = _poisson(rng, expected)
        for _ in range(n_events):
            # The realised descendant fraction is stochastic; only take it if it crosses half.
            frac = rng.betavariate(max(share * 4.0, 0.1), max((1 - share) * 4.0, 0.1))
            if frac > 0.5:
                variants[rng.randrange(seq_len)] = rng.choice(BASES)
    return variants


def _poisson(rng: random.Random, lam: float) -> int:
    if lam <= 0:
        return 0
    if lam > 30:  # normal approximation; we never need this branch at realistic rates
        return max(0, int(rng.gauss(lam, lam**0.5)))
    import math

    limit = math.exp(-lam)
    k, p = 0, rng.random()
    while p > limit:
        p *= rng.random()
        k += 1
    return k


def simulate(cfg: SimConfig, out_dir: str | Path) -> dict:
    """Write reads + truth. Returns a summary dict."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(cfg.seed)

    clones = [_rand_seq(rng, cfg.seq_len) for _ in range(cfg.n_clones)]
    variants: list[tuple[int, int, str, str]] = []   # clone, 0-based pos, ref base, alt base
    clone_weights = None
    if cfg.variant_af > 0.0:
        clones, variants, clone_weights = _variant_clones(cfg, rng, clones[0])
    qual_char = chr(33 + cfg.mean_qual)

    reads_path = out / "reads.fq.gz"
    reads2_path = out / "reads2.fq.gz"
    read_len = cfg.read_len or cfg.seq_len
    truth_reads = out / "truth_reads.tsv"
    truth_mols = out / "truth_molecules.tsv"
    truth_cons = out / "truth_consensus.fa"
    # The clone sequences the molecules were drawn from, as a reference. Map-first UMI tools
    # (UMI-tools, fgbio) group on (position, UMI) and so need an alignment before they can be
    # compared with us at all; this is what they align to.
    truth_clones = out / "clones.fa"
    with open(truth_clones, "w") as cf:
        for i, seq in enumerate(clones):
            cf.write(f">clone{i}\n{seq}\n")

    # With variants, a caller needs ONE reference and the answer at known positions. Both are
    # written even though `clones.fa` holds the same sequence: pointing a caller at every variant
    # clone as well would let each variant read align to its own reference and be called nothing.
    reference = out / "reference.fa"
    truth_variants = out / "truth_variants.tsv"
    if variants:
        reference.write_text(f">clone0\n{clones[0]}\n")
        with open(truth_variants, "w") as vf:
            vf.write("chrom\tpos\tref\talt\texpected_af\n")
            for _, pos, ref_base, alt in variants:
                vf.write(f"clone0\t{pos + 1}\t{ref_base}\t{alt}\t{cfg.variant_af:.6g}\n")

    n_reads = 0
    n_collisions = 0
    seen_umis: dict[str, int] = {}

    fq2 = gzip.open(reads2_path, "wt") if cfg.paired else None
    with gzip.open(reads_path, "wt") as fq, open(truth_reads, "w") as tr, open(
        truth_mols, "w"
    ) as tm, open(truth_cons, "w") as tc:
        tr.write("read_id\tmolecule_id\tclone_id\tumi_true\tumi_observed\tn_errors\n")
        tm.write("molecule_id\tclone_id\tumi_true\tn_reads\trt_errors\tearly_pcr_errors\n")

        for mol in range(cfg.n_molecules):
            clone = (rng.randrange(cfg.n_clones) if clone_weights is None
                     else rng.choices(range(cfg.n_clones), weights=clone_weights, k=1)[0])
            umi = _draw_umi(rng, cfg.umi_len, cfg.umi_base_freqs)
            if umi in seen_umis:
                n_collisions += 1
            seen_umis[umi] = seen_umis.get(umi, 0) + 1

            # RT errors happen before amplification, so they are in every read of this molecule.
            template, n_rt = _mutate(clones[clone], cfg.rt_error, rng)
            early = _early_pcr_variants(cfg, rng, cfg.seq_len)
            consensus = list(template)
            for pos, base in early.items():
                consensus[pos] = base
            consensus_seq = "".join(consensus)

            # Lognormal MIG sizes: the empirical distribution is strongly right-skewed, and a
            # Poisson would make the 1-5 read regime -- the one the retention rule is about --
            # vanish.
            mu = max(0.0, cfg.coverage)
            size = max(1, int(rng.lognormvariate(_log_mu(mu, cfg.coverage_cv), cfg.coverage_cv)))

            tm.write(f"{mol}\t{clone}\t{umi}\t{size}\t{n_rt}\t{len(early)}\n")
            tc.write(f">mol{mol} clone={clone} umi={umi} reads={size}\n{consensus_seq}\n")

            for _ in range(size):
                read_seq, n_err = _mutate(consensus_seq, cfg.seq_error, rng)
                umi_obs, _ = _mutate(umi, cfg.umi_error, rng)
                rid = f"r{n_reads}"
                if fq2 is None:
                    full = f"{umi_obs}{cfg.adapter}{read_seq}"
                else:
                    # Independent sequencing of the two ends of the same fragment.
                    tail, _ = _mutate(consensus_seq[len(consensus_seq) - read_len:],
                                      cfg.seq_error, rng)
                    full = f"{umi_obs}{cfg.adapter}{read_seq[:read_len]}"
                    mate = _revcomp(tail)
                    fq2.write(f"@{rid} UMI:{umi_obs}\n{mate}\n+\n{qual_char * len(mate)}\n")
                fq.write(f"@{rid} UMI:{umi_obs}\n{full}\n+\n{qual_char * len(full)}\n")
                tr.write(f"{rid}\t{mol}\t{clone}\t{umi}\t{umi_obs}\t{n_err}\n")
                n_reads += 1

    if fq2 is not None:
        fq2.close()

    return {
        "reads": str(reads_path),
        "reads2": str(reads2_path) if cfg.paired else "",
        "truth_reads": str(truth_reads),
        "truth_molecules": str(truth_mols),
        "truth_consensus": str(truth_cons),
        "clones": str(truth_clones),
        "reference": str(reference) if variants else "",
        "truth_variants": str(truth_variants) if variants else "",
        "n_variants": len(variants),
        "n_reads": n_reads,
        "n_molecules": cfg.n_molecules,
        "n_distinct_umis": len(seen_umis),
        "n_umi_collisions": n_collisions,
        # Paste-ready barcode pattern for `migec checkout` on this corpus.
        "pattern": "N" * cfg.umi_len + cfg.adapter.lower(),
    }


def _variant_clones(cfg: SimConfig, rng: random.Random, reference: str):
    """One reference plus `n_clones - 1` single-substitution variants of it, each at `variant_af`.

    Positions are spread evenly and kept 10 bases clear of both ends, because an aligner soft-clips
    there and a variant nobody can align to is a measurement of the aligner rather than the caller.
    """
    n_var = cfg.n_clones - 1
    if n_var < 1:
        raise ValueError("variant_af needs n_clones >= 2: one reference and at least one variant")
    total = cfg.variant_af * n_var
    if not 0.0 < total < 1.0:
        raise ValueError(
            f"variant_af {cfg.variant_af} x {n_var} variants is {total:.3f} of the library; "
            f"it has to leave room for the reference"
        )
    margin = 10
    usable = cfg.seq_len - 2 * margin
    if usable < n_var:
        raise ValueError(f"{n_var} variants do not fit in {usable} usable positions")
    clones, variants = [reference], []
    for i in range(n_var):
        pos = margin + (i * usable) // n_var
        ref_base = reference[pos]
        alt = rng.choice([b for b in BASES if b != ref_base])
        clones.append(reference[:pos] + alt + reference[pos + 1:])
        variants.append((i + 1, pos, ref_base, alt))
    weights = [1.0 - total] + [cfg.variant_af] * n_var
    return clones, variants, weights


def _log_mu(mean: float, cv: float) -> float:
    """mu of a lognormal with the given arithmetic mean and sigma=cv."""
    import math

    return math.log(max(mean, 1e-9)) - 0.5 * cv * cv
