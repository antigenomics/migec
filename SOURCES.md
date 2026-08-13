# SOURCES — antigenomics/migec

Every dataset this repo ships, consumes or benchmarks against, where it came from, and how to
regenerate it.

**Experimental** = measured/sequenced. **Derived** = computed by us from something else. The two
are never conflated in a table row.

## Shipped in this repo

| Artifact | Origin | Provenance | Regenerate |
|---|---|---|---|
| `tests/synthetic/_sim.py` output | none — generated | derived | `SimConfig(seed=...)`; every truth file is a pure function of the seed |
| `tests/cpp/doctest.h` | [doctest](https://github.com/doctest/doctest) 2.4.11, MIT | vendored | copy from upstream release |

Nothing else is committed. Test corpora live on HuggingFace (below) rather than in git.

## Benchmark data

Nothing here has been fetched yet — this table is the plan of record, and rows gain a "fetched"
date as they land.

### MIGEC (Shugay et al., Nat Methods 2014, doi:10.1038/nmeth.2960)

| Item | Value |
|---|---|
| Public spike-ins | BioProject `PRJNA239303`, runs `SRR1200517`–`SRR1200520`, AMPLICON paired |
| Ground truth | Supplementary Table 1a: 12 control clonotypes (5 TRA, 5 TRB, 2 IGH) with expected frequencies |
| Truth source | `41592_2014_BFnmeth2960_MOESM376_ESM.pdf`, extractable with `pdftotext -layout` |
| Provenance | experimental |
| Chemistry | MiSeq 2×150; SMART adapter `GTGGTATCAACGCAGAG` |

| Item | Value |
|---|---|
| Experiment 1 (internal) | `/projects/cdr3_ngs/2012/12_alvaro_ab_bcr_nnn/` (IGH) and `/projects/cdr3_ngs/2012/08_alvaro_nnnb/` (TCR) on aldan3 |
| Staged copies | `/projects/tcr_bcr_rnaseq/data_migec_exp1/{IGH_P41,TCR_Project25}_R{1,2}.fastq.gz` |
| Access | `aldan3 ls`, `aldan3 pull` (see `~/vcs/code/aldan3-client`) |
| Provenance | experimental, unpublished |

⛔ **Experiment 1 raw reads must not leave the cluster.** Only derived summaries (histograms,
error-rate tables, consensus statistics) may be published or uploaded to HuggingFace.

⛔ `scratch/spikein/S1_R2_2M.fq` on aldan3 is **corrupt past record 1,742,617**. Do not use it.

### MAGERI (Shugay et al., PLoS Comput Biol 2017, doi:10.1371/journal.pcbi.1005480, PMID 28475621)

| Item | Value |
|---|---|
| Error-model datasets | SRA `PRJNA352143` — UMI-tagged sequencing of a known template with 9 polymerases |
| Duplex sequencing | SRA `SRR1799908`; primer patterns `NNNNNNNNNNNNtgact` / `agtcaNNNNNNNNNNNN` |
| HIV protease amplicons | SRA `SRP052322`; patterns `NNNNNNNNNcagtttaacttttgggccatccattcc` / `ctatcggctcctgnnnn` |
| Companion repo | https://github.com/mikessh/mageri-paper (error model PDFs, analysis scripts) |
| Provenance | experimental |

The patterns above are quoted verbatim from the paper's Methods and are directly reusable as
`checkout` test cases.

Reference values worth keeping (MAGERI Methods, for comparison rather than reimplementation):
UMIs below Phred 20 discarded; MIG pairs differing by 1 or 2 substitutions with size ratios above
20× and 400× treated as error children; MIG size threshold at the square root of the distribution
peak; consensus core 30 bases with ±5 offset; `CQS = (40/3)·(4f − 1)`; per-substitution error
rates fitted as Beta, counts as Beta-Binomial, `Q = −10 log10 P`, capped at 100.

### 10x Genomics

| Item | Value |
|---|---|
| Datasets | four "Connect-generated GEX+VDJ" sets: human PBMC, mouse PBMC, human melanoma, mouse splenocytes (Cell Ranger 6.0.1) |
| Whitelists | barcode whitelists shipped with Cell Ranger (`10XGenomics/cellranger`) |
| Reference calls | the published `filtered_feature_bc_matrix` — used as a comparator, not re-run |
| Provenance | experimental (10x), reference calls derived |

**X1 (read-start dispersion) used `pbmc_1k_v3`**, Cell Ranger 3.0.0, GRCh38-3.0.0:

| Item | Value |
|---|---|
| BAM | `https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_possorted_genome_bam.bam` (4.79 GB) + `.bai` (4.6 MB) |
| Fetch | not fetched — the server serves HTTP range requests (verified 206, 2026-08-13), so `pysam.AlignmentFile(url)` reads only the regions asked for |
| Regions | `11:65497688-65508073` (MALAT1), `7:5527151-5530601` (ACTB), `15:44711477-44718877` (B2M) — chosen for expression, since an unexpressed locus has no UMI with more than one read |
| Contig names | plain `1`, `2`, … `X` — **not** `chr1`. Cell Ranger's GRCh38-3.0.0 reference is Ensembl-styled |
| Regenerate | `python scripts/read_start_dispersion.py --bam <url> --region … ` |
| Provenance | experimental (10x); the dispersion statistics are derived |

### HIV-1 Primer ID — the X2 quality-floor control

| Item | Value |
|---|---|
| Run | `SRR1763769`, 2,122,456 read pairs, study `PRJNA272736` |
| Fetch | `curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR176/009/SRR1763769/SRR1763769_2.fastq.gz` (248 MB; R1 is not needed — the Primer ID is on the cDNA primer in **R2**) |
| Layout | 9 nt Primer ID, then `CAGTTTAACTTTTGGGCCAT`; recovered from the data by per-cycle entropy, not from the protocol |
| Paper | Zhou, Jones, Mieczkowski & Swanstrom, *J Virol* 89(16):8540–8555, 2015, [doi:10.1128/JVI.00522-15](https://doi.org/10.1128/JVI.00522-15) — reports a residual error rate of ~1 in 10,000 |
| Regenerate | `python scripts/quality_floor.py --reads SRR1763769_2.fastq.gz --out x2/ --window 180` (X2, the quality floor) and `python scripts/permutation_nulls.py --reads SRR1763769_2.fastq.gz --out x3/ --cycles 32 --window 180` (X3, the three permutation nulls) |
| Provenance | experimental (ENA); the floor and its interval are derived |
| ⚠ | ENA's metadata gives every run in the study the same title, so it does **not** identify which runs are controls. This is HIV plasma — a quasispecies — so the estimator restricts to monomorphic positions rather than assuming clonality. The library is also 49.6% occupied on its 9 nt barcode, which `checkout` flags as saturated; the measured floor is an upper bound. |

### Calib (github.com/vpc-ccg/calib)

| Item | Value |
|---|---|
| Use | comparator for UMI grouping accuracy — it clusters on barcode *and* sequence, we (today) on barcode alone |
| Get it | `git clone https://github.com/vpc-ccg/calib && cd calib && make` → `calib`, `calib_cons` |
| Run | `calib -f R1.fq -r R2.fq -l <barcode_len> -o prefix` → `prefix.cluster` |
| `.cluster` format | 9 TSV columns: `cluster_id, node_id, read_id, f_name, f_seq, f_qual, r_name, r_seq, r_qual` (verified against the upstream README, 2026-08-13) |
| Truth used here | **our** simulator, `tests/synthetic/_sim.py`, which writes `truth_reads.tsv` (`read_id`, `molecule_id`). Calib's own simulator emits no read→molecule map |
| Compared by | `scripts/compare_calib.py` — adjusted Rand index, plus split and merge fractions separately |
| Storage | ⛔ do not store simulated reads — record the exact command and seed here instead |
| Provenance | derived (simulated) |

## HuggingFace — `isalgo/umi_data`

Does not exist yet. When created: a git + git-lfs mirror at `~/hf/umi_data`, written by
committing and pushing **in the mirror**, one commit per change set. Never through the HTTP API —
that writes remotely only, leaves the mirror silently stale, and lands one commit per call.

Planned layout:

```
umi_data/
  ci/            small slices for CI, subsampled by WHOLE UMIs (all reads of N UMIs)
  truth/         ground-truth tables (spike-in clonotypes, simulated molecule tables)
  whitelists/    barcode whitelists + a README recording their upstream and license
  SOURCES.md
```

⛔ Not in this dataset: aldan3 Experiment 1 raw reads, and anything regenerable by a one-line
command (record the command here instead of storing gigabytes in LFS).
