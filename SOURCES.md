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
| `assets/pipeline.svg` | `assets/pipeline.dot` | derived | `dot -Tsvg assets/pipeline.dot -o assets/pipeline.svg` |
| `assets/benchmark_threads.tsv` + `.json` | this machine, 2026-08-13 | **experimental** (a timing measurement) | `python scripts/benchmark_threads.py --reads 2000000 -o assets/` |
| `assets/SRR1763769.mig.tsv`, `assets/assemble.coverage.tsv` | `isalgo/umi_data` CI fixture | derived | `migec refine ci/SRR1763769_umi0.5pct.fq.gz -o r/ && migec assemble r/CTRL.fq.gz -o a/` |
| `assets/*.svg`, `assets/*.gp` | the tables beside them | derived | `migec plot assets/ -o assets/` |
| `assets/ctdna_titration.tsv` (+ `_runs.tsv`) | 100 runs of `PRJNA788522` + `PRJNA507366` | derived (migec over experimental data) | `python scripts/sra_fetch.py get <runs> -o simsen/` then `python scripts/ctdna_titration.py --reads simsen/ --out ctdna/ --design design.tsv` |
| `assets/ctdna_panel.bed` | inferred from coverage, 3 deepest `PRJNA788522` runs aligned to GRCh38 | **derived** | `sbatch scripts/ctdna_infer_panel.sbatch`: minimap2 -> `bedtools genomecov` above a depth floor -> `merge -d 50` -> named against Ensembl 110 |
| `assets/ctdna_per_target.tsv` | per-TARGET molecule counts, 72 runs | **derived** | `sbatch scripts/ctdna_per_site.sbatch` then `python scripts/ctdna_persite.py --molecules m.tsv --variants v.tsv --design d.tsv --out out/` |
| `assets/ctdna_minreads.tsv` | calls at `--min-reads` 1/3/5, 12 runs x 4 certified arms | **derived** | `sbatch scripts/ctdna_minreads.sbatch` |

Note: the three `scripts/ctdna_*.sbatch` are the cluster half of the ctDNA work and are committed
because a `Regenerate` cell that names a file nobody has is not provenance. They run on aldan3
against `ROOT=/projects/tcr_bcr_rnaseq/migec_ctdna` and build their own micromamba environment,
because the cluster's system python is 3.8 and migec needs >= 3.10; they install migec from PyPI
rather than from source, which doubles as a check that the published Linux wheels work.

Note: `ctdna_titration.tsv` is kept beside `ctdna_per_target.tsv` deliberately. The first divides a
library total by an amplicon count inferred from consensus prefixes; the second counts molecules
actually aligned to each target. The gap between them is what a reference buys, and deleting the
first would hide it: it reported 5 amplicons where there are 6 intervals (one of them off-target),
and a per-panel average where the weakest target holds **0.09-0.64x** the mean.

Result tables and figures are output, not data, so they live here next to the script that made
them rather than in `isalgo/umi_data`. Test corpora go the other way: HuggingFace, never git.

## Error-rate constants

The pre-amplification floor (`--pre-amp-error`, formerly `--rt-error`) is a property of the
protocol, so its values are cited rather than fitted. Nothing here was measured by us except the X2
row.

Never: only an RNA library has a reverse transcription step. On a DNA library the same floor is
supplied by library-preparation damage plus the first PCR cycle, and the two damage chemistries are
sourced below. The class names `rt` / `medium` / `high` are historical brackets on the *rate*, not
claims about mechanism.

| Value | What it is | Source | Provenance |
|---|---|---|---|
| 1e-4 per base | V(D)J RT reaction; single-UMI bases get Q40, bases covered by >=2 UMIs get Q60 | 10x Genomics Cell Ranger V(D)J algorithm documentation | vendor-stated |
| 1.54e-4 [1.36e-4, 1.74e-4] | the same floor, fitted from `SRR1763769` as the intercept of `e_out(c) = p_floor + a/c` | X2, `docs/quality_floor.rst`, `scripts/quality_floor.py` | **experimental**, ours; an upper bound (49.6% barcode occupancy inflates it) |
| 7.37e-5 (0.00737%) | TruSight Oncology 500 v2 error rate | Illumina product documentation | vendor-stated |
| Taq 4.3e-5 +/- 1.8; Pfu 2.8e-6; Phusion 2.6e-6; Pwo 2.4e-6, per bp per template duplication | polymerase fidelity, no RT | McInerney P, Adams P, Hadi MZ. *Error Rate Comparison during Polymerase Chain Reaction by DNA Polymerase.* Mol Biol Int 2014:287430. doi:10.1155/2014/287430, PMID 25197572 | published, third party |
| 0.3-6.6e-5 per base per cycle over nine polymerases; **linear-amplification errors 5 +/- 1x the per-cycle PCR rate** | why the FIRST cycle is the one that matters — it is copied into every read of the molecule | Shagin DA, Shagina IA, Zaretsky AR, Barsova EV, Kelmanson IV, Lukyanov S, Chudakov DM, Shugay M. *A high-throughput assay for quantitative measurement of PCR errors.* Sci Rep 2017;7:2718. doi:10.1038/s41598-017-02727-8, PMID 28578414 | published, ours |
| `C>A`/`G>T` transversions at low allele fraction, read-orientation biased | oxidation of guanine to 8-oxoG during **acoustic shearing**, in extracts carrying reactive contaminants. The pre-amplification floor of a DNA library, and the reason `--min-reads` cannot remove it: damage predates the barcode. Detected by orientation bias, not by family size | Costello M, Pugh TJ, Fennell TJ, Stewart C, Lichtenstein L, Meldrim JC, Fostel JL, Friedrich DC, Perrin D, Dionne D, Kim S, Gabriel SB, Lander ES, Fisher S, Getz G. *Discovery and characterization of artifactual mutations in deep coverage targeted capture sequencing data due to oxidative DNA damage during sample preparation.* Nucleic Acids Res 2013;41(6):e67. doi:10.1093/nar/gks1443, PMID 23303777 | published, third party |
| `C>T`/`G>A` transitions | cytosine (and 5-methylcytosine) deamination to uracil/thymine; the other DNA-library damage chemistry, dominant in FFPE material | Do H, Dobrovic A. *Sequence artifacts in DNA from formalin-fixed tissues: causes and strategies for minimization.* Clin Chem 2014;61(1):64-71. doi:10.1373/clinchem.2014.223040, PMID 25421801 | published, third party |
| 10,000 reads per barcode | coverage cap into the consensus (never into the count) | 10x Genomics: "Very high coverage (greater than 10,000 reads) of transcripts can be problematic because it degrades computational performance and adds little information." | vendor-stated |

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

### Ageing cohort, 5'-RACE bulk TCR beta (Britanova et al., J Immunol 2014, doi:10.4049/jimmunol.1302064)

| Item | Value |
|---|---|
| What | One HiSeq lane, ten donors multiplexed by a 4 nt sample tag, 149,588,907 read pairs, 16 nt UMI |
| Layout | Taken from the lane's own metadata, never reconstructed from a preset: the tag length, the adapter and the UMI shape differ between lanes of the same cohort |
| Provenance | experimental, published cohort; the raw reads are held internally and do not leave it |
| Derived | `assets/shallow_repertoire.tsv` (per-donor correction and consensus summary), written up in `docs/validation.rst` |
| Why | The real 1-3 reads/UMI library. Every shallow-regime claim was otherwise measured on simulated data or on a deep amplicon |

**Never: Experiment 1 raw reads must not leave the cluster.** Only derived summaries (histograms,
error-rate tables, consensus statistics) may be published or uploaded to HuggingFace.

Never: `scratch/spikein/S1_R2_2M.fq` on aldan3 is **corrupt past record 1,742,617**. Do not use it.

### MAGERI (Shugay et al., PLoS Comput Biol 2017, doi:10.1371/journal.pcbi.1005480, PMID 28475621)

| Item | Value |
|---|---|
| Error-model datasets | SRA `PRJNA352143` — UMI-tagged sequencing of a known template with 9 polymerases |
| Duplex sequencing | SRA `SRR1799908`; primer patterns `NNNNNNNNNNNNtgact` / `agtcaNNNNNNNNNNNN` |
| HIV protease amplicons | SRA `SRP052322`; patterns `NNNNNNNNNcagtttaacttttgggccatccattcc` / `ctatcggctcctgnnnn` |
| Companion repo | https://github.com/mikessh/mageri-paper (error model PDFs, analysis scripts) |
| Jar | `gh release download 1.1.1 --repo mikessh/mageri -p mageri.zip` (v1.1.1, 2016-12-20); runs on JDK 11 |
| Head to head | `python scripts/compare_mageri.py --out out/ --jar mageri.jar --molecules 20000 --clones 200 --min-count 1` -> `assets/mageri.tsv` |
| Provenance | experimental; `assets/mageri.tsv` is derived (both pipelines over one simulated library with known truth) |

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

**`sc5p_v2_hs_PBMC_1k` VDJ-T — the cell-barcode dataset** (fetched 2026-08-13):

| Item | Value |
|---|---|
| Fetch | `curl -O https://cf.10xgenomics.com/samples/cell-vdj/5.0.0/sc5p_v2_hs_PBMC_1k/sc5p_v2_hs_PBMC_1k_t_fastqs.tar` (569 MB, 1.1 GB unpacked) |
| Why this one | the **VDJ-T** library, not GEX: 569 MB against tens of GB, and it is the part that matters here. The 3' GEX libraries carry the same barcode structure at 100x the size |
| Layout | 5' v2: **R1 = 16 nt cell barcode + 10 nt UMI, exactly 26 nt and nothing else**; R2 = 90 nt cDNA; I1/I2 index reads present |
| Pattern | `^XXXXXXXXXXXXXXXXNNNNNNNNNN`, or `--preset 10x-v2`, or the slice `cell:0:16,16:26`. Purely positional — there is no constant sequence to anchor on, which is why the pattern grammar had to allow a scored-nothing pattern at a fixed offset |
| Measured | 3,155,166 reads, **100% assigned**, 221,024 barcodes at 14.28 reads each, effective UMI length 9.97/10. refine on **R2** (the mate carrying the cDNA): 305,702 molecules, **813 cells** called by OrdMag, clonality 0.0014 |
| Note: | refine and assemble take **R2**, not R1: trimming the pattern leaves R1 empty, so the payload evidence is vacuous there and the reported clonality comes out 1.0 saying so |
| Provenance | experimental (10x public) |
| Regenerate the fixture | `migec checkout R1 R2 --preset 10x-v2 --sample PBMC -o co/` then `migec subsample co/PBMC_R2.fq.gz -o sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz --keep 1` |

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
| Identity | **HXB2 2,328-2,595 on the minus strand**, 268 nt, 2 mismatches (0.7% divergence) — covers *pol*: the 3' end of **protease** and the start of **RT**. Placed by `scripts/diagnose.py` against `K03455.1`, fetched from NCBI |
| Note: headers | SRA-normalised to `@SRR1763769.N N/2`, so lane/tile/x/y are **gone**. Optical duplicates and index hopping are unmeasurable on this file — a property of the download, not a clean result |
| Provenance | experimental (ENA); the floor and its interval are derived |
| Note: | ENA's metadata gives every run in the study the same title, so it does **not** identify which runs are controls. This is HIV plasma — a quasispecies — so the estimator restricts to monomorphic positions rather than assuming clonality. The library is also 49.6% occupied on its 9 nt barcode, which `checkout` flags as saturated; the measured floor is an upper bound. |

### ctDNA UMI benchmark — Maruzani et al. 2024 (checked 2026-08-13)

| Item | Value |
|---|---|
| Paper | Maruzani R, Brierley L, Jorgensen A, Fowler A. *Benchmarking UMI-aware and standard variant callers for low frequency ctDNA variant detection.* **BMC Genomics** 25, 2024-09-03. [doi:10.1186/s12864-024-10737-w](https://doi.org/10.1186/s12864-024-10737-w), PMC11370058 |
| Synthetic template | `SRR10296599` — cfDNA from a healthy Han Chinese female, Roche ctDNA panel (17 genes), NextSeq 550, 6,230,802 spots, 2x96 nt, `PRJNA577992` |
| Tumour samples | `SRR15081468/70/72/77/80/82/93/94` — 8 pre-treatment metastatic breast cancer samples, COMET trial, custom 54-gene panel, HiSeq 2500, ~3.2 M spots, 2x~105 nt, `PRJNA745047` |
| Their benchmark | 303 COSMIC variants spiked at 0.005–0.075 VAF across 200x/450x/850x; six callers compared (Mutect2, bcftools, LoFreq, FreeBayes vs UMI-aware UMI-VarCal, UMIErrorCorrect) |
| Provenance | experimental (SRA), human clinical; their spiked truth is derived |

**Never: Neither accession carries a recoverable UMI, so neither is a migec input.** Established from the
data, not from the text:

- Both are **aligned BAM submissions** (`NCBI:align:db:alignment_sorted`), so the Illumina headers
  are stripped and the read lengths are already trimmed (101 vs 110 nt between mates).
- `migec suggest` finds **no pattern** in either mate of either run: near-uniform composition with
  no constant anchor after it, which is what diverse payload looks like rather than a barcode.
- `vdb-dump -T SEQUENCE` shows `CMP_LINKAGE_GROUP` **empty** — that is where an `MI`/`RX` tag would
  survive, so no UMI was preserved.
- The paper says so too: the UMIs of the synthetic dataset were **generated in silico**, at 9 nt
  with Phred fixed to 37, assigned by Poisson (λ=1/2) to reads *sharing start and end positions*.

Note: That assignment rule bakes in the co-terminal assumption **X1 measured as false**
(`docs/fragmented.rst`: 7.8% of 10x groups overall, 0.3% at ≥6 reads). cfDNA has preferred cut
sites so it is less wrong for a capture panel than for 3' GEX, but any comparison against their
numbers inherits it.

**What these are good for:** the *design* — spiking COSMIC variants at known VAF into a real cfDNA
background at controlled depth is the ground truth a consensus-quality claim wants — and two more
comparators for M5, `UMI-VarCal` and `UMIErrorCorrect`.

### SiMSen-Seq — the ctDNA data that DID keep its UMIs (fetched 2026-08-13)

Found by screening SRA rather than by assuming, after the Maruzani entry above established that
their two accessions had not. Both studies are the Gothenburg SiMSen-Seq protocol on **commercial
cell-free DNA reference material with certified mutant allele frequencies**, and both deposited
the reads untrimmed, so the **12 nt inline UMI is still there**. `migec suggest` recovers it from
base composition alone with no knowledge of the protocol — a 12 nt near-uniform run at cycles
0-11, then the constant 16 nt spacer (`ATGGGAAAGAGTGTCC`), which is exactly UMIErrorCorrect's
`-ul 12 -sl 16`.

| Item | Value |
|---|---|
| Pattern | `--bc-pattern '0:12'`, single-end, positional (no `--max-offset`) |
| Panel | multiplex amplicon, measured from the consensus prefixes: **5 amplicons** for `PRJNA788522`, **3** for `PRJNA507366` (whose own labels say `3plx`) |
| Cluster copy | `/projects/tcr_bcr_rnaseq/migec_ctdna/fastq` on aldan3, 100 runs, 6.9 GB (2026-08-13). A convenience cache — the fetch command above regenerates it |
| Provenance | experimental (SRA), commercial reference material; the certified VAFs are the vendor's |
| Fetch | `python scripts/sra_fetch.py get <runs> -o data/` (NCBI S3, ~7 MB/s) |
| Analyse | `python scripts/ctdna_titration.py --reads data/ --out out/ --design design.tsv` |
| Designs | `curl -sG 'https://www.ebi.ac.uk/ena/portal/api/search' -d result=read_run -d query=study_accession=<PRJ> -d fields=run_accession,sample_alias -d format=tsv -d limit=0` |

**`PRJNA788522` — the titration.** Österlund T, Filges S, Johansson G, Ståhlberg A.
*UMIErrorCorrect and UMIAnalyzer: Software for Consensus Read Generation, Error Correction, and
Visualization Using Unique Molecular Identifiers.* **Clin Chem** 2022;68(11):1425–1435.
[doi:10.1093/clinchem/hvac136](https://doi.org/10.1093/clinchem/hvac136), PMID 36031761.

72 runs, 41.8 M reads, Illumina MiniSeq, single-end 151 nt. The design is in `sample_alias` as
`<input>ng_<vaf>_<depth>x_rep_<n>`: mutant allele frequency **0% (`WT`), 0.125%, 0.25%, 1%** plus
an undiluted `cell_line` arm, crossed with DNA input **5 / 20 / 80 ng** and depth **3.3 / 10 / 30x
reads per UMI**, three replicates each. Note: the `WT` arm is a **true negative** — its variant
frequency is zero by construction, so anything a caller reports there is its own false-positive
rate on real chemistry rather than on a simulation.

**`PRJNA507366` — the polymerase panel.** Filges S, Yamada E, Ståhlberg A, Godfrey TE. *Impact of
Polymerase Fidelity on Background Error Rates in Next-Generation Sequencing with Unique Molecular
Identifiers/Barcodes.* **Sci Rep** 2019;9(1):3503.
[doi:10.1038/s41598-019-39762-6](https://doi.org/10.1038/s41598-019-39762-6), PMID 30837525.

28 runs, 81.4 M reads. Note: the design is **not** in `sample_alias` (every row reads
`SeraCare_Reference_Material`) — it is in **`library_name`** (the enzyme: Phusion, Accuprime,
Accuprime_hifi, Platinum, Platinum_Hifi, Platinum superfi, three replicates each) and in
`experiment_title` (`Wildtype` against `0.031% VAF` / `0.0625% VAF` / `0.125% VAF`). Reading the
alias alone would have called this study undesigned. It extends the frequency range **4x below**
PRJNA788522's floor, and it is the only public dataset here that varies the polymerase while
holding template and protocol fixed — which is the comparison `--rt-error auto` needs and the one
McInerney 2014's published fidelities can be checked against.

### MRD and duplex — the read-structure screen of 2026-08-14

Two studies were proposed as the MRD and exome arms. Neither is usable **as deposited**, and both
verdicts are from the data rather than from the methods section. Screening for a public substitute
then found two TwinStrand Duplex Sequencing studies that did keep their tags. Every accession
touched is listed, negatives included, so that the next screen does not repeat this one.

| study | bioproject | runs | assay | UMI as deposited | verdict |
|---|---|---|---|---|---|
| `SRP475624` | `PRJNA1028782` | 234 | MAESTRO-Pool MRD, cfDNA + WGS, duplex UMI adapters | not checkable | controlled access, dbGaP `phs003447` |
| `SRP578416` | `PRJNA1250480` | 8 | xGen AML capture panel, UMI adapters | absent | UMIs existed, SRA did not keep them |
| `SRP598252` | `PRJNA1278818` | 82 | duplex capture, normal bladder | not checkable | controlled access, consent `GRU` |
| `SRP613942` | `PRJNA1312025` | 51 | TwinStrand duplex capture, normal human skin | present, inline, both mates | **usable** |
| `SRP677910` | `PRJNA1425166` | 12 | duplex, HEK293 LIG1 mutagenesis | present, inline, both mates | **usable** |
| `SRP579806` | `PRJNA1252794` | 36 | duplex capture, retroviral and lentiviral vectors | present, inline, both mates, non-uniform | usable, second choice |

| Item | Value |
|---|---|
| Provenance | experimental (SRA); every verdict below is derived |
| Probe | `python scripts/sra_fetch.py probe <runs>` |
| Peek | `curl -r 0-4000000 <ena fastq url>`, `gunzip -c`, then `migec suggest <peek>.fq --cycles 30` on **each mate separately** |

**`SRP475624` — MAESTRO-Pool, controlled.** SRA study title: *MAESTRO-Pool enables highly parallel
and specific mutation-enrichment sequencing for minimal residual disease detection in cohort
studies* (submitting centre `NCI-PHS003447`; the paper is not cited here because it has not been
retrieved). 9 melanoma patients, 98 plasma samples screened for 22,333 mutations, Kapa Hyper
Prep with custom dual-index **duplex** UMI adapters from IDT, NovaSeq 6000, 215 Targeted-Capture
runs plus 19 WGS, 26,613,420,959 spots, 959 GB. ENA answers `Protected file(s). Go to dbGap` for
every run; the consent group is `DS-MSC-MDS`. It is the ideal MRD arm and it needs an approved
dbGaP application before any of it can be fetched. Recorded here so the next screen does not
rediscover it.

**`SRP578416` — the UMIs existed and SRA dropped them.** Targeted DNA-Seq of human AML treated with
TP-0903 (Ohio State), xGen DNA Lib Prep MC UNI kit, xGen AML Cancer Hybridization Panel plus a
three-gene spike-in, NovaSeq 6000, 8 runs, 65,630,386 spots, 8.8 GB. The design description says
"Picard Tools was used to perform UMI-consensus calling", so the library carried a UMI. It is not
in the public record:

- The submissions are **aligned BAMs** (`NCBI:align:db:alignment_sorted`, `bam-load.3.0.8`), and
  `vdb-dump -T SEQUENCE -C LINKAGE_GROUP` is **empty** — that is where an `RX` would have survived.
- Read names are rewritten to ordinals (`NAME` = `1`, `2`, `3`), so nothing is hiding in the header.
- Spots hold **two** biological reads and no technical read; per-cycle composition has no UMI run in
  either mate, and the mates of a pair overlap as plain genomic sequence (one fragment: 161 bp, R1
  `[0,150)`, R2 `[11,161)`).
- The original `U_21_0148.bam` is listed on `s3://sra-pub-src-18` with `access_type="Use Cloud Data
  Delivery"` and `free_egress="-"`. Anonymous `GetObject` returns 403; `--request-payer requester`
  returns `AccessDenied` on `GetObject`, `HeadObject` and `ListObjectsV2`; and
  `prefetch --type all` resolves only the 1,339 MB normalised run. **NCBI Cloud Data Delivery**, an
  account-bound request, is the only remaining route.

This is the same failure mode as the Maruzani entry above, and it is the case that motivated BAM
input: had the pre-consensus BAM been deposited rather than an SRA-normalised run, it would have
been a one-command migec input.

**`SRP613942` — Mutomics: mutations in normal human skin.** 51 runs, 7,048,584,922 read pairs,
988.2 GB, NovaSeq 6000, 2x151. "Library prepped via targeted capture using custom probe panel from
TwinStrand Biosciences and Duplex Sequencing kit". Sun-shielded against sun-exposed skin, blood and
saliva, with a UV dose series (300 / 1000 J/m2).

**`SRP677910` — HEK293 LIG1 WT against K845N.** 12 runs, 3,768,336,083 read pairs, 501.7 GB,
NovaSeq X Plus, 2x151. Duplex sequencing of mutation frequency after menadione at 0 / 8 / 20 uM,
48 h recovery, two replicates per cell. The smaller of the two and the one with a dose-response
ground truth, so it is the first target.

Both carry the same recoverable layout, read off the data with no knowledge of the kit: an **8 nt
inline tag followed by a fixed `T`, on both mates** — `migec suggest` calls cycles 0-7 UMI and
cycle 8 variable, and cycle 8 is `T` in 77.8% of `SRP613942` reads and 80.4% of `SRP677910` reads
against 4.4-8.5% for the other three bases. In migec's grammar that is `^NNNNNNNN.` as master and
the same as slave, which is the `duplex` preset's shape at this kit's lengths. Note: cycle 0 is
semi-degenerate by design (T at 0.09 against C at 0.37), which is a property of the TwinStrand
adapter and not a defect — `checkout` measures the effective barcode length rather than assuming
`4^L`.

Never: **these give a single-strand consensus, not a duplex one.** The two strands of a molecule
carry the tag pair swapped, and migec concatenates master and slave into one 16 nt identifier, so
the two strands are two different UMIs and are consensed separately. That is the correct behaviour
for this codebase (`docs/roadmap.rst`, "deliberately not doing"), and it is also what makes these
studies a *measurement*: TwinStrand's own duplex calls are the comparator for what single-strand
consensus does and does not remove.

**`SRP579806`** (quantifying the mutational landscape of retroviral and lentiviral vectors in gene
therapy patients, Duplex Sequencing library prep and ADA hybridization probe capture, NovaSeq X
Plus, 2x150, 36 runs) also carries an inline tag on both mates, but the prefix distribution is
strongly non-uniform (`TTGCGAAG` alone takes 3 of 8 sampled reads) and cycles 5-7 are only partly
constant, so the tag boundary is not clean from composition alone. Second choice.

**`SRP598252`** (somatic mutation in normal bladder, duplex sequencing over a bladder-cancer gene
panel, 82 runs) is controlled under consent `GRU`: ENA answers `Protected file(s). Go to dbGap`,
and `fastq-dump` exits 3.

### Calib (github.com/vpc-ccg/calib)

| Item | Value |
|---|---|
| Use | comparator for UMI grouping accuracy — it clusters on barcode *and* sequence, we (today) on barcode alone |
| Get it | `git clone https://github.com/vpc-ccg/calib && cd calib && make` → `calib`, `calib_cons` |
| Run | `calib -f R1.fq -r R2.fq -l <barcode_len> -o prefix` → `prefix.cluster` |
| `.cluster` format | 9 TSV columns: `cluster_id, node_id, read_id, f_name, f_seq, f_qual, r_name, r_seq, r_qual` (verified against the upstream README, 2026-08-13) |
| Truth used here | **our** simulator, `tests/synthetic/_sim.py`, which writes `truth_reads.tsv` (`read_id`, `molecule_id`). Calib's own simulator emits no read→molecule map |
| Compared by | `scripts/compare_calib.py` — adjusted Rand index, plus split and merge fractions separately |
| Storage | Never: do not store simulated reads — record the exact command and seed here instead |
| Provenance | derived (simulated) |

### UMI-tools (github.com/CGATOxford/UMI-tools)

| Item | Value |
|---|---|
| Use | comparator for UMI grouping accuracy — the map-first approach: align the raw reads, group on *(position, UMI)* |
| Cite | Smith, Heger & Sudbery, *Genome Res* 27:491 (2017), doi:10.1101/gr.209601.116 |
| Get it | `uv pip install --no-build-isolation umi_tools` (v1.1.6 here). Never: it needs `setuptools` present and `distutils`, so on Python 3.12+ install `setuptools` first and pass `--no-build-isolation` |
| Run | `umi_tools group -I sorted.bam --group-out groups.tsv --umi-separator _` after moving the barcode into the read name |
| Truth used here | **our** simulator, `tests/synthetic/_sim.py` — `truth_reads.tsv` plus `clones.fa` as the reference to align to |
| Compared by | `scripts/compare_grouping.py`, scored by `scripts/compare_calib.py`'s adjusted Rand index |
| Provenance | derived (simulated) |

### fgbio (fulcrumgenomics.github.io/fgbio)

| Item | Value |
|---|---|
| Use | the other map-first comparator; `GroupReadsByUmi -s adjacency` is the standard for consensus calling from duplex and single-strand UMIs |
| Get it | `brew install fgbio` (v4.1.1 here). Never: it needs a **JDK 17 or newer** — a JDK 11 default gives `UnsupportedClassVersionError` from htsjdk, not a version message. `brew install openjdk@21` and point `JAVA_HOME` at it |
| Run | `fgbio FastqToBam --read-structures <n>M+T`, align carrying `RX`, then `fgbio GroupReadsByUmi -s adjacency -e 1` → `MI` tag per read |
| Truth used here | as UMI-tools above |
| Compared by | `scripts/compare_grouping.py` |
| Provenance | derived (simulated) |

### `assets/grouping_tools.tsv`

| Item | Value |
|---|---|
| What | grouping accuracy, wall clock and peak RSS for migec, UMI-tools and fgbio over a clone-diversity and a depth sweep |
| Regenerate | `python scripts/compare_grouping.py --out DIR --molecules 20000 --clones {1,20,200,20000} --coverage {1.2,2.5,5,10} --umi-error 3e-3 --tsv out.tsv` |
| Provenance | derived (computed from simulated reads); the drawn conclusion is in `docs/grouping.rst` |

### MIGEC 1.2.9 (the Groovy implementation this repo replaced)

| Item | Value |
|---|---|
| Use | the implementation comparison: same barcode dialect, same sheet, same library, both pipelines end to end |
| Get it | `gh release download 1.2.9 --repo antigenomics/migec -p 'migec-1.2.9.zip'` → `migec-1.2.9.jar`. Runs on JDK 11 |
| Run | `java -jar migec-1.2.9.jar Checkout -cute sheet.txt R1.fq.gz . co/` then `Assemble -c --filter-collisions -m <n> co/S1_R0.fastq.gz . asm/` |
| Never | `--min-count` must be given to BOTH -- v1 defaults to 5 and migec 2 to 1, and v1 names its output `.t5.` because of it. Leaving each at its own default compares defaults |
| Source | this repository, branch `legacy-v1` / tag `v1-final`; the release predates the rewrite |
| Compared by | `scripts/compare_migec_v1.py` → `assets/migec_v1.tsv`, drawn in `docs/validation.rst` |
| Provenance | derived (simulated reads); the jar itself is the published 2018 artefact |

### `assets/migec_v1.tsv`

| Item | Value |
|---|---|
| What | consensuses, exactness, wall clock and peak RSS for MIGEC 1.2.9 and migec 2 at matched MIG-size thresholds |
| Regenerate | `python scripts/compare_migec_v1.py --out DIR --jar migec-1.2.9.jar --molecules 20000 --clones 200 --coverage 8 --min-count {1,5} --tsv out.tsv` |
| Provenance | derived (computed) |

### `assets/mageri.tsv`

| Item | Value |
|---|---|
| What | consensuses, exactness, wall clock and peak RSS for MAGERI 1.1.1 and migec 2 at matched MIG-size thresholds, one row per (tool, min_count) |
| Regenerate | `for mc in 1 2 5; do python scripts/compare_mageri.py --out DIR$mc --jar mageri.jar --molecules 20000 --clones 200 --min-count $mc --tsv out_$mc.tsv; done`, then concatenate |
| Never | MAGERI's threshold is `forceOverseq`/`defaultOverseq` in its preset XML, not a flag. The script rewrites the exported preset and then reads back the threshold MAGERI REPORTS it used, refusing to score if they differ |
| Note | the `migec-2+minimap2` row folds `minimap2 -ax sr -y` onto the same reference into migec's clock, because MAGERI aligns as part of its run. It costs 0.05 s of 0.35 |
| Provenance | derived (computed) |

### `assets/collision_split.tsv`

| Item | Value |
|---|---|
| What | how often `assemble`'s linkage sub-clustering separates two molecules that truly collided on a barcode, against depth |
| Regenerate | `python scripts/collision_split.py --out DIR --coverage 5 20 40 80 160 --tsv out.tsv` |
| Never | "collided" is the TRUE barcode, never the observed one -- an observed coincidence caused by a barcode error is what `refine` corrects, and counting it makes the collision rate grow with the read count |
| Provenance | derived (computed from simulated reads); drawn in `docs/grouping.rst` |

### Example QC figures in `assets/`

| Item | Value |
|---|---|
| Use | the panels the README embeds: barcode rank, MIG size spectrum, rank/Zipf, consensus quality, and the two barcode-error-against-depth panels |
| Regenerate | `python scripts/example_figures.py` -- simulates, runs all three stages, copies the tables, redraws |
| Tables | `assets/PBMC.cell_rank.tsv`, `PBMC.sizes.tsv`, `PBMC.umi_errors.tsv`, `assemble.quality_by_depth.tsv`, `assemble.coverage.tsv` |
| Shape | 120 cells x 40-200 molecules, Pareto(1.1) depth, plus 4,000 ambient barcodes with 1-3 molecules each; 0.2% barcode and 1% payload error, seed 3 |
| Why that shape | a uniform simulator draws every panel as a straight line. The Pareto tail reaching ~19,000 reads on one molecule is also what makes the barcode-error panels show anything: saturation of the `3L` neighbour shell is only visible on a molecule deep enough to have spawned all of it. The ambient population is what makes the knee; the Pareto tail is what gives the rank curve slope; the per-base error is what gives emitted quality a spread to draw quartiles over |
| Provenance | **derived (simulated)** -- never presented as a measurement. `assets/SRR1763769.mig.tsv` beside it is experimental; see the Primer ID entry |
| Storage | the TSVs are committed, the FASTQ is not (`--keep` writes it locally) |

### minibwa (github.com/lh3/minibwa)

| Item | Value |
|---|---|
| Use | downstream aligner in the `docs/downstream.rst` contract table, alongside minimap2 and bwa |
| Get it | `git clone https://github.com/lh3/minibwa && cd minibwa && make` → `./minibwa` |
| Version run | `0.7-r424-dirty`, commit `f0e1174` (2026-08-10), built on darwin/arm64 |
| Run | `minibwa index ref.fa` then `minibwa map -y ref.fa cons.fq.gz` |
| Result | 600/600 consensus records kept `RX`, `CB`, `MI`, `BC`, `cD`; `samtools sort` → `quickcheck` valid |
| Note | the comment flag is **`-y` on `map`** (the minimap2 spelling) and **`-C` on the legacy `mem` subcommand** (bwa's). Each rejects the other's flag with a non-zero exit, so the tags are never dropped silently |
| Checked by | `tests/unit/test_downstream.py::test_an_aligner_carries_the_tags_into_the_sam[minibwa]`, skipped when not on `PATH` |
| Provenance | derived (run here against our own synthetic 10x-shaped fixture) |

### Reference genome and panel definitions (aldan3, 2026-08-14)

Not fetched here — already on the cluster, and the panel is **inferred from the data** rather than
taken from a vendor file.

| Item | Value |
|---|---|
| Genome | `/projects/cdr3_common/reference/genome/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa` (3.15 GB) + `.fai` |
| Annotation | `Homo_sapiens.GRCh38.110.chr.gtf` (1.46 GB), Ensembl 110 |
| Note: contig naming | **Ensembl style — `>1`, not `>chr1`.** Every BED, region string and tool must match, or it silently intersects nothing |
| Note: existing indices | `bwamem2_index/` is **empty**, and the `.bwameth.c2t.*` index is bisulfite-converted and unusable for ordinary alignment. Build a minimap2 index |
| Panel | derived: align the consensus, `bedtools genomecov` above a depth floor, `merge -d 50`, then name the intervals against the GTF. `scripts/ctdna_infer_panel.sbatch` |
| Provenance | experimental (Ensembl reference); the inferred panel is **derived** |

**Why inferred rather than vendor-supplied:** the amplicon count in `assets/ctdna_titration.tsv`
came from a consensus-prefix tally, which assumes the panel is evenly covered. It is not — observed
shares 20.4 / 17.4 / 16.9 / 15.9 / 7.6%, so `molecules / n_amplicons` overstates the weakest target
by **2.6x**. Real coordinates turn a per-target average into a per-target count.

Two production pipelines read as design references for the panel handling, neither run here:

| Repo | What it is good for |
|---|---|
| [AWGL/TSO500_post_processing](https://github.com/AWGL/TSO500_post_processing) | Illumina TSO500 post-processing. Ships `hotspot_variants/*.bed`, `hotspot_coverage/*combined.bed`, `vendorCaptureBed_100pad_updated.bed` and `TSO_extra_padding_chr.interval_list` — note the padding beyond Illumina's +/-2 bp, and the `chr`-prefixed naming, which is the opposite convention to the Ensembl reference above |
| [ikmb/exome-seq](https://github.com/ikmb/exome-seq) | Exome capture. `--kit xGen_v2 \| Agilent_v7 \| xGen_pan_cancer`; `--assembly GRCh38_no_alt` recommended for short reads; `--baits`/`--targets` as Picard interval lists, `--panel cardio\|cancer\|liver\|...`, and **`--amplicon_bed` for amplicon primer positions** — the same object this job infers |
| [AstraZeneca-NGS/reference_data](https://github.com/AstraZeneca-NGS/reference_data/tree/master/hg38) | **Ready-made hg38 capture BEDs**, no vendor login: `Exome-Agilent_V2/V4/V5/V6` (plus `_UTR` and `_Padded` variants), `Exome-IDT_V1.bed`, `Exome-MedExome.bed`, `Exome-NGv3.bed`, `Exome-AZ_V2.bed`, `Exome-Agilent-OneSeq.bed`, and `CDS-canonical.bed` (5.7 MB). Also a `tricky_regions/` directory. Fetch one file with `curl -sL https://raw.githubusercontent.com/AstraZeneca-NGS/reference_data/master/hg38/bed/<name>.bed` |

Note: these are **exome capture** targets, so they are the right object for
`notebooks/exome_capture.py` and the wrong one for the SiMSen-Seq amplicon panels above — an
amplicon panel's coordinates are set by its primers, not by a capture bait set, which is why they
are inferred here. `CDS-canonical.bed` is still useful against the inferred panel, as a check on
whether the amplicons land in coding sequence.

### UMIErrorCorrect / UMIAnalyzer — Österlund et al. (Clin Chem 2022, doi:10.1093/clinchem/hvac136)

Read 2026-08-13. Not a data source: a **comparator and a design reference** for the map-first order.

| Item | Value |
|---|---|
| What it is | a Python pipeline + R Shiny app for UMI ctDNA panels: `preprocess.py` → `run_mapping.py` (bwa mem) → `umi_error_correct.py` → `call_variants.py` |
| Order | aligns **raw reads first**, then groups on *(chromosomal position, UMI)* with edit distance ≤ 1 — the opposite order to migec, and the reason `docs/downstream.rst` has a section weighing the two |
| Layout grammar | `-ul` UMI length, `-sl` spacer length (SiMSen-Seq: 12 and 16). Our positional slices and `--read-structure` cover the same ground |
| Depth reported | **3.3 and 10 reads per UMI** on SiMSen-Seq ctDNA — independent corroboration that 1–3 reads/UMI is the ordinary regime, not the exotic one |
| Consensus cutoff | group size ≥ 3 by default; we default `--min-reads 1` and report the cutoff sweep instead |
| Variant caller | beta-binomial background per position, Q ≥ 20 cutoff. **Not implemented here** — it needs a reference and an alignment, which is a variant caller and not one of migec's five commands |
| Public data cited by it | PRJNA788522, PRJNA507366 (SiMSen-Seq); PRJNA577992 (Roche Avenio, QiaSeq); PRJEB31811 (Archer) |
| Provenance | reference (published method), not data |

### UMI RNA-seq — Fennell et al. / NCGR (Sci Rep 2018, doi:10.1038/s41598-018-31064-7)

The layout behind the `smarter-umi` preset, and a worked example of why a deposited FASTQ has to be
checked rather than assumed.

| | |
|---|---|
| Pipeline | [`ncgr/UMI-analysis`](https://github.com/ncgr/UMI-analysis) — Perl/C, single-end plant scRNA-seq |
| Layout | **10 nt inline UMI at the read start**, then the `GGG`/`GGGG` the Clontech SMARTer template switch leaves. Read off `fastq_qual_filter in.fq good.fq bad.fq log 30 0 10` (offset 0, length 10) and `fastq_umi_clipper`, which moves those 10 bases into the read header |
| Accession | SRP150352, 20 single-end HiSeq 2000 runs, ENA `filereport?accession=SRP150352&result=read_run` |
| Fetch | `curl -O https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR729/00<n>/<run>/<run>.fastq.gz` |
| Provenance | experimental |

**Never: the SRA copy cannot be reprocessed — the UMI is gone.** `fastq_umi_clipper` writes the UMI
into the FASTQ *header*, and SRA rewrites headers to `@<run>.<n> <n>/1`, so it is discarded. The
reads themselves are 36 nt with a flat base composition from cycle 0 (checked on SRR7295905,
SRR7295928 and SRR7295906, 2026-08-13): no G enrichment at cycles 11-13, where the template-switch
`GGG` would be if the 10 nt UMI were still in front of it. `migec suggest` says so unprompted —
*"the only near-uniform run sits after the last constant sequence, with nothing to anchor it. That
is what diverse payload looks like"* — which is the intended behaviour and is why the preset is
documented from the pipeline's own source rather than fitted to this data.

## HuggingFace — `isalgo/umi_data`

A git + git-lfs mirror at `~/hf/umi_data`, written by committing and pushing **in the mirror**,
one commit per change set. Never through the HTTP API — that writes remotely only, leaves the
mirror silently stale, and lands one commit per call. Note: The repo is **public**.

```
umi_data/
  ci/            small slices for CI, subsampled by WHOLE barcodes (all reads of N barcodes)
  README.md
  SOURCES.md     a copy of this file, plus a section on what is and is not shipped
```

Never: **sequences and metadata only** — `.txt`, `.md`, `.tsv.gz`, `.json`, `.fastq.gz`, `.fa.gz`,
`.sam`, `.bam`. No reports, figures, logs or pipeline output; a derived results table is *output*,
not data, even when it is a TSV, and lives in this repo next to the script that made it. A
`results/` directory was published there once and has been removed.

Still to add: `truth/` (spike-in clonotypes) and `whitelists/` (10x barcode lists with their
upstream and license), neither of which has been fetched yet.

**Published at [huggingface.co/datasets/isalgo/umi_data](https://huggingface.co/datasets/isalgo/umi_data)**
(2026-08-13): `ci/SRR1763769_umi0.5pct.fq.gz` and `ci/sc5p_v2_hs_PBMC_1k_t_cells1pct.fq.gz`, each
all the reads of a fraction of the *barcodes*, built with `migec subsample`. Written through the
local git+git-lfs mirror at `~/hf/umi_data`, never the HTTP API.

Never: Not in this dataset: aldan3 Experiment 1 raw reads, and anything regenerable by a one-line
command (record the command here instead of storing gigabytes in LFS).
