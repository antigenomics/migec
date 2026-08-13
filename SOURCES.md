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

Result tables and figures are output, not data, so they live here next to the script that made
them rather than in `isalgo/umi_data`. Test corpora go the other way: HuggingFace, never git.

## Error-rate constants

The pre-amplification floor (`--rt-error`) is a property of the protocol, so its values are cited
rather than fitted. Nothing here was measured by us except the X2 row.

| Value | What it is | Source | Provenance |
|---|---|---|---|
| 1e-4 per base | V(D)J RT reaction; single-UMI bases get Q40, bases covered by >=2 UMIs get Q60 | 10x Genomics Cell Ranger V(D)J algorithm documentation | vendor-stated |
| 1.54e-4 [1.36e-4, 1.74e-4] | the same floor, fitted from `SRR1763769` as the intercept of `e_out(c) = p_floor + a/c` | X2, `docs/quality_floor.rst`, `scripts/quality_floor.py` | **experimental**, ours; an upper bound (49.6% barcode occupancy inflates it) |
| 7.37e-5 (0.00737%) | TruSight Oncology 500 v2 error rate | Illumina product documentation | vendor-stated |
| Taq 4.3e-5 +/- 1.8; Pfu 2.8e-6; Phusion 2.6e-6; Pwo 2.4e-6, per bp per template duplication | polymerase fidelity, no RT | McInerney P, Adams P, Hadi MZ. *Error Rate Comparison during Polymerase Chain Reaction by DNA Polymerase.* Mol Biol Int 2014:287430. doi:10.1155/2014/287430, PMID 25197572 | published, third party |
| 0.3-6.6e-5 per base per cycle over nine polymerases; **linear-amplification errors 5 +/- 1x the per-cycle PCR rate** | why the FIRST cycle is the one that matters — it is copied into every read of the molecule | Shagin DA, Shagina IA, Zaretsky AR, Barsova EV, Kelmanson IV, Lukyanov S, Chudakov DM, Shugay M. *A high-throughput assay for quantitative measurement of PCR errors.* Sci Rep 2017;7:2718. doi:10.1038/s41598-017-02727-8, PMID 28578414 | published, ours |
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
