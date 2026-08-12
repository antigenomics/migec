# MIGEC v2 — `assemble` (consensus) design

## 0. Position in the pipeline

`refine` hands `assemble` a **group-sorted binary stream** (`.migb`), sorted by `(sample_id, cell_id, umi)`, so `assemble` is a pure streaming consumer: no re-sort, no random access, no index.

```
record: u8 flags | u16 len1 | u16 len2 | packed2bit(seq1) | qual1 | packed2bit(seq2) | qual2
group header: u32 sample_id | u64 cell_code | u64 umi_code | u32 n_reads | f32 p_umi
```
`p_umi` = posterior that this UMI group is a mis-assigned/spurious child, computed in `refine`; consumed here by the quality derate (§3). Fallback path reads grouped FASTQ with `UMI:`/`RX:` headers for interop.

**seqtree is not used in `assemble`.** It indexes a reference set; here every group is 3–30 sequences and building a trie per group costs more than the whole consensus. seqtree belongs in `checkout`/`refine` (barcode correction). Flagging this to prevent architecture-driven overuse.

---

## 1. Consensus algorithm per (sample, cell, umi)

**Chosen: modal-read draft + ungapped quality-weighted offset placement + per-column log-likelihood accumulation.** Rejected: POA (O(L²·n) with a graph allocation per group — at 10⁸ groups this is 100–1000× the budget), MIGEC's centre-anchored 21nt core (arbitrary window, discards signal, and its hard "drop reads with >5 core mismatches" is length-blind).

Per group, in order:

**(a) Pack + dedup.** 2-bit pack with a separate N-mask (`PackedSeq`). Hash identical sequences into a small open-addressing map → unique sequences with multiplicity `c_u`. Quality is *not* collapsed; identical sequences keep their per-read quals (they are averaged correctly by the LL accumulator anyway, so we may keep only the element-wise max-qual representative when `--fast-dedup`; default off).

**(b) Draft.** `draft` = unique sequence with max `c_u`; ties broken by max Σqual, then lexicographic (determinism). Justification: with per-base error e≈10⁻³ and L=150, P(read error-free) = exp(−L·e) ≈ 0.86, and errors are scattered, so no single erroneous variant can outvote the mode for n≥3. For n≤2 the draft choice is irrelevant (the column model decides).

**(c) Offset placement.** For each read pair, best shift `s ∈ [−S, +S]`, `S = --max-offset` (**default 0**: amplicon data has no offset, so this collapses to a no-op and costs nothing). Score is the ungapped log-likelihood of the read given the draft:

```
score(s) = Σ_{j in overlap} [ r_j == d_{j+s} ? log(1 − e_j) : log(e_j / 3) ],  e_j = 10^(−q_j/10)
```
Require `overlap ≥ --min-overlap` (30). Accept `s* = argmax` only if `score(s*) − score(s_2nd) ≥ --offset-margin` (10 nats); otherwise the read is left **unplaced** (counted, not silently dropped). Implementation: 2-bit XOR + popcount Hamming over 2S+1 shifts (L=150 → 5×64-bit words), then the LL is evaluated only at the best 2 Hamming shifts.

**(d) Outlier drop.** Drop the read *pair* (never a single mate — that would desynchronise the FASTQ) if `nmm_total > max(3, --max-mismatch-frac × overlap_total)`, default frac 0.15. Dropping happens **after** sub-clustering (§2), so a genuine second molecule is not thrown away as an outlier.

**(e) Column accumulation.** Frame = `[min(start_i), max(end_i)]` over placed reads (equals the draft frame when S=0). Accumulate, per position j and base b∈{A,C,G,T}:

```
LL[j][b] = Σ_i ( r_ij == b ? log(1 − e_ij) : log(e_ij/3) )
```
N bases and non-covered positions contribute 0.

*Hot-loop identity (this is what makes the throughput target reachable):* write `m_ij = log(1−e_ij)` and `x_ij = log(e_ij/3) = −0.23026·q − 1.0986`. Then
```
LL[j][b] = C_j − acc[b][j],   C_j = Σ_i x_ij  (one running sum per column)
acc[b][j] += (x_ij − m_ij)    only for b = r_ij
```
so the inner loop is **one indexed float add per read per position**, with `x−m` from a 64-entry quality LUT. `argmax_b LL[j][b] = argmin_b acc[b][j]`.

**(f) Base call.** `b* = argmin_b acc[b][j]`. Indels are **not** modelled by default (Illumina indel rate ~10⁻⁶/base, 3 orders below substitution). `--indel-mode banded` adds a banded Myers bit-parallel edit-distance placement used only for reads that fail the ungapped threshold — **recommend cutting from v1**, add when we benchmark IonTorrent/ONT.

Cost: O(n·L) float adds + O(n·(2S+1)·L/64) popcounts. ≈1.5 µs for a 10-read × 150 bp group.

---

## 2. Multiple consensuses per UMI

**Key constraint:** frequency alone cannot separate an early-PCR error from a UMI collision. An error in PCR cycle 1 is present in ~50% of reads. The only reliable discriminator is the **number of divergent positions**: two genuinely different molecules differ at many positions; independent early PCR errors reaching high frequency at D positions simultaneously has probability ~(per-base early-error rate)^D.

**(a) Split threshold, derived not hard-coded.** Per-position probability that a single molecule shows a PCR-derived minor allele at fraction ≥ f (Luria–Delbrück tail approximation):
```
p_pos ≈ 2 · e_pcr / f          e_pcr = --pcr-error-rate (default 1e-4, per-base early-cycle)
E     = L · p_pos              expected such positions per MIG
d_split = min{ D : Σ_{k≥D} e^(−E) E^k / k! < α },  α = --split-alpha (1e-3)
```
L=150, e_pcr=10⁻⁴, f=0.10 → E=0.30 → d_split = **4**. Recomputed per run from the actual read length and CLI parameters; printed in the log.

**(b) Clustering (greedy, count-ordered, no EM).** EM/mixture over the PWM is cut: K is unknown, it needs multiple passes, and it is far more expensive than the signal justifies.
```
sort unique seqs by count desc
seeds = [top]
for each u in order:
    k = argmin_k hamming(u, seed_k)   (at u's best offset)
    if hamming < d_split: assign to k
    else: new seed (cap K ≤ --max-consensus-per-umi, default 4; overflow → nearest, counted)
```
O(U·K·L/64), U ≤ n. **Fast path**: skip entirely if `U == 1` or `n < 2·--min-minor-reads`.

**(c) Split acceptance (likelihood ratio / BIC).** Build each candidate cluster's consensus, let D = number of positions where the two consensuses differ, n = total reads:
```
LL_1 = Σ_j max_b LL_all[j][b]                         (one molecule)
LL_2 = Σ_j max_b LL_S1[j][b] + Σ_j max_b LL_S2[j][b]  (two molecules)
Λ    = 2 (LL_2 − LL_1)
ΔBIC = Λ − D · ln(n)                                   D extra free base parameters
```
Accept the split **iff all four hold**: `D ≥ d_split` **and** `n_minor ≥ --min-minor-reads` (3) **and** `n_minor/n ≥ --min-minor-frac` (0.10) **and** `ΔBIC > 0`. Failing pairs are merged back in reverse count order until the partition is stable. The LL terms make the test quality-aware automatically: a minor cluster supported only by Q15 bases will not clear ΔBIC.

**(d) Pairs are the clustering unit.** For paired data, cluster on the concatenation R1‖R2 so a split always assigns whole pairs and the two output FASTQs stay record-aligned.

**(e) Indexing.** Clusters sorted by read count desc → index `ci = 0,1,2,…`. `ci=0` is the dominant molecule. Emitted in the header as `cI:i:<ci>` and in the read name (§6). `cK:i:<K>` carries the total number of consensuses for that UMI, so a downstream filter `cK==1` recovers exactly the old MIGEC "no collision" behaviour.

**(f) Chimera flag (optional, cut if time).** If the minor consensus differs from the major only in one contiguous block anchored at an end, set `chimera_flag=1` in the stats table. Flag only, never a filter.

---

## 3. Quality recalibration

**(a) Posterior.** Uniform prior over {A,C,G,T}:
```
P(b | obs) = exp(LL[j][b]) / Σ_b' exp(LL[j][b'])
p_cons(j)  = 1 − P(b* | obs) = Σ_{b≠b*} exp(Δ_b) / (1 + Σ_{b≠b*} exp(Δ_b)),  Δ_b = LL[b] − LL[b*] ≤ 0
```
Evaluated with a 4096-entry `exp` LUT over Δ ∈ [−80, 0] (3 lookups per column; a real `expf` here would cost ~4.5×10¹⁰ calls over a 10⁸-group run).

**(b) Error floors, combined additively in probability space** (independent, all small):
```
p_final(j) = p_cons(j) + p_rt + p_umi · δ
Q(j)       = min( --max-qual, round( −10 · log10 p_final(j) ) ),  clamped to ≥ 2
```

- `p_rt` = `--rt-error-rate`, **default 1e-5 → Q50 ceiling** (`1e-6` → Q60). *Why:* an error made by the reverse transcriptase on the first-strand cDNA, or by the polymerase in the first extension before the UMI-tagged molecule is amplified, is present in **every read of the MIG**. Consensus averages away sequencing error and late-cycle PCR error but has exactly **zero** power against a first-cycle error — the MIG is internally consistent and wrong. So no MIG, however deep, may report a base error probability below the first-cycle error rate.
- `p_umi · δ` = the MIG-identity derate. `p_umi` comes from `refine` (posterior that this UMI is a corrupted child of a larger parent that we deliberately kept rather than discarded). If it is a child, some of its reads are leakage from the parent molecule, making the column a mixture; the resulting per-base error rate is the library divergence δ = `--library-divergence` (default 0.05, or estimated as the mean pairwise divergence between consensuses in the sample). Worked example: a kept 4-read orphan with `p_umi = 0.3`, δ = 0.05 → floor 0.015 → **Q18**. Exactly the "medium quality, sequence not lost" behaviour requested. Disable with `--no-umi-derate`.
- Ties: if `Δ_2nd > −--tie-eps` (default 0, i.e. exact tie), emit `N` at Q2. Positions with depth 0 inside the frame → `N` at Q2.

**(c) Small MIGs, explicitly.**
- **n = 1**: consensus = the read. `p_final(j) = 10^(−q_j/10) + p_rt + p_umi·δ`, i.e. the read's own quality, floored. Emitted when `--min-reads 1` (default), tagged `cD:i:1`.
- **n = 2, agreeing**: `p_cons ≈ (e₁e₂/3)/(…)`, ≈ Q(q₁+q₂), immediately truncated by `p_rt` for typical Q30+Q30 → Q50. Correct: two reads of the same molecule give no evidence about RT error.
- **n = 2, disagreeing**: the posterior resolves to the higher-quality base with `p_cons ≈ 10^(−|q₁−q₂|/10)/(1 + 10^(−|q₁−q₂|/10))`; equal quality → `N` Q2.

**(d) Encoding.** Sanger offset 33; `--max-qual` default 60 (`|`). Note for benchmarking: many Illumina-era readers assume ≤Q41; `--max-qual 41` is the compat setting and is what we should use for the arda/bwa-meme comparisons unless verified otherwise.

---

## 4. Cell doublet detection

**What is actually computable here:** we have reads, CB, UMI, and consensus sequences — no gene×cell count matrix and no alignment. That rules out Scrublet/DoubletFinder/scDblFinder-style simulated-doublet classification entirely. Two statistics are genuinely computable and one is genuinely good.

**v1 implements (a) `chain_multiplicity`** — the classic immune-repertoire signal, and the one that matters for MIGEC's actual user base:
- Per CB, group its consensuses by `seq_hash` (or by clustering at `d_split`, `--doublet-key {exact,cluster}`); count molecules (distinct UMIs) per variant.
- A variant is *major* if it has `m ≥ max(2, --doublet-min-frac × M)` molecules, M = total molecules for that CB (default frac 0.10).
- Call: `n_major > --doublet-max-variants` (default 2 — a T cell has ≤2 TRB alleles, a B cell ≤2 IGH).
- Ambient guard (a third variant can be ambient contamination, not a doublet): for the candidate variant with global cross-cell molecule frequency `f_g`, test
  `p = P(Binom(M, λ·f_g) ≥ m)`, λ = `--ambient-inflation` (1.0); call the doublet only if `p < --doublet-alpha / n_cells` (Bonferroni, α=0.01).

**v1 also emits (b) `hop_rate`** (diagnostic, not a call): index-hopping / ambient estimate = fraction of groups whose `(umi, seq_hash)` key also occurs in a different CB. Computed on a **1/64 hash slice** (`hash & 63 == 0`) so the table is ~25 MB rather than ~1.6 GB; `--hop-stats/--no-hop-stats`.

**v1 emits but does not model (c)**: per-CB `n_umi`, `n_reads`, `n_consensus`, `multi_consensus_rate`, ranked → the barcode-rank table (`cells.tsv`), which is what the 10x knee plot is drawn from.

**Honest limits, and the delegation:**
- For whole-transcriptome 3′ GEX, the only signal available at this stage is "this CB has an outlying UMI count" (> 2× median of knee-selected cells). Precision is poor. **We output the barcode-rank table and delegate doublet calling to scDblFinder downstream.** We should not pretend to do GEX doublet detection.
- CellPlex/HTO/deMULTIplex2 demultiplexing is *tag counting*, which is a `checkout` concern (the tag is a barcode in the read), not an assemble concern. **Cut from v1.** Assemble does emit the per-(CB × tag) molecule count matrix if a tag whitelist was declared in checkout (`tags.mtx`), so a later `migec demux-tags` — or an external tool — can fit the negative-binomial model without re-running anything.

Output: `<sample>.cells.tsv` (cell, n_umi, n_reads, n_consensus, n_major_variants, top1_frac, top2_frac, top3_frac, ambient_p, doublet_call) and, under `--tag-doublets`, `DB:i:1` in the FASTQ header.

---

## 5. Optional contig assembly

**The offset layout already is the contig.** Because every read is placed at an offset and the frame is the union of placed reads, a MIG whose reads tile a longer molecule produces a longer consensus for free whenever `--max-offset > 0`. The only genuine gap is a read that overlaps no part of the *draft* and therefore cannot be placed at all.

**Minimal fix (greedy consensus extension, ~20 lines):** after building consensus #1, re-run placement of the still-unplaced reads against the *current consensus* (now longer than the draft); repeat until no frame growth or `--contig-iters` (3).

- **Trigger:** `--contigs` set, `--max-offset > 0`, and unplaced reads exist after pass 1.
- **Ordering constraint:** sub-clustering (§2) runs *first*, and extension never crosses sub-cluster boundaries — otherwise two collided molecules get glued into a fake contig.
- **Failure modes:** (1) internal repeats/homopolymers → mis-placement; mitigated by the `--offset-margin` second-best rule, which leaves ambiguous reads unplaced rather than wrong. (2) coverage gaps → a run of depth-0 columns emitted as `N` (or split into separate records under `--contig-split-gaps`). (3) 5′RACE MIGs with a single deep read stack and one stray read produce a 1× extension tail — mitigated by `--trim-depth` (default: trim terminal columns with depth < min(2, n)).

**Recommendation: keep the union-frame behaviour (free), ship `--contigs` off by default, and do not benchmark it in v1.** De Bruijn / OLC over a MIG is cut — 10–30 reads is far too few for a k-mer graph to beat layout. Note that R1/R2 **merging** is a different thing and stays in `checkout` (where MIGEC's `--overlap` already lives).

---

## 6. Output spec

**FASTQ**, one file per sample per mate: `<sample>.R1.fastq.gz`, `<sample>.R2.fastq.gz` (single file `<sample>.fastq.gz` if merged/single-end). Mates are emitted in lockstep, one record each, identical read names.

```
@<sample>:<cb>:<umi>:<ci>\tRX:Z:<umi>\tQX:Z:<umi_qual>\tCB:Z:<cb>\tCR:Z:<cb_raw>\tMI:Z:<sample>:<cb>:<umi>:<ci>\tcD:i:<n_reads>\tcI:i:<ci>\tcK:i:<n_consensus_this_umi>\tcM:f:<minor_frac>\tcU:f:<p_umi>\tcW:i:<frame_width>\tcS:i:<n_dropped_reads>
```
- Standard SAM tags where they exist (`RX`/`QX` UMI + UMI quality, `CB`/`CR` corrected/raw cell barcode, `MI` molecular identifier) — this is what fgbio, umi_tools, Cell Ranger and arda already read. MIGEC-specific counters use lowercase two-letter tags (SAM reserves lowercase for local use).
- All variable fields live in the **comment**, never in the name, so R1 and R2 names are byte-identical.
- No cell barcode → `CB`/`CR` dropped and name is `<sample>:<umi>:<ci>`.
- `--comment-sep {tab,space}`: tab matches `bwa -C` / `minimap2 -y` passthrough into BAM; verify against bwa-meme during benchmarking (flagging this as a compat detail to confirm, not assume).
- `--header-style migec` reproduces `@MIG.<i> R<n> UMI:<umi>:<count>` for legacy scripts.

**Per-consensus stats `<sample>.mig.tsv.zst`** — one row per emitted consensus:
`sample, cell, umi, ci, n_consensus, reads_total, reads_assigned, reads_dropped, reads_unplaced, cons_len, frame_start, frame_end, depth_median, depth_min, mean_qual, frac_q30, minor_frac, delta_bic, n_diff_to_major, p_umi, derate_q, chimera_flag, contig_gaps, seq_hash`
(`seq_hash` lets the doublet/hop passes work off this table without re-reading FASTQ.)

**Assemble log** — keep 12 of MIGEC's 25 columns, drop the ×2 per-mate duplication (mates are now assembled jointly, so per-mate counters are noise):
- **Keep:** `SAMPLE_ID, INPUT_FASTQ1, INPUT_FASTQ2, OUTPUT_ASSEMBLY1, OUTPUT_ASSEMBLY2, MIG_COUNT_THRESHOLD, MIGS_TOTAL, MIGS_GOOD, READS_TOTAL, READS_GOOD, READS_DROPPED_WITHIN_MIG, MIGS_DROPPED_OVERSEQ`
- **Drop:** `SAMPLE_TYPE`, `MIGS_GOOD_FASTQ1/2`, `READS_GOOD_FASTQ1/2`, `READS_DROPPED_WITHIN_MIG_1/2`, `MIGS_DROPPED_OVERSEQ_1/2`, `READS_DROPPED_OVERSEQ_1/2`, `MIGS_DROPPED_COLLISION_1/2`, `READS_DROPPED_COLLISION_1/2` — collision is no longer a drop, it is a split.
- **Add:** `CONSENSUS_TOTAL, MIGS_SPLIT, MIGS_SINGLETON, MEAN_MIG_SIZE, MEAN_CONSENSUS_Q, D_SPLIT, RT_ERROR_RATE, CELLS, DOUBLET_CELLS, HOP_RATE_EST, WALL_SEC`
Written as TSV (`#` header) **and** as `assemble.json` per sample for programmatic consumption.

---

## 7. Performance

**Threading.** One reader thread parses `.migb` and cuts on group boundaries into **blocks of ~64 groups / ~1 MB**, pushing into a bounded MPMC queue (capacity 4×`--threads`). N worker `std::thread`s pop a block and process every group with a thread-local, preallocated `GroupWorkspace` (packed reads, offsets, `acc[4][W]` floats, unique-seq map, output byte buffer) — nothing allocated in the loop. Workers **bgzf-compress their own output block** so compression is off the single-threaded path; a small reorder buffer keyed by block sequence number makes output byte-deterministic regardless of thread count. Two queues are the only synchronisation.

**Memory.** Per group: reads n·L/4 bytes + 4 floats × frame width (150 → 2.4 KB) + the unique map. Per thread steady-state ~1 MB. `--max-reads-per-mig` (default 10000) caps a jackpot UMI; excess reads are subsampled highest-Σqual-first and counted as `reads_downsampled`. **Target total RSS < 1 GB at 16 threads for any input size** (fully streaming).

**Throughput target.** ≈1.5 µs per 10-read × 150 bp group single-thread → **≥ 500k groups/s at 16 threads**; 10⁸ MIGs (10⁹ reads) in a few minutes of CPU-bound work, with upstream gzip decompression the real bottleneck. Regression gate in `tests/benchmark`: ≥ 3M reads/s/thread of pure consensus work.

**pybind11 must not be in the loop.** Python builds an `AssembleConfig` POD and makes exactly one call into `migec_core::assemble(cfg, progress*)` under `py::gil_scoped_release`; the progress callback fires at most 1 Hz and is the only place the GIL is reacquired. No Python callback per group, per read, or per file.

SIMD is deliberately not in v1 — the accumulate loop is already a single indexed float add; measure before vectorising.

---

## 8. C++ files / classes and CLI

```
include/migec/packed_seq.hpp   PackedSeq, pack/unpack, hamming_shift() (XOR+popcount)
include/migec/mig_stream.hpp   MigRecord, MigGroup, MigBlock, MigBlockReader
include/migec/consensus.hpp    ConsensusConfig, Consensus, ConsensusBuilder, GroupWorkspace
src/consensus.cpp              select_draft(), place_reads(), accumulate_ll(), call_bases()
include/migec/subcluster.hpp   split_group(), d_split_from_poisson(), delta_bic()
src/subcluster.cpp
include/migec/quality.hpp      qual_lut, exp_lut, posterior_perr(), apply_floors(QualCaps)
include/migec/contig.hpp       extend_frame()            [--contigs only]
src/contig.cpp
include/migec/doublet.hpp      CellAccumulator, call_doublets(), binom_sf(), hop_rate()
src/doublet.cpp
include/migec/fastq_writer.hpp BgzfBlockWriter, format_header()
src/fastq_writer.cpp
include/migec/assemble.hpp     AssembleConfig, AssembleStats, assemble(cfg, progress)
src/assemble.cpp               reader thread, worker pool, reorder buffer
src/_bindings.cpp              one assemble(**kwargs) -> dict
tests/cpp/{test_consensus,test_subcluster,test_quality}.cpp   (doctest)
python/migec/assemble.py       thin wrapper -> polars DataFrame of mig.tsv
python/migec/cli/assemble.py   typer command
```

```
migec assemble [OPTIONS] INPUT...
  -o/--output-dir DIR   -s/--sample-sheet FILE   --threads INT(all)   --json FILE

GROUPING   --min-reads INT(1; 0=auto from histogram)  --max-reads-per-mig INT(10000)
           --max-offset INT(0)  --min-overlap INT(30)  --offset-margin FLOAT(10.0)
           --max-mismatch-frac FLOAT(0.15)  --trim-depth/--no-trim-depth(on)
           --indel-mode {none,banded}(none)            [banded: cut from v1]
SPLITTING  --split/--no-split(on)  --max-consensus-per-umi INT(4)
           --pcr-error-rate FLOAT(1e-4)  --split-alpha FLOAT(1e-3)
           --min-minor-reads INT(3)  --min-minor-frac FLOAT(0.10)
QUALITY    --rt-error-rate FLOAT(1e-5)  --max-qual INT(60)
           --umi-derate/--no-umi-derate(on)  --library-divergence FLOAT(0.05)
           --min-cons-qual INT(2)  --tie-eps FLOAT(0.0)
CONTIGS    --contigs/--no-contigs(off)  --contig-iters INT(3)  --contig-split-gaps(off)
DOUBLETS   --doublets/--no-doublets(auto: on iff CB present)  --doublet-key {exact,cluster}
           --doublet-max-variants INT(2)  --doublet-min-frac FLOAT(0.10)
           --doublet-alpha FLOAT(0.01)  --ambient-inflation FLOAT(1.0)
           --hop-stats/--no-hop-stats(on)  --tag-doublets(off)
OUTPUT     --header-style {tags,migec}(tags)  --comment-sep {tab,space}(tab)
           --compress {bgzf,zstd,none}(bgzf)  --stats/--no-stats(on)
```

---

## 9. Cut list (things I would not build in v1)

| Item | Decision | Replacement |
|---|---|---|
| POA / partial-order alignment | cut | offset placement |
| EM / mixture model for sub-clustering | cut | greedy count-ordered + ΔBIC |
| seqtree inside assemble | cut | direct popcount Hamming |
| Indel-aware (banded Myers) placement | defer | ungapped only; revisit for ONT/IonTorrent |
| De Bruijn / OLC contigs | cut | greedy consensus extension, off by default |
| GEX (non-VDJ) doublet calling | cut | emit barcode-rank table, delegate to scDblFinder |
| CellPlex / HTO demultiplexing | cut from assemble | belongs in `checkout`; emit `tags.mtx` only |
| Exact ambient/hop-rate table | cut | 1/64 sampled estimator |
| Chimera flag | optional | cut if schedule slips |
| Per-position / context-specific error models | cut | single `--rt-error-rate` |
| MIGEC's 25-column log | cut to 12 | plus 11 new, plus JSON |
| SIMD in the accumulate loop | defer | measure first |