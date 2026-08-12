# 1. ERRORS FOUND

**E1 — sort: hash-bucketing by `(cell,umi)` makes the UMI correction in `refine` impossible.**
Claim (sort §2.2, §3): "No external merge sort… hash-partition at checkout write time… `migec sort` is not a user-facing command", and (refine §7) "the correction map is applied during the merge phase of the external sort — one pass, no re-sort."
Why wrong: `b = mix64(key) & (nbuckets-1)` sends a UMI and its Hamming-1 child to *different, uncorrelated* buckets, and sort has no merge phase for refine to hook into. After the correction map is applied, a merged child's reads sit in bucket `mix64(child_key)` while the parent's are in `mix64(parent_key)`; the group is permanently split across files. Every downstream group-size, consensus and molecule count is then wrong, and the error is invisible (each half looks like a well-formed MIG).
Correction — two fixes, both cheap:
- Partition on the **top bits of the key**, not a hash: `b = key >> (KEYBITS - log2(nbuckets))`. The 2-bit-packed key is already near-uniform (that is the design's own argument for why hash buckets balance), so range partitioning is equally balanced, costs the same, and additionally delivers the **globally sorted on-disk output by sample/cell/UMI that the user explicitly asked for** and that hash order destroys. Buckets emitted in index order = key order.
- Make the partition key **correction-invariant**: when a CB is present and `--umi-scope cell` (refine's default), partition on `cell` alone — cell barcodes are whitelist-corrected in checkout, so the key is final and all UMI corrections are bucket-local. Without a CB, partition on the top bits of `umi`: a substitution in position `j` moves the key by `< 4^(L-j)`, so all 1-mismatch neighbours except those hitting the top `log2(nbuckets)/2` bases stay in-bucket; handle the boundary by replicating each bucket's *distinct-barcode summary* (not its payload) to its two neighbours.

**E2 — assemble: `--max-offset 0` default, and the assumption that reads sharing a UMI share coordinates.**
Claim (assemble §1c): "default 0: amplicon data has no offset, so this collapses to a no-op and costs nothing."
Why wrong: (a) 10x GEX and 5' VDJ reads sharing a `(CB,UMI)` come from *fragmented* cDNA — they tile the molecule at different start points and frequently do not overlap at all. A depth-0-tolerant "consensus" over them is meaningless, and the frame `[min(start), max(end)]` becomes a chimeric N-padded string. (b) MIGEC Experiment 1 and all SMART/5'RACE data have variable start offsets — that is exactly why MIGEC v1 scanned offsets −5..+5 around a 21 nt core. Shipping `--max-offset 0` silently corrupts both of the pipeline's headline datasets (D2, D7–D9).
Correction: make the placement mode explicit and library-dependent, `--mode {amplicon,fragmented}`:
- `amplicon`: `--max-offset 5` (not 0), ungapped placement as designed.
- `fragmented`: reads of one `(CB,UMI)` are partitioned into **overlap components** (union–find over `overlap ≥ --min-overlap` at the best offset); each component emits its own consensus record, indexed by the same `ci` mechanism as a UMI collision. Components with `n=1` emit the read verbatim. This is exactly the contig path in §5, so it costs no new machinery — but it must be *on* for 10x, not `--contigs` off-by-default.

**E3 — refine §3.2: the per-base consensus error formula is not an error probability.**
Claim: `p_cons(i) = 1 − (n_maj(i)+α)/(c_x+4α)`, α=0.5, described as "the same shape [as MIGEC] but calibrated".
Why wrong: this is the Dirichlet posterior mean of the *minor-allele frequency*, not P(consensus base is wrong). With 10 reads all agreeing at Q30 it returns `1 − 10.5/12 = 0.125` → **Q9**, when the correct posterior error is ~1e-9 (floored at the RT cap). At `c=1` it returns 0.5 → Q3 for *any* read regardless of its Phred; the read's own quality never enters. It is quality-blind, depth-miscalibrated, and would fail bench G3 immediately.
Correction: delete it. Use assemble's column posterior, which is right:
```
LL[j][b] = Σ_i ( r_ij==b ? log(1−e_ij) : log(e_ij/3) )
p_cons(j) = 1 − exp(LL[j][b*]) / Σ_b exp(LL[j][b])
          = Σ_{b≠b*} e^{Δ_b} / (1 + Σ_{b≠b*} e^{Δ_b}),  Δ_b = LL[j][b] − LL[j][b*]
```
`refine` should not compute base qualities at all.

**E4 — refine §3.2 / assemble §3b: the "MIG-identity derate" has the contamination backwards.**
Claim (assemble): "`p_umi` = posterior that this UMI is a corrupted child of a larger parent… if it is a child, some of its reads are leakage from the parent molecule, making the column a mixture", floor `p_umi · δ`, δ=0.05.
Why wrong: if UMI *x* is genuinely an error child of *y*, then **all** of *x*'s reads are reads of molecule *y* — a pure, uncontaminated MIG of *y*'s sequence. The consensus sequence is correct; only the *molecule count* is wrong. The mixture case is the opposite one: *x* is a **real** molecule that additionally received a few leaked error-reads from *y*. So the derate is applied in inverse proportion to the risk it claims to cover. Worked example in the design ("`p_umi`=0.3 → Q18") degrades a base that is almost certainly right.
Correction: the quantity that matters is the expected **foreign-read fraction**, and it is a property of the neighbourhood, not of the error posterior:
```
φ_x = min(1, Σ_{y ≠ x} c_y · ε_j(x|y) / c_x)      # expected leaked reads / own reads
```
and even that need not enter the base quality: foreign reads that disagree show up as a minority allele and the column posterior already down-weights the base, while assemble's sub-clustering splits them. **Recommend cutting `γ_mig` / `p_umi·δ` from the base-quality path entirely** (removes two tunables, `--library-divergence` and `--no-umi-derate`) and instead emitting the molecule-identity posterior as a *tag* — `pE:f:<P(error-child)>` — for downstream molecule-count weighting. This also removes the double counting between refine's `γ_mig` and assemble's `p_umi·δ`, which currently apply the same derate twice.

**E5 — checkout §6 / refine §4: effective UMI space computed with Shannon entropy; birthday collisions require collision entropy.**
Claims: checkout `umi_effective_bits = Σ_p H_p`, `S = 2^(Σ H_p)`; refine `L_eff = (2L − Σ_j D_j)/2`, `s = n/4^L_eff`.
Why wrong: the birthday functional is `P(two independent molecules share a UMI) = Σ_u p_u²`, i.e. Rényi entropy of order 2, not order 1. Since `H_2 ≤ H_1`, using Shannon **overestimates the usable space and underestimates collisions**, in the direction that silently loses molecules. (Refine's own `P1` correctly uses `m_j = Σ_a p_j(a)²`, so the two metrics in the same file are mutually inconsistent.)
Correction (position-independent model):
```
m_j     = Σ_a p_j(a)²                      # per-position collision prob (=1/4 uniform)
P_coll  = Σ_u p_u² = Π_j m_j
S_eff   = 1 / Π_j m_j                      # report this, not 2^(Σ H_p)
L_eff   = −log4( Π_j m_j ) = −Σ_j log4(m_j)
E[#colliding pairs among M molecules]   = C(M,2) · Π_j m_j
E[#UMIs carrying ≥2 molecules]          = Σ_u (1 − (1+M p_u) e^{−M p_u})
                                        ≈ (M²/2) · Π_j m_j     for M·max p_u ≪ 1
```
Keep `Σ_j H_p` only as a display statistic. Note further that position independence itself is an assumption (oligo synthesis produces position-correlated bias); `Π_j m_j` is therefore still a *lower bound* on the true collision rate — measure it (Experiment 1).

**E6 — refine §1.2/§1.3: the error model is sequencing-only, so PCR-derived UMI children are kept as real molecules.**
Claim: `eps_j(x|y) = κ·p̄_j(x)/3` from Phred; `Lerr(c_x|y) = ZTPoisson(c_x; c_y·eps_j)`.
Why wrong: a UMI substitution introduced by the polymerase in PCR cycle 1–3 is present in ~50%/25%/12% of the descendants and carries **high** Phred in every read. With `eps_j ~ 3e-4`, `Poisson(c_x; c_y·3e-4)` assigns essentially zero probability to `c_x = 0.4·c_y`, so the posterior sends it to the "independent molecule" branch and it is emitted as a second molecule. This is the dominant residual error in UMI counting, and it is precisely what umi_tools' `directional` rule (`c_y ≥ 2c_x − 1`) exists to catch — the design's own comparison table dismisses that rule as "a special case of ours", when in fact it covers a case ours does not.
Correction: make the error component a two-part mixture (sequencing + polymerase), which restores the count-ratio evidence the directional rule uses:
```
ρ    = 3L · ε_pol · n_cycles                       # P(this neighbour is polymerase-derived)
f    = c_x / (c_x + c_y)                           # child's fraction of the family
LD(f) ≈ 2 ε_pol / f²    (branching-process / Luria–Delbrück tail density; ∫ over f)
Lerr(c_x | y) = (1−ρ)·ZTPois(c_x; c_y·ε_seq,j) + ρ·LD(c_x/(c_x+c_y))
```
Default `ε_pol` from the same source as the assemble split threshold (one constant, two uses). This also fixes the ε estimator: the neighbour method-of-moments in §1.2 currently attributes *all* excess D1 pairs to sequencing error, so `ε_nbr` and hence `κ` absorb the polymerase term and are biased high.

**E7 — refine §5.1: `f_err = ZTNB(c; r, ε·μ)` is the wrong distribution for error-MIG sizes.**
Why wrong: an error child's size is `Poisson(ε·C)` where `C` is the *parent's* size drawn from `f_real` — large parents make large error children, so `f_err` is a Poisson mixture over `f_real`, which is far heavier-tailed than an NB with mean `ε·μ`. Under-modelling the tail is exactly the direction that makes `m_hi` (the FDR threshold) too small.
Correction (exact, one convolution over the size histogram, same cost):
```
f_err(c) = [ Σ_C f_real(C) · Pois(c; ε·C) ] / [ 1 − Σ_C f_real(C) e^{−εC} ],   c ≥ 1
```
Still three free parameters `(w, r, μ)`; ε comes from the estimator. Weight `C` by 3L parents-per-child if enumerating children rather than pairs.

**E8 — refine §1.2: `D1_obs` double-counts sibling pairs, inflating the estimated error rate.**
Claim: `E[children](ε) = Σ_i 3L(1 − e^{−c_i ε})`, solved against `D1_err = D1_obs − D1_ind`.
Why wrong: `D1_obs` counts *unordered distance-1 pairs among distinct barcodes*. Two children of the same parent that differ at the **same** position by different bases (`AAAA→CAAA`, `AAAA→GAAA`) are at Hamming distance 1 **from each other** and are counted as an extra pair. For deep MIGs (`c_i ε ≫ 1`) all three alternatives at every position are observed, adding `3L` sibling pairs per parent — the same order as the child term, so ε is overestimated by up to ~2×.
Correction:
```
E[D1_obs](ε) = C(n,2)·P1
             + 3L · Σ_i (1 − e^{−c_i ε})          # parent–child pairs
             + 3L · Σ_i (1 − e^{−c_i ε})²         # sibling pairs (C(3,2)=3 per position)
```
Bisect on `log ε` against this. (Grandchild pairs contribute at O(ε²) and can be neglected.)

**E9 — refine §5.1: the over-sequencing rule inverts on deep data.**
Claim: "over-sequenced iff `μ ≥ 5 and w < 0.5`", where `w` = mixture weight of the error component over *distinct UMIs*.
Why wrong: `w` grows monotonically with depth. Expected distinct children per parent `= 3L(1 − e^{−cε})`; at L=12, ε=3.3e-4, `c=20` → 0.24 → `w≈0.19`; at `c=1000` → 10.1 → `w≈0.91`. So a *massively* over-sequenced library is classified as not over-sequenced. Correction: drop `w` from the test. Over-sequenced iff `μ ≥ 5`; additionally report the monotone, interpretable `frac_reads_in_MIGs_with_c≥5` and `mean_reads_per_umi = R/n_corrected`.

**E10 — assemble §2: `d_split` tests polymorphism, not linkage; and the ΔBIC penalty is far too weak.**
Claims: `p_pos ≈ 2ε_pcr/f`, `E = L·p_pos`, `d_split = min{D : P(Pois(E) ≥ D) < α}` = 4; accept split iff `ΔBIC = 2(LL₂−LL₁) − D·ln(n) > 0`.
Why wrong: (a) The Poisson test asks "how likely are D polymorphic positions in one molecule?" but the split hypothesis requires those D minor alleles to be **on the same reads**. D independent PCR subclones at similar frequency almost never co-segregate; conversely a single early subclone with D linked errors is a different (and much rarer) event. The discriminating statistic is linkage, not the count of polymorphic sites. (b) BIC: `n` is taken as the number of *reads* and the penalty as `D` parameters. The correct observation count is bases (`n_reads·L`) and the second cluster adds `L+1` free parameters plus `n_reads` latent assignments. As written, a *single* read carrying 5 real sequencing errors gives `Λ ≈ 2·5·7 ≈ 70` against a penalty `5·ln(20) ≈ 15` → **the spurious split is accepted**. Only `--min-minor-reads 3` prevents it, i.e. the ΔBIC term is doing no work.
Correction: cut `d_split` (it is redundant once linkage is tested) and replace the penalty:
```
ΔBIC = 2(LL₂ − LL₁) − [ (L+1)·ln(n_reads·L) + n_reads·ln 2 ]
```
and add the explicit linkage/coherence test that is the actual discriminator — per divergent position *j*, with mean per-base error `ē`:
```
p_j = P( Binom(n_reads, ē/3) ≥ n_minor )        # all n_minor reads carry the same minor base at j
p_link = Π_{j∈D} p_j                            # requires the SAME n_minor reads at every j
accept iff  p_link < α / (#MIGs tested)         # Bonferroni over the run
```
Then calibrate α empirically by the within-MIG re-partition null (Experiment 3) rather than trusting the Poisson derivation.

**E11 — assemble §4a: the doublet rule needs chain assignment the pipeline does not have.**
Claim: "`n_major > --doublet-max-variants` (default 2 — a T cell has ≤2 TRB alleles, a B cell ≤2 IGH)."
Why wrong: `n_major` as defined counts **all** distinct consensus variants in the cell, with no locus assignment — because the pipeline deliberately ends before alignment, there is no V/J call and no way to know which consensus is TRA vs TRB vs an off-target transcript. A normal T cell in a 5' VDJ library routinely yields ≥4 major variants (≥1 TRB, 1–3 TRA — allelic exclusion at TRA is leaky — plus off-target). The rule as written calls essentially every cell a doublet. The ambient guard makes it worse in the direction that matters: `p = P(Binom(M, λ f_g) ≥ m)` uses the sample-wide frequency `f_g`, so doublets involving an **expanded clone** (the most common kind) have large `f_g` and are systematically missed.
Correction: define loci without alignment by clustering the sample's consensus set (they already compute `seq_hash`; cluster at, e.g., ≥90% identity over the overlap) and count *distinct lineages*, not distinct variants; call a doublet only when ≥2 lineages each hold ≥20% of the cell's molecules **and** the lineage pair is not the sample's dominant pair. Replace the global `f_g` with a *rank-conditioned* ambient rate estimated from below-knee barcodes only. Given the ambiguity, bench's recommendation is right: ship this as unvalidated, validate on `synthetic/tenx_sim` only.

**E12 — index hopping is claimed in two places and detected in neither.**
Claims: checkout §3 `ambig` bucket = "≥2 sample rows matched with ΔS < Δ (index hopping)"; assemble §4b `hop_rate` = fraction of groups whose `(umi, seq_hash)` also occurs in another CB.
Why wrong: (a) A hopped read carries a *valid, unambiguous* index — it matches exactly one sample row with a perfect score. The `ambig` bucket catches sequencing errors in the index, not hopping, and will report a hop rate of ~0 for a library with real hopping. (b) The `hop_rate` estimator is confounded by clonal expansion: for a clone present at 1e5 molecules across cells, the number of same-sequence same-UMI pairs *by chance* is `C(1e5,2)/4^12 ≈ 300`, all in different CBs, all counted as hops. In exactly the RepSeq/VDJ datasets this project targets, the estimator measures clonality, not hopping.
Correction: (a) Estimate hopping the only way it is estimable — from **unused index combinations**. With dual indices `i7×i5` and a sheet declaring `K` used combinations out of `|I7|·|I5|`, the read counts landing in unused combinations give a direct rate: `hop_rate ≈ (Σ_{unused} reads) · (|I7|·|I5| − K) / [ K · Σ_{used} reads ]`… simply: report the unused-combination read fraction, normalised by the number of unused cells. This requires `checkout` to keep the *full i7×i5 contingency table*, which it currently discards — one small addition. (b) For the CB-level estimator, restrict to `seq_hash` observed in ≤ k cells, or drop it.

**E13 — checkout §2: the mismatch term of the LLR carries a spurious `1/(4−m)`.**
Claim: `Σ_i [b_i ∉ S_i] · w_i · log2( 4 e_i / (3(4−m_i)) )`, with quoted values −11.0 bits at q30, −4.5 at q10.
Why wrong: with true base `t` uniform over `S_i` and a specific observed base `b ∉ S_i`, `P(b | adapter) = Σ_t (1/m_i)(e_i/3) = e_i/3` — the probability is *not* spread over the `4−m_i` outside bases; each of them individually has probability `e/3`. The null is `1/4`. So the ratio is `4e_i/3`, and the extra `1/(4−m_i)` overweights every mismatch by a factor of 3 (m=1) to 2 (m=2).
Correction:
```
S(o) = Σ_{b_i∈S_i} log2( 4[(1−e_i)/m_i + (m_i−1)e_i/(3 m_i)] )
     + Σ_{b_i∉S_i} w_i · log2( 4 e_i / 3 )
```
Corrected numbers for `m=1`: match +2.00 bits, mismatch at q30 **−9.55** bits (not −11.0), q10 **−2.91** (not −4.5), q2 (**N**) −0.60. The threshold `S_min ≈ 20.5` bits should then be re-derived; and it should count only the offsets actually scanned, `(L_read − L_pattern + 1)·N_patterns`, not `L_read·N_patterns`. Note also that `α` is called `--fdr` but is a per-read family-wise rate, and the null "read is i.i.d. uniform ACGT" is badly violated by real reads (shared primers, composition bias) — calibrate `S_min` empirically against shuffled/decoy patterns rather than analytically.

**E14 — checkout §4 / refine §1.6: the whitelist posterior drops the `(1−e)` odds and the "not in the whitelist" hypothesis.**
Claim: `P(c|obs) ∝ π_c · e_i/3`.
Why wrong: candidates differ from `obs` at *different* positions `i`, so `P(obs|c) = (e_i/3)·Π_{j≠i}(1−e_j)`; the common factor `Π_j(1−e_j)` cancels only after dividing by `(1−e_i)`. The design drops that, which biases toward correcting at high-quality positions exactly where the correction is least likely — and the bias is largest at low `q`, which is where corrections actually occur.
Correction:
```
P(c | obs) ∝ π_c · e_i / (3(1 − e_i))          # odds, not probability
π_c = (n_c + 1)/(Σ n + |W|)
```
Additionally, the posterior must be normalised against a background hypothesis "the true barcode is not in the whitelist" with weight `π_0 · Π_j(1−e_j)`; without it, a read from an undeclared sample or a hopped index is always assigned to *some* whitelist entry with posterior 1.0. Accept iff `max_c P ≥ 0.975` **after** including `π_0`.

**E15 — the RT/PCR quality cap is asserted, not measured, and the three designs disagree with each other by 10×, inverted.**
Claims: assemble `--pcr-error-rate 1e-4`, `--rt-error-rate 1e-5`; repo simulator `pcr_error_rate=1e-5`, `rt_error_rate=1e-4`; refine `--rt-error 1e-5 → Q_max 50`; repo's synthetic assertion "must not claim a quality above `−10log10(rt + pcr·cycles)`" → **Q36** with its own defaults, contradicting the Q50 cap that assemble emits.
Why wrong: the cap is the single most consequential constant in the project (it sets every emitted quality above ~Q40) and it is currently a guess that two designs write in opposite directions. Also, the correct floor is not `ε_RT` alone: an error at PCR cycle *k* reaches the consensus whenever its descendant fraction exceeds ½, contributing `≈ 2ε_pol` on top of the RT term.
Correction:
```
p_floor = ε_RT + 2·ε_pol            # not ε_RT alone
Q_max   = −10 log10(p_floor)
```
Make `--rt-error auto` the default and **derive it from the data**: on a clonal control (bench D5) the residual consensus error at large MIG size *is* `p_floor`, so fit `e_out(c) = p_floor + a/c` over MIG-size bins and take the intercept. Refuse to emit any Q above the measured value. Set a single constant in one header consumed by refine, assemble and the simulator.

**E16 — assemble: nothing normalises read orientation within a MIG, yet checkout deliberately creates mixed-strand MIGs.**
Claim: checkout `--rc-search` "searches each pattern's reverse complement in each mate; strand recorded as `ZS:A:+/−`" and `--umi-canonical duplex` emits `UMI = min(α∥β, β∥α)` so that "the two single-strand families remain recoverable". Assemble never mentions strand.
Why wrong: both features guarantee that a single `(cell,umi)` group contains reads in *both* orientations. Assemble's draft/placement/column accumulation is orientation-blind, so the reverse-complement reads mismatch at ~75% of positions, are dropped by the `nmm > 0.15·overlap` rule (silently losing half the MIG) or, worse, are split off as a second "molecule" and emitted as a spurious consensus.
Correction: sequences must be strand-normalised **before** the group is written (checkout reverse-complements and reverses the quality string when `ZS:A:−`, keeping the flag for provenance) — the `.mig` `flags` field already reserves `b0 rc1, b1 rc2`, so make it *descriptive of what was already applied*, not a to-do. If duplex is to be supported properly, assemble must first build the α and β **single-strand consensuses separately** and only then compare them; see M4.

**E17 — sort: `--max-reads-per-umi` cannot be applied where it is placed.**
Claim: sort §3 and its CLI put `--max-reads-per-umi 100000` (reservoir downsample) on `migec checkout`, at the IO write layer.
Why wrong: at write time the records are unsorted and streaming across multiple writer threads; the count for a UMI is unknown until after grouping, so no per-UMI reservoir can exist. It is also stated to be seeded/deterministic, which is impossible with a nondeterministic multi-writer block order.
Correction: move the cap into assemble, applied **after** sub-clustering (the design's own open question, option (a)) — which is also the statistically correct answer, since capping before the split biases the group-size ratio that the split test depends on. Downsample deterministically by `(highest Σqual, then src_index)`, not by a reservoir.

**E18 — memory estimates understated.**
(a) refine §2.1: "~28 MB at 1M, ~400 MB at 10M". Per barcode is `packed 8 + count 4 + meanq L=12` = 24 B; the hash at the stated 2× load factor is `2·16·n` = 32 B/barcode. Total ≈ **56 MB at 1M, 560 MB at 10M** — the text computes `16n` for a table it declares to be at 2× load. (b) assemble §7: "Target total RSS < 1 GB at 16 threads for any input size (fully streaming)" is contradicted by sort §2.3, which decompresses **an entire bucket into one contiguous arena** sized `target_bucket_bytes = --mem/--threads` — 1 GB per thread by construction, so 16 GB at 16 threads. Correction: state the real budget, `RSS ≈ threads × bucket_arena + O(1)`, and set `nbuckets` from it; delete the "<1 GB for any input" claim.

**E19 — checkout and sort contradict each other on the on-disk format, and both quote unverified compression numbers.**
checkout `.mgc` uses 2-bit packed sequence + binned quality and claims "~2.5–3× smaller than gzip FASTQ"; sort `.mig` measures 2-bit packing as **worse** (227 vs 197 B/pair) and stores raw ASCII; assemble `.migb` again uses `packed2bit(seq)`. Three formats, three record layouts, one pipeline. Additionally sort's "raw + zstd-1 → ~72 B for 300 bases" is 1.92 bits/base, *below* the i.i.d. ACGT entropy — achievable only by exploiting cross-read redundancy, which is real for amplicon libraries and largely absent after hash-bucketing scrambles read order (E1) or for GEX. Correction: one format; store raw ASCII per sort's arithmetic; **lay blocks out column-major** (all `seq1`, then all `seq2`, then all `qual1`, then all `qual2`) rather than interleaving sequence and quality per record — interleaving two very different distributions costs zstd 10–20%; and quote the ratio as measured-per-dataset, not as a design constant.

**E20 — nominal Phred is used everywhere, while the one design that measures a calibration table has no consumer.**
checkout §6 emits `checkout.qual_calibration.tsv` (`ê(q)` from mismatches against constant segments). refine and assemble both use `e = 10^(−q/10)` directly. On NovaSeq/NextSeq RTA3 there are only ~4 distinct Q values, so `ê(Q37)` spans an order of magnitude by cycle and by context — and the LLR's whole design intent ("high-quality mismatches are lethal, low-quality ones near-free") collapses when 95% of bases carry the same Q. Correction: make `ê(q)` (optionally `ê(q, cycle)`) a first-class artefact written by checkout, carried in the `.mig` file header, and consumed by (i) the LLR, (ii) `ε_j` in refine, (iii) `LL[j][b]` in assemble. Refine's `κ` clip to `[0.2, 20]` becomes a per-`q` table instead of a scalar. This is a small change and it is the difference between passing and failing bench G3 on NovaSeq data.

**E21 — repo §8: `--order first` UMI subsampling is not unbiased.**
Claim: "`--order first` = the first N distinct (deterministic, no seed needed)".
Why wrong: a UMI with 100 reads is ~100× more likely to make its first appearance early than a UMI with 1 read, so first-appearance order **oversamples large MIGs** — which destroys exactly the MIG-size distribution the example notebooks are built to display. (bench correctly rejects lexicographic-first for a different reason but proposes the right fix.) Correction: use bench's hash partition everywhere and delete `--order first`:
```
keep(umi) ⟺ blake2b_64(umi_seq, key="umi_data-v1") mod 10000 < K
```
Unbiased, streaming, order-independent, reproducible without a UMI table. Keep repo's exit assertion (`mean_reads_per_umi` preserved) as the guard.

**E22 — refine cuts EmptyDrops; bench gates cell calling against a reference that used it.**
refine §6: OrdMag only, "EmptyDrops — cut". bench G6: `jaccard_cells ≥ 0.98` against Cell Ranger 6.0.1's `filtered_feature_bc_matrix`. Cell Ranger ≥3 calls cells with OrdMag **plus** an EmptyDrops-style second pass that rescues low-RNA barcodes; OrdMag alone systematically misses that rescued set, so the gate is unreachable by construction. Also `per_read_cb_concordance ≥ 0.995` requires reproducing Cell Ranger's abundance prior exactly. Correction: either keep the delegation and set the gate on the OrdMag-comparable subset (`jaccard` against Cell Ranger's OrdMag-only call, or `recall of Cell Ranger cells ≥ 0.98` with the FP set reported separately), or lower to `jaccard ≥ 0.90`, `per_read_cb_concordance ≥ 0.99`, and report the disagreement breakdown by barcode rank. As written G6 will fail for a reason unrelated to our code quality.

**E23 — sort: `src_index` is `uint32`, and it is load-bearing for determinism.**
`uint32` caps at 4.29e9 read pairs; a NovaSeq X run exceeds this. On overflow the tiebreak key duplicates and the "byte-identical at 1 vs 8 threads" guarantee (bench G8) fails nondeterministically — the worst possible failure mode. Correction: `uint64`, or `uint32` plus an explicit hard error at 2^32 records with the actual count in the message.

**E24 — refine §1.7: per-cell UMI scope over-merges, because gene assignment is out of scope.**
Claim: "default `cell` when cell barcodes are present… per-cell `n` is 10³–10⁴ (correction is aggressive and correct)".
Why wrong: Cell Ranger corrects UMIs within `(cell, gene)`; we have no gene. With a 10-nt UMI (10x 5' v2, the design's own D7 chemistry) and 5,000 UMIs in a cell:
```
λ_ind = n·3L/4^L = 5000·30/1,048,576 = 0.14      (0.29 at n=10,000)
```
i.e. 14–29% of UMIs in a cell have a *genuine* 1-mismatch neighbour from a different transcript, and the count evidence alone will merge a good fraction of them. Correction: within a cell, require **sequence compatibility** before merging — which means the merge decision cannot be made in `refine` before `assemble` sees the reads. Concretely: let refine emit *candidate* merges with their posteriors, let the grouping engine co-locate candidate parent/child (E1's range partition does this for free when the substitution is not in the top bases), and let assemble accept or reject the merge using the same `ΔBIC`/linkage test it already runs for splits. That inverts the current philosophy (refine merges conservatively at τ=0.95, assemble splits) into the correct one: **merge liberally, split on sequence evidence** — with the bonus that the two stages then share one test instead of having two mutually inconsistent ones.

**E25 — bench G4 prescribes a quality that the corrected model cannot produce.**
Gate: "MIGs with 3–5 reads and no Hamming-1 parent … are emitted at `Q ≤ 30`". Three agreeing Q30 reads give `p_cons ≈ 1e-9`, floored at `p_floor` → Q40–Q50. `Q ≤ 30` is achievable only via the γ derate, which E4 shows is wrong. The gate therefore hard-codes the defect. Correction: replace with a *calibration* criterion, which is falsifiable and model-free — for the 3–5-read bucket, `ê(Q) ≤ 2·10^{−Q/10}` for every Q bucket with `n_Q ≥ 1000` (i.e. G3 restricted to that bucket), plus retention `≥ 95%` and `e_out` reported, not bounded a priori.

---

# 2. MISSING COVERAGE

**M1 — Paired-read overlap/merge is assigned in a circle and implemented nowhere.**
checkout §10 cuts it: "belongs in assemble (where contigs are built anyway) or is delegated to fastp/bbmerge". assemble §5 says: "R1/R2 **merging** is a different thing and stays in `checkout` (where MIGEC's `--overlap` already lives)." Neither builds it. MIGEC v1 had it, and it matters for 2×100 on a ~250 bp amplicon (D1, D3, D4 are all in this regime) and for the D5 single-end 600 bp merged design.
Proposal: implement it in **assemble**, not checkout, as a special case of the machinery already specified — after strand normalisation, place R2's reverse complement against R1 with the existing `place_reads()` offset search over `[−S, L]`, accept if `overlap ≥ --min-overlap` and `score(s*) − score(s₂) ≥ --offset-margin`, and merge into one frame. This is ~15 lines on top of §1(c), reuses the LL column accumulator (higher-quality base wins automatically, and the merged base quality is the correct combined posterior rather than MIGEC's max), and removes four flags and a whole second matcher from checkout. Emit `<sample>.fq.gz` single-file output when merged, which the output spec already anticipates.

**M2 — Chimeras / PCR recombination.**
Only a per-MIG "flag only, cut if time" in assemble §2f, which detects a chimera *within* one MIG. The dominant chimera source in multiplex amplicon RepSeq is **cross-molecule PCR recombination**, producing a clean, high-confidence consensus that is a hybrid of two abundant clonotypes — indistinguishable at the MIG level and a primary false-positive source in bench §5.3 (`FP_count`, gate G5).
Proposal: a sample-level post-pass over `<sample>.mig.tsv` (no new IO): for each consensus with `n_MIGs` below a threshold, test whether a breakpoint `k` exists such that `prefix[0:k]` matches an abundant consensus A and `suffix[k:]` matches an abundant consensus B, both at ≥99% identity, with `min(count_A, count_B) ≫ count_self`. Report `chimera_parents` and `breakpoint` as columns; **flag, never filter** (consistent with the project's keep-everything stance). ~60 lines, and it is the difference between passing and failing G5.

**M3 — `N` in a cell/sample barcode is discarded rather than corrected.**
checkout defaults `--umi-max-n 0` and sort asserts "an `N` in a barcode never reaches here". Cell Ranger treats an `N` in the CB as a free wildcard and corrects it against the whitelist. Discarding is a systematic, chemistry-dependent read loss and will show up directly in bench's `per_read_cb_concordance`.
Proposal: for a whitelisted tag, an `N` at position `j` expands to 4 candidates at that position with `e_j = 0.75` (uniform prior over bases), scored by the same posterior as E14; accept at the same threshold. For a non-whitelisted UMI, keep the read, record `umi_has_N` (the flag already exists), and let the neighbour search treat `N` as matching all four bases (a UMI with one `N` is exactly a 1-edit query — this is the case where seqtree's variable-cost search is genuinely the right tool).

**M4 — Duplex (α/β) support is extraction-only; no single-strand→duplex consensus.**
checkout provides `--umi-canonical duplex` and bench gate G9 is scoped to "the UMI histogram… within 2%", i.e. extraction only. But the entire point of a duplex tag is that the two strand families are consensused **separately** (SSCS) and then compared to give the duplex consensus (DCS), which is what suppresses damage and early-cycle errors — the very thing the flat `p_floor` cannot fix.
Proposal: either (a) state explicitly in docs and CLI help that v1 supports duplex *tags* but emits single-strand consensuses (and then do not use duplex data to justify the error-suppression claim), or (b) implement it as a 30-line addition on top of the existing sub-clustering: the strand flag partitions the group before consensus, both SSCS records are emitted with `ZS:A:+/−`, and where both exist emit a third DCS record with per-base `p = p_α·p_β/(p_α·p_β + (1−p_α)(1−p_β)·⅓)` at agreeing positions and `N` at disagreeing ones (fgbio's `aD/aM/aE` tags for the duplex counters). Recommend (a) for v1, with (b) explicitly in ROADMAP — but the *decision* must be made, because right now the design implies (b) and delivers (a).

**M5 — "on-disk sort by sample/cell/UMI barcode" — the user's explicit request — is deleted.**
sort §7: "`migec sort` — **not exposed.**" and the output is in hash order. See E1: range partitioning on the top key bits restores true key order at zero cost. Expose `migec sort` as a documented command producing key-ordered `.mig` (it is then a no-op wrapper over the partitioner + per-bucket `pdqsort`, ~20 lines of CLI). This is also what makes `migec view --group CELL:UMI` a seek rather than a scan, and what lets two runs be diffed.

**M6 — seqtree, explicitly requested, is cut by all four algorithmic designs.**
checkout cuts it for whitelists ("3L hash probes are faster"); refine cuts it for UMIs ("the single biggest simplification available"); assemble cuts it ("architecture-driven overuse"); repo keeps it as a Python runtime dependency that nothing on a hot path calls. Each individual argument is correct, but the net result is that a user-mandated dependency is vestigial.
Proposal: name the three places where it is genuinely the right primitive and commit to them: (i) whitelist lookup at `max_substitutions = 2` and for `N`-containing barcodes (M3), where enumeration is `9L²/2 ≈ 1152` probes at L=16 and the trie wins; (ii) variable-length / indel-bearing UMIs (`--umi-indels`, Ion Torrent, ONT, homopolymer-adjacent anchors), where enumeration cannot express the query at all; (iii) the adapter/primer **discovery** step in `migec suggest` §6.3 — building an index over the top conserved 12-mers and searching reads against it is exactly `Index::build` + `search_batch`, and currently that step is specified as a bespoke greedy extension. Otherwise, tell the user it was evaluated and dropped, with the measurement.

**M7 — Sample-barcode error correction without a whitelist.**
Requirement: "guesses and corrects sample/cell/UMI barcodes using an error-rate model". Sample and cell barcodes are only ever corrected against a whitelist (checkout §4, refine §1.6). For a MIGEC-style `barcodes.txt` the sample set *is* a whitelist, so this is mostly covered — but the "identify missing/absent barcodes from a list" requirement is met only for sample/cell, and there is no path for the case of a declared-but-unlisted barcode (a sample sheet typo, a swapped i5/i7). Proposal: reuse the E14 posterior with `π_0` — a read whose best sample-barcode posterior is dominated by the background hypothesis goes to `undef-sample` and its raw barcode is aggregated into `checkout.unassigned_barcodes.tsv` (top 100 by count). A swapped index pair or a typo then appears as a single high-count row, which is the actual diagnostic users need and which no current output provides.

**M8 — Separating sequencing error from PCR error in the reported "error rate".**
Requirement: "estimates error rates". checkout estimates a *sequencing* error from constant segments; refine estimates a lumped `ε_nbr` from UMI neighbours which (per E6) silently mixes sequencing and polymerase error. These are different constants used for different purposes (`ε_seq` → LLR and consensus posterior; `ε_pol` → quality floor and split threshold), and confusing them propagates into `κ`, `m_hi`, `d_split` and `Q_max`.
Proposal: report both, and separate them by the one signal that distinguishes them — **quality dependence**. Sequencing errors correlate with Phred; polymerase errors do not. Fit `ê(q) = ε_pol_effective + a·10^{−q/10}` over the constant-segment mismatch counts stratified by reported `q`; the intercept is the quality-independent (polymerase + RT + damage) component and the slope calibrates `κ`. One regression over a table checkout already builds. Write both to `qc_summary.tsv` as `err_seq`, `err_quality_independent`.

**M9 — No sample-level QC gate for the case where the pattern is simply wrong.**
`migec suggest` is excellent but is a separate opt-in command; a user who supplies a wrong offset gets a low match rate and no diagnosis. Proposal: when `reads_matched/reads_total < --min-match-rate` (default 0.5), `checkout` automatically runs the §6 profiling on the first 10⁵ reads and prints the suggested pattern in the error message. Zero new machinery, and it converts the most common user error from a silent 10%-yield run into a one-line fix.

**M10 — Molecule-count correction for undetectable collisions.**
The multi-consensus-per-UMI feature can only ever split collisions between *different* sequences. Two molecules of the same clonotype sharing a UMI (guaranteed in RepSeq and in any expanded clone) are undetectable, so molecule counts are systematically biased low, by `≈ C(M,2)·Π_j m_j / M = (M−1)·Π m_j /2` per molecule — at M=1e6 and L=12 that is ~3%. Nothing in the design corrects or reports it. Proposal: report the model-based correction alongside the raw count, `M̂ = S_eff·(−ln(1 − M_obs/S_eff))` (the Poisson/rarefaction inverse, with `S_eff` from E5), plus the collision-adjusted per-clonotype count; and make bench's `mol_count_rel_err` compare against `M̂`, not `M_obs`.

---

# 3. The three cheapest experiments that would falsify the riskiest assumptions

**X1 — Read-start dispersion within `(CB,UMI)` on one 10x GEX+VDJ run. Falsifies E2 (the ungapped, co-terminal-reads assumption) and the entire consensus premise for single-cell.**
One pass over ~5M read pairs of D7. Extract `(CB, UMI)`, map R2 to the transcriptome (or simply cluster reads by their first 30 nt), and tabulate, per `(CB,UMI)` group with `n ≥ 3`: the spread of read start positions and the fraction of read *pairs within a group* whose best ungapped overlap is `< 30 nt`. Prediction under the current design: near 0. Prediction under E2: the majority. Cost: one afternoon, no new code beyond a 40-line script over `.mig`. If the dispersion is high — and it will be for 3' GEX — `--max-offset 0` and single-consensus-per-UMI must be replaced by the overlap-component path before any 10x number in the plan (G6) is meaningful. Run this **first**; it is the single assumption whose failure invalidates the most of the design.

**X2 — Emitted-quality calibration on a clonal control, stratified by MIG size, with the derate on and off. Falsifies E3, E4, E15, and gate G4 simultaneously.**
Take D5 `SRR1763769` (8E5, clonal) or D1's spike-in, mask real variant sites *from the data itself* (any position whose across-MIG consensus allele frequency exceeds 20% is real, not error — do not rely on a variant metadata file the cell line does not have), then compute
```
ê(Q) = n_errors(Q)/n_bases(Q)   vs   10^(−Q/10),   per MIG-size bin {1,2,3–4,5–9,10–19,20–49,≥50}
```
and fit `e_out(c) = p_floor + a/c` over the bins; the intercept **is** the RT/PCR floor and settles the 1e-4-vs-1e-5-vs-1e-6 disagreement between the three designs empirically. Run the same pass twice, with `γ_mig`/`p_umi·δ` enabled and disabled: if the derate makes calibration *worse* (it will, per E4 — it lowers Q on MIGs whose bases are fine), that is the falsification. Cost: one pipeline run plus a pileup script; it produces bench's `quality_calib.tsv` as a by-product, so it is work that has to happen anyway — just do it before committing to the constants rather than after.

**X3 — Empirical nulls from one deep real dataset (no simulation, no truth needed). Falsifies E5, E6, E8 and E10 in a single pass.**
Three permutations over the same `.mig` file for one deep MIGEC/MAGERI run:
1. **Split-half UMI collision** → falsifies the uniform/independent-UMI model behind every birthday formula. Partition reads into two halves by `src_index` parity (or use two sample barcodes from the same run), count UMIs present in both halves whose consensus sequences *differ* → a direct, model-free estimate of `Σ_u p_u²`. Compare against `4^{−L}`, checkout's `2^{−Σ H_p}`, and E5's `Π_j m_j`. If the measured value exceeds `Π_j m_j` materially, position-independence is dead and the collision term must be measured per run, not modelled.
2. **Size-preserving UMI shuffle** → calibrates `D1_obs`. Reassign the observed UMI multiset to the observed MIG-size multiset at random (preserving both marginals, destroying only the parent–child structure), recount `D1_obs`. That empirical null replaces `C(n,2)·P1` and directly exposes the sibling-pair double count of E8: if the shuffled `D1_obs` already accounts for most of the observed excess, the error-rate estimator is measuring nothing.
3. **Within-MIG read re-partition** → calibrates the split test. For MIGs with `n ≥ 20` from a clonal template, randomly split the reads in half and run the full §2 split-acceptance test on each half-pair. Any accepted split is a **false positive by construction**. Measure the false-split rate as a function of `d_split`, `ΔBIC`, `--min-minor-reads` and `--min-minor-frac`, and set them from that curve. This is the only honest way to set `d_split` and it costs one pass; the Poisson derivation in §2(a) can then be kept as documentation rather than as the source of the threshold.

All three are single passes over one already-available dataset, need no ground truth and no cluster time, and each one kills a formula that is currently load-bearing across multiple designs.