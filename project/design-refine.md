# AREA: barcode error correction, statistics, `refine`

Scope note: this stage consumes the sorted binary record store produced by checkout+sort (records keyed `(sample_id, cell_packed, umi_packed)`) and emits (a) a correction map, (b) QC TSVs, (c) a rewritten/merged record store + FASTQ headers for `assemble`. It never sees read sequences except for the optional overlap in checkout.

---

## 1. The correction model

### 1.1 State

Per **correction scope** (see §1.7) we have `n` distinct barcodes `x_1..x_n` of length `L`, read counts `c_i`, and a per-position mean error probability stored as a phred-like `uint8` (`meanq[i*L + j]`).

### 1.2 The three signals, fused into one likelihood ratio

For an ordered candidate pair (child `x`, parent `y`), `d_H(x,y)=1` at position `j`, parent base `a`, child base `b`:

**(b) Phred → per-event error probability.**

```
eps_j(x|y)  =  kappa * pbar_j(x) / 3
pbar_j(x)   =  (1/c_x) * sum_{r in reads(x)} 10^(-Q_{r,j}/10)      [stored as meanq]
```
`/3` because a miscall to one *specific* alternative base. `kappa` is the calibration factor from (a).

**(a) Error rate from the data — three estimators, one output `kappa`.**

- `eps_phred = mean over all barcode bases of 10^(-Q/10)` (assumes calibrated phred).
- **Neighbour method-of-moments (the primary one, no external truth needed).** Let `D1_obs` = number of unordered distance-1 pairs among the `n` distinct barcodes (counted for free by the same neighbour scan, §2). Let `P1` be the probability that two *independent* barcodes are at distance 1 under the observed composition model (§1.3). Then

```
D1_ind   = (n choose 2) * P1                                  # independent-collision component
D1_err   = max(0, D1_obs - D1_ind)                            # error component
```
  and the expected number of *observed distinct* error children is
```
E[children](eps) = sum_i  3L * (1 - exp(-c_i * eps))
```
  Solve `E[children](eps) = D1_err` for `eps` by bisection on `log eps` over `[1e-6, 1e-1]` (monotone increasing, ~30 iterations, cost `O(n)` per iteration → do it on the size *histogram*, not on `n` items: `O(#distinct sizes)`).
  Output `eps_nbr`, and `kappa = eps_nbr / eps_phred` (clip to `[0.2, 20]`, warn outside).
- **Whitelist estimator** (10x, or any designed barcode set): from the exact-match pass, `eps_hat(j, a->b) = reads at distance 1 from a whitelist entry with substitution (j,a->b) / reads exactly matching`. Gives the full `L x 4 x 4` table for free; use it directly instead of `kappa * phred` when `--whitelist` is given.

**(c) Birthday-paradox prior.** Under the observed composition model `p_j(a)` (the PWM of §4), define per-position agreement probability
```
m_j = sum_a p_j(a)^2                      (= 1/4 when uniform)
pi(x) = prod_j p_j(x_j)                   (= 4^-L when uniform)
P1  = sum_j (1 - m_j) * prod_{k != j} m_k (= 3L / 4^L when uniform)
```
Expected number of *independent-molecule* pairs at Hamming distance 1 among `n` observed barcodes:
```
E_pairs = (n choose 2) * P1  ~=  n^2/2 * 3L/4^L        [uniform case, as stated]
```
Per-UMI expected independent 1-mm neighbours: `lambda_ind = (n-1)*P1`.

### 1.3 Prior odds (presence), likelihood ratio (counts), posterior

Two *presence* priors:
```
A_err(x|y) = 1 - exp(-c_y * eps_j(x|y))      ~= c_y * eps_j    # >=1 erroneous read of y lands exactly on x
A_ind(x)   = 1 - exp(-n * pi(x))             ~= n * pi(x)      # >=1 independent molecule carries x
```
Two *count* likelihoods, both zero-truncated:
```
Lerr(c_x | y) = Poisson(c_x ; c_y * eps_j) / (1 - exp(-c_y*eps_j))
Lind(c_x)     = f_real(c_x)                                    # the NB MIG-size pmf of §5
```

**Single-parent posterior odds (the formula to implement, in logs):**
```
log O = log A_err(x|y) + log Lerr(c_x|y) - log A_ind(x) - log Lind(c_x)
      = log c_y + log eps_j + log P(c_x; c_y*eps_j) - log n - log pi(x) - log f_real(c_x)
```
`P(error) = sigmoid(log O)`. Merge iff `P(error) >= tau` (default `tau = 0.95`).

### 1.4 Many candidate parents

`x` has up to `3L` 1-mm neighbours; keep those with `(c_y, packed_y) > (c_x, packed_x)` (strict order → the merge graph is a forest, no cycles, no tie deadlock). Mixture responsibilities:

```
w_k = A_err(x|y_k) * Lerr(c_x|y_k)     for each candidate parent y_k
w_0 = A_ind(x)     * Lind(c_x)

P(error | .) = (sum_k w_k) / (w_0 + sum_k w_k)
parent(x)    = argmax_k w_k
```
Summing over parents is the correct behaviour: several plausible parents *increases* the error posterior. Report `n_candidate_parents` and the margin `w_(1)/w_(2)` in the correction map; if `P(error) >= tau` but `w_(1)/w_(2) < 3`, tag `decision = merge_ambiguous` (still merged — reads are not lost — but flagged).

Transitivity (distance-2 chains) comes free from union-find over accepted merges; do **not** run a distance-2 search.

### 1.5 Recovering and correcting MIGEC's gate

MIGEC enabled collision filtering only when `n < 0.05 * 4^(L-1)`. Substituting into `lambda_ind = n * 3L/4^L`:
```
lambda_ind* = 0.05 * 4^(L-1) * 3L/4^L = 0.0375 * L
```
i.e. MIGEC's gate is *"turn correction off once the expected number of independent 1-mm neighbours per UMI exceeds 0.0375·L"* — 0.30 at L=8, **0.45 at L=12**, 0.60 at L=16. The `L`-dependence is an artefact of writing the gate in terms of `4^(L-1)` and has no statistical meaning.

**Corrected general rule:** there is no global gate. `A_ind(x) = n*pi(x)` enters every pair decision, so saturation suppresses merging automatically and *per-barcode* (a barcode in a high-probability composition region is protected more than one in a rare region — MIGEC could not express this). If a single reported number is wanted, publish
```
lambda_ind = (n-1) * P1        # composition-aware, not 4^-L
```
and emit `WARN_SATURATED` when `lambda_ind > 0.05`, which at L=12 uniform is `n > 0.05*4^12/(3*12) = 23,301` — 9x stricter than MIGEC's implied 209,715. That is the substantive change: MIGEC was correcting far into the saturated regime.

### 1.6 Whitelist mode (designed barcode sets, 10x CB)

Same equation, different `A_ind`: the "independent" hypothesis becomes "`x` is itself a real whitelist barcode", and the parent prior is the whitelist abundance:
```
P(true = y | obs = x, Q) ∝ (n_y + 1) * prod_{j: x_j != y_j} 10^(-Q_j/10)/3
```
`n_y` = reads exactly matching `y` in the exact-match pass. Correct if `max_y P >= 0.975` (Cell Ranger's constant; expose as `--whitelist-tau`). Barcodes not in and not within 1 of the whitelist are routed to `undef-cb` and counted, never silently dropped.

### 1.7 Correction scope

`--umi-scope {sample,cell}`; default `cell` when cell barcodes are present, else `sample`. This sets `n` in the birthday term, and it matters enormously: per-cell `n` is `10^3–10^4` (correction is aggressive and correct), per-sample `n` is `10^6–10^7` (correction is conservative). Cell barcodes are always corrected at `--cell-scope sample`.

---

## 2. The algorithm at scale

### 2.1 Recommended: 2-bit hash + neighbour enumeration (and I would cut seqtree here)

`L <= 32` barcodes pack into a `uint64`. Build an open-addressing hash `packed -> index` (2x load factor). For each barcode `i`, enumerate its `3L` substitution neighbours by XOR on 2-bit lanes and probe.

- **Complexity:** `O(n * 3L)` probes; 360M probes at `n=1e7, L=12`; embarrassingly parallel over disjoint index ranges, no synchronisation (writes go to per-thread `vector<Merge>`). Expect a few seconds on 8 threads.
- **Memory at 1M / 10M distinct:** hash `16*n` bytes (160 MB at 10M) + `packed` 8B + `count` 4B + `meanq` `L`B per barcode → **~28 MB at 1M, ~400 MB at 10M** (L=12). No adjacency list is ever materialised: the decision is made inside the per-query loop and only accepted merges are emitted.
- Same loop counts `D1_obs` for the §1.2 estimator.

**Why not seqtree here:** for fixed-length barcodes with `max_substitutions=1`, brute neighbour enumeration is exact, needs no index build, uses ~4x less memory than the trie, and is faster. Using seqtree would be reuse-for-its-own-sake. Flagging this as the single biggest simplification available.

### 2.2 Where seqtree *is* used

`seqtree::Index::build(refs, Alphabet::Nucleotide)` + `search_batch(queries, params, threads)` (which releases the GIL and is lock-free on the shared index) for:
1. **Whitelist search** — build once over the 737K/3M 10x whitelist, `Index::save()`/`load()` it to `~/.cache/migec/`, `SearchParams{engine=SeqTm, max_substitutions=1, mode=TopHit, max_hits=8}`. Trie ~3M nodes at 737K refs of 16nt → ~60–100 MB; amortised over the whole run.
2. **`--umi-indels`** — `max_insertions=1, max_deletions=1, max_total_edits=1`, variable-length UMIs (Ion Torrent / nanopore / homopolymer-adjacent anchors). Neighbour enumeration cannot do this.
3. Checkout's primer/adapter fuzzy matching (other agent's area).

At `n=1e6` L=12 the trie would be ~3.3M nodes (`sum_d 4^d (1-exp(-n/4^d))`) → ~100 MB; at `n=1e7` L=16 → ~30M nodes → ~1 GB. Quoted for the record; not on the hot path.

### 2.3 Threading plan

Shared immutable `BarcodeTable` (SoA, read-only) + shared hash. `std::thread` over `[lo,hi)` index ranges, `py::gil_scoped_release` at the binding boundary, zero synchronisation, per-thread output vectors concatenated at the end. `--threads 0` → `hardware_concurrency()`. The union-find pass over accepted merges is single-threaded (`n <= 1e7` → milliseconds).

### 2.4 Comparison with alternatives

| Method | Rule | Why not chosen |
|---|---|---|
| umi_tools `directional` | edge `y->x` iff `d=1` and `c_y >= 2*c_x - 1` | Special case of ours with a flat count likelihood and a constant prior: no phred, no `n`, no `L`, no composition. Implemented as `--legacy-directional` for benchmarking only. Its `O(k^2)` per-group all-pairs adjacency also does not scale past ~1e4 UMIs per group. |
| Calib | minhash-LSH over (barcode ‖ read prefix), then clustering | Solves a *joint* barcode+sequence clustering problem, so it needs read sequences at correction time and breaks the stage separation (refine must run before assemble, streaming). LSH is approximate with recall knobs; enumeration is exact and knob-free. |
| MIGEC v1 | enumerate `4L` neighbours, child if `c_x < ratio * c_y`, gated by `n < 0.05*4^(L-1)` | We keep the enumeration (it is the right primitive) and replace the fixed ratio + global gate with §1.3. |

---

## 3. The low-coverage retention rule

MIGEC dropped every MIG with `count < 5`. We drop **nothing for count alone**; a barcode leaves the output only by being *absorbed* into a parent (its reads are not lost either way).

### 3.1 Retention rule (evaluated per barcode `x` after §1.4)

```
1. P(error) >= tau            -> MERGE into parent(x). Reads join the parent MIG.
2. tau_lo <= P(error) < tau   -> MERGE into parent(x), decision = "merge_weak",
                                 and the parent MIG's confidence is multiplied by
                                 (1 - P(error)) of the absorbed child... no: the parent
                                 keeps gamma; the child is recorded in the map.
3. P(error) <  tau_lo  and  c_x >= m_hi   -> KEEP as MIG, gamma = 1        (status "hi")
4. P(error) <  tau_lo  and  c_x <  m_hi   -> KEEP as MIG, gamma < 1        (status "lq")
5. c_x < m_lo (default m_lo = 1)          -> DROP  [unreachable at the default]
```
`tau = 0.95`, `tau_lo = 0.5`, `m_hi` from the FDR fit of §5 (typically 3–8), `m_lo = 1`.
The user's case — 3–5 reads, no plausible parent — lands on rule 4: **kept, tagged `lq`, derated.**

### 3.2 The derate, and how it propagates

The per-MIG confidence is exactly the mixture posterior of §5 combined with the neighbourhood posterior:

```
gamma_mig(x) = (1 - P(error | .)) * P_real(c_x)

P_real(c) = (1-w) * f_real(c) / [ (1-w)*f_real(c) + w*f_err(c) ]      # from the §5 mixture
```
Both factors are already computed; there is no new tunable constant. Clip `gamma_mig` to `[1e-5, 1]`.

Per-base consensus error, from the count-weighted PWM of the MIG with a Dirichlet-`alpha` prior (`alpha = 0.5`), then the MIG confidence, then the chemistry floor:

```
p_cons(i)  = 1 - (n_maj(i) + alpha) / (c_x + 4*alpha)
p_final(i) = 1 - (1 - p_cons(i)) * gamma_mig
Q(i)       = min( Q_max, round(-10 * log10(p_final(i))) ),  Q_max = -10*log10(eps_RT)
```
`--rt-error 1e-5` → `Q_max = 50` (`1e-6` → 60). MIGEC clamped to 40 and used a linear rescale `max(2, (maxFreq/count - 0.25)/0.75 * 40)`; the Dirichlet form is the same shape but calibrated and non-saturating at `c=1`.

Propagation: `gamma_mig` is written into the consensus FASTQ header as `cq:i:<round(-10*log10(1-gamma))>` so a downstream tool can filter MIG-level confidence without recomputing anything, and it is already folded into every base Q.

### 3.3 What changes in the output vs MIGEC `--min-count 5`

Every MIG with `1 <= c < 5` and no plausible parent now appears, at `Q ~ 15–30` instead of absent. Quantify in `qc_summary.tsv` with `reads_in_lq` / `reads_in_hi` / `reads_merged` (typically 20–50% of distinct UMIs and 5–15% of reads move from "discarded" to "kept, derated"). `--min-reads 5 --no-derate` reproduces MIGEC exactly for regression testing.

---

## 4. UMI quality / non-uniformity

Computed in one pass over `BarcodeTable` (weighted by reads and, separately, by distinct UMIs — report both; distinct-UMI weighting is the one that matters for `pi(x)`).

| Statistic | Formula | Actionable? | Warning |
|---|---|---|---|
| PWM | `p_j(a)` over distinct UMIs | feeds `pi(x)`, `P1` | — |
| Per-position KL | `D_j = sum_a p_j(a) * log2(4*p_j(a))` bits, in `[0,2]` | yes | `WARN_POS_BIAS` if any `D_j > 0.1` bit |
| Effective length | `L_eff = (2L - sum_j D_j)/2` | **yes — the key one** | `WARN_LOW_DIVERSITY` if `L_eff < 0.9*L` |
| Saturation | `s = n / 4^L_eff` | **yes** | `WARN_SATURATED` if `lambda_ind = (n-1)*P1 > 0.05` |
| Obs/exp 1-mm ratio | `rho = D1_obs / ((n choose 2) * P1)` | **yes** | `rho ~= 1` → `WARN_NO_ERROR_SIGNAL` (correction will do nothing); `rho < 1` → `WARN_DESIGNED_SET` (switch to `--whitelist`) |
| Per-position mean Q, N-rate | from `meanq`, N counter | **yes** | `WARN_UMI_WINDOW` if mean Q < 20 or N-rate > 1% at any position → the UMI window is probably misplaced; points the user at checkout's `--suggest` |
| Total entropy vs `2L` bits | identical to `sum_j D_j` | — | do not report twice |
| GC / AT skew | `(G+C)/tot - 0.5` | no | ship as PWM columns, not a metric |
| Multinomial chi-square GOF | `chi2_j = c * sum_a (o-e)^2/e`, 3 df | **no — cut.** Astronomically significant at any real `n`; effect size (`D_j` in bits) is what matters | — |

---

## 5. MIG-size distribution and threshold selection

### 5.1 UMI barcodes: two-component zero-truncated mixture

Fit on the **size histogram** (not the barcode list): `O(#distinct sizes)`, trivial cost.

```
f(c) = w * f_err(c) + (1-w) * f_real(c),    c >= 1

f_real(c) = ZTNB(c; r, mu)   = NB(c; r, mu) / (1 - NB(0; r, mu))    # Poisson-Gamma, PCR overdispersion
f_err(c)  = ZTNB(c; r, mu_err),  mu_err := eps_hat * mu             # CONSTRAINED, not free
```
Constraining `mu_err` to the estimated error rate times the real mean removes the identifiability problem, leaving **three free parameters** `(w, r, mu)`. Fit by EM (E-step: responsibilities per size bin; M-step: weighted method-of-moments for `mu`, 1-D Newton for `r`), `<= 100` iterations, converge on `|dlogL| < 1e-6`.

**Threshold at target FDR:**
```
m_hi = min { m : [ sum_{c>=m} w*f_err(c) ] / [ sum_{c>=m} f(c) ] <= alpha },   alpha = 0.05
```
This replaces MIGEC's `round(2 ** (argmax_{bins 3..10} smoothed_log2_hist / 2))` geometric-midpoint heuristic, which had no error-rate input at all. Per §3, `m_hi` is the *full-confidence* line, not a discard line.

**Over-sequencing call** replaces `sum(bins 0..2) < 2*sum(bins 3..16)` with: over-sequenced iff `mu >= 5 and w < 0.5`. Also report `mean_reads_per_umi = R/n` and `frac_reads_in_migs_with_c>=5`.

`P_real(c)` from the same fit is the derate of §3.2 — one model, two uses.

### 5.2 Cell barcodes: OrdMag, no distribution fitting

Different generative process (no PCR-family structure, 4–5 log dynamic range, ambient plateau). Use Cell Ranger's OrdMag:

```
sort counts descending
t   = quantile_0.99( counts[0 : N_exp] )        # N_exp = --expect-cells, default 3000
thr = max(1, round(t / 10))
is_cell(b) = count(b) >= thr
```
Plus, for the plot only, the **inflection** (min of `d log(count) / d log(rank)`) and the **knee** (max curvature of the same log-log curve, over a monotone-smoothed count vector). Do **not** fit an NB to cell barcodes.

---

## 6. The 10x cell-barcode path

**Implement in v1 (all of it is one small module):**
- Whitelist load (`737K-august-2016.txt`, `3M-february-2018.txt`, gz), exact-match pass, then 1-mm correction with the §1.6 posterior — same code path as UMI correction, different prior source.
- Barcode-rank curve TSV + knee + inflection + OrdMag call → `is_cell` boolean.
- Ambient profile TSV (all below-threshold barcodes with their read/UMI counts).

**Delegate / cut:**
- **EmptyDrops** — cut. It needs the gene x cell count matrix, which we deliberately never build (the pipeline ends at consensus FASTQ). We emit `ambient.tsv`; DropletUtils / Cell Ranger does the Dirichlet-multinomial LR + Monte-Carlo p-values downstream.
- **SoupX / CellBender ambient removal** — entirely out of scope.
- **Cellplex / hashtag demultiplexing** — out of scope; that's a sample barcode read from a separate library.
- **Doublet calling** — not here. Refine contributes one input to the consensus stage's doublet step: `<S>.cb_umi_graph.tsv` (`cell, n_umis, n_umis_multigroup, frac_multigroup`), consumed by the multi-consensus-per-UMI logic. Gate behind `--emit-cb-umi-graph`.

---

## 7. Output artifacts

All plots are rendered by `python/migec/plots.py` (matplotlib) **from the committed TSVs**; no plotting code in C++, and `migec plot DIR` regenerates every figure from artifacts alone.

| File | Columns | Plots generated from it |
|---|---|---|
| `<S>.mig_size_hist.tsv` | `reads_per_umi, n_umis, n_reads, cum_umis, cum_reads, f_err_fit, f_real_fit, p_real, fdr_at_m` | `mig_size_distribution.png` (log-log hist + two fitted components + `m_hi` line), `mig_saturation.png` (cum_reads vs cum_umis) |
| `<S>.correction_map.tsv.gz` | `child_umi, child_cell, child_reads, parent_umi, parent_cell, parent_reads, mm_pos, child_base, parent_base, child_q_at_pos, log_A_err, log_L_err, log_A_ind, log_L_ind, log_odds, posterior, n_candidate_parents, margin, decision` (`merge`/`merge_weak`/`merge_ambiguous`/`keep`) | `correction_by_position.png` (mm_pos x decision), `correction_posterior_hist.png`, `correction_count_ratio.png` (child_reads vs parent_reads, coloured by decision) |
| `<S>.umi_pwm.tsv` | `pos, count_a, count_c, count_g, count_t, count_n, freq_a, freq_c, freq_g, freq_t, kl_bits, mean_q, frac_n` | `umi_pwm.png` (stacked bars + KL line + mean-Q line) |
| `<S>.qc_summary.tsv` (one row) | `sample, total_reads, distinct_umis_raw, distinct_umis_corrected, umi_len, l_eff, saturation, lambda_ind, d1_obs, d1_expected_ind, rho, err_rate_phred, err_rate_neighbour, kappa, nb_r, nb_mu, mixture_w, m_hi_fdr, reads_hi, reads_lq, reads_merged, umis_hi, umis_lq, umis_merged, mean_reads_per_umi, warnings` | `qc_summary` is also the multi-sample comparison table |
| `<S>.barcode_rank.tsv` (CB only) | `rank, cell, reads, umis, is_cell, ordmag_threshold, knee_rank, inflection_rank` | `barcode_rank.png` |
| `<S>.ambient.tsv` (CB only, opt) | `cell, reads, umis` | — (input to DropletUtils) |
| `<S>.barcode_table.tsv.gz` (opt, default **off**) | `sample, cell, umi, reads, reads_absorbed, n_children, parent, status, posterior_error, p_real, gamma_mig, mean_q, min_q` | — (can be 10M rows) |
| `refine.filelist.txt` | manifest mirroring checkout's | — |

**Write-back.** The correction map is applied **during the merge phase of the external sort** — one pass, no re-sort: the merge reads the `(old_key -> new_key)` map (a `uint64->uint64` hash, `16*|merges|` bytes) and rewrites `cb_packed`/`umi_packed` as records stream out, so the store emerges already grouped by corrected key.

**FASTQ headers** — SAM-tag style so tags survive into BAM via `minimap2 -y` / `bwa -C` and are read directly by arda and `umi_tools --extract-umi-method=tag`:

```
@<readname> BC:Z:<sample_bc> CB:Z:<corrected_cb> CR:Z:<raw_cb> CY:Z:<cb_qual> UB:Z:<corrected_umi> UR:Z:<raw_umi> UY:Z:<umi_qual> cq:i:<mig_confidence_phred>
```
Consensus FASTQ read name: `<sample>:<cell>:<umi>:<group_idx>` with `RD:i:<reads_in_mig> cq:i:<...>`. `--header-style migec` restores the legacy ` UMI:<seq>:<qual>`.

---

## 8. C++ files, classes, and CLI

```
include/migec/
  barcode.hpp        struct Barcode{uint64 packed; uint8 len;}; pack2bit/unpack2bit;
                     sub_neighbours(uint64,L,cb)  -- 3L XOR-lane enumeration
  barcode_table.hpp  class BarcodeTable  (SoA: packed[], reads[], meanq[n*L], flags[])
                     built streaming from the sorted store; count(i), mean_err(i,j)
  composition.hpp    class CompositionModel  from the PWM:
                     log_pi(uint64), p_dist1(), kl_bits(pos), l_eff(), lambda_ind(n)
  error_model.hpp    struct ErrorModel{double kappa; array<double,64> q2p; vector<double> per_pos;}
                     eps(q,pos); estimate_from_phred / estimate_from_neighbours (MoM bisection)
                     / estimate_from_whitelist
  mig_size_model.hpp struct MigSizeModel{w,r,mu,mu_err; log_f_real, log_f_err, p_real,
                     threshold_at_fdr(alpha)};  fit(hist, eps)  -- EM
  corrector.hpp      struct CorrectionParams{tau,tau_lo,indels,threads,legacy_directional};
                     struct Merge{uint32 child,parent; float posterior; uint8 pos; uint8 n_cand;};
                     class BarcodeCorrector::run(table, err, size_model, comp, params)
                       -> vector<Merge>   (hash path; seqtree path when indels/variable length)
  whitelist.hpp      class Whitelist  (gz load, exact hash, 1-mm enumeration,
                     correct(obs,quals)->optional<uint32> via the CR posterior; seqtree Index
                     save/load cache for the 3M list)
  cell_calling.hpp   ordmag_threshold(sorted_counts, n_expected);
                     knee_and_inflection(sorted_counts)
  stats_report.hpp   write_mig_size_hist / write_correction_map / write_umi_pwm /
                     write_qc_summary / write_barcode_rank / write_ambient
src/  matching .cpp;  src/_bindings.cpp exposes refine(), MigSizeModel, ErrorModel,
                      CompositionModel, ordmag_threshold  (for tests + plot driving)
```
Union-find is ~20 lines inside `corrector.cpp`, not its own header.

### CLI

```
migec refine INPUT... -o DIR
  --umi-len INT                     (auto-detected)      --cell-len INT
  --umi-scope {sample,cell}         (default: cell if CB present else sample)
  --whitelist PATH | 10x-v2 | 10x-v3       --whitelist-tau 0.975
  --error-rate auto|FLOAT           --error-rate-from {neighbour,phred,whitelist}
  --qual-calibration / --no-qual-calibration      # kappa
  --tau 0.95  --tau-low 0.5
  --fdr 0.05                        # sets m_hi
  --min-reads 1                     # hard floor; 5 + --no-derate == MIGEC legacy
  --no-derate                       # hard cut at m_hi instead of gamma
  --rt-error 1e-5                   # Q cap = -10*log10(rt_error)
  --legacy-directional              # umi_tools rule, benchmarking only
  --umi-indels                      # seqtree path, default off
  --expect-cells 3000  --cell-caller {ordmag,knee,none}
  --emit-correction-map / --no-emit-correction-map   (default on)
  --emit-barcode-table              (default off)     --emit-ambient  --emit-cb-umi-graph
  --header-style {sam,migec}
  --threads 0  --tmp-dir PATH

migec stats INPUT... -o DIR         # all estimators, NO rewriting; TSVs + suggested thresholds
  --umi-len --sample --plots/--no-plots --threads

migec plot DIR                      # regenerate every figure from committed TSVs
```

---

## What I would cut from v1

1. **seqtree on the UMI 1-substitution path.** Use the 2-bit hash + `3L` enumeration (§2.1): exact, faster, ~4x less memory, no index build. Keep seqtree for whitelists, indel mode, and checkout primers. Biggest simplification available.
2. **EmptyDrops.** Emit `ambient.tsv`, delegate to DropletUtils. It needs a count matrix we never build.
3. **Chi-square GOF and GC-skew as metrics.** Redundant with per-position KL; ship them as PWM columns only.
4. **Distance-2 direct search.** Transitive union-find over accepted 1-edit merges covers it.
5. **Indel error-rate estimation.** `--umi-indels` uses a flat `--indel-rate 1e-4`; don't fit it.
6. **Iterative correction rounds.** One pass with raw counts. Note it in ROADMAP; add `--rounds` only if benchmarks show a gap.
7. **`--emit-barcode-table` default on.** Off by default (10M rows).
8. **NB fitting for cell barcodes.** OrdMag only.

## Open question for the architect

`--umi-scope` default materially changes results (per-cell `n ~ 1e3` vs per-sample `n ~ 1e7` in the birthday prior). I default to `cell` when CBs are present, but for MIGEC-style bulk RepSeq there is no CB and the per-sample `n` is large enough that `WARN_SATURATED` will fire on real data — which is correct, and is precisely the regime MIGEC v1 was silently over-correcting in. Confirm that firing that warning (rather than suppressing correction) is the desired behaviour.