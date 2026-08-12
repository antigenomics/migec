# Design record — 2026-08-12

The v2 design was produced by six independent subsystem designs followed by two adversarial
critiques. These are the raw outputs, kept because the *rejected* options and the corrected
mathematics are worth more than the summary: `CLAUDE.md`'s rules are terse, and this is where each
one's justification lives.

| File | Content |
|---|---|
| `design-checkout.md` | pattern grammar, bit-parallel matching, LLR acceptance, whitelists, sample sheets, placement suggestion |
| `design-refine.md` | barcode error model, birthday prior, correction algorithm, MIG-size model, QC statistics |
| `design-io-interop.md` | record format, partitioning, grouping iterator, FASTQ IO, downstream interop (verified against bwa/minimap2/fgbio/arda) |
| `design-assemble.md` | consensus algorithm, multi-consensus splitting, quality recalibration, doublets, contigs |
| `design-repo.md` | archive procedure, file tree, CMake, pyproject, CI, docs, notebooks |
| `design-benchmarks.md` | dataset inventory, HF layout, comparison matrix, metrics, acceptance gates |
| `review-scope.md` | the cut list (25 components removed or deferred), conflict resolutions, the minimal command surface, milestones |
| `review-algorithms.md` | 25 numbered errors in the designs above, each with the corrected formula; missing coverage; the three falsifying experiments |

## How to read this

`review-algorithms.md` supersedes the six designs wherever they disagree — it is a correction
pass over them, and several of its findings invalidate formulas that appear in the design files
unchanged. In particular the designs' versions of the following are **wrong** and must not be
copied out of them:

- the per-base consensus error formula in `design-refine.md` §3.2 (a Dirichlet posterior mean of
  the minor-allele frequency, not an error probability — returns Q9 for ten agreeing Q30 reads)
- the MIG-identity quality derate in both `design-refine.md` and `design-assemble.md` (applies the
  penalty in inverse proportion to the risk it claims to cover)
- effective barcode space via Shannon entropy in `design-checkout.md` §6 and `design-refine.md` §4
  (the birthday functional is Rényi-2 collision entropy; Shannon underestimates collisions)
- the hash-bucketing scheme in `design-io-interop.md` (makes correction impossible to apply)
- the `d_split` Poisson test in `design-assemble.md` §2 (tests polymorphism where linkage is what
  discriminates, and the BIC penalty is weak enough that a single read with five errors splits)

`review-scope.md` supersedes all of them on *what gets built*: thirteen proposed CLI commands
became five, three competing on-disk formats became one.
