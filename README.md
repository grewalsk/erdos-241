# Erdős Problem #241 — maximum size of B₃ sets

**Status snapshot: 2026-06-12 (all workstreams complete; see `results.md`)**

> Quick map: [`results.md`](results.md) (labeled outcomes) · [`writeup/literature_memo.md`](writeup/literature_memo.md) ·
> [`writeup/modular_theorem.md`](writeup/modular_theorem.md) (THEOREM candidate) ·
> [`writeup/upper_bound_program.md`](writeup/upper_bound_program.md) (Green + certified kernel) ·
> [`writeup/assessment.md`](writeup/assessment.md) (final assessment).
>
> Headlines: (1) exact f(N) to 159 = OEIS A387704, b-file extended by 9 terms; (2) modular THEOREM
> candidate k³≤2m+3k² (adversarially verified); (3) **certified** improvement of Green's constant
> 1.51825→1.51587 N^{1/3} by optimizing his kernel; (4) Bose–Chowla shown asymptotically locally
> maximal — no family beats ratio 1.

## The problem

A finite set A of positive integers is a **B₃ set** if every integer has at most one
representation as a+b+c with a ≤ b ≤ c, a,b,c ∈ A (repetition allowed; representations
are multisets). Let f(N) = max size of a B₃ subset of {1,…,N}.

**Erdős #241 (Bose):** Is f(N) ~ N^{1/3}?

Known: Bose–Chowla (1962) gives f(N) ≥ (1+o(1))N^{1/3}; the best upper bound is
Green (2001), commonly quoted as (3.5N)^{1/3}(1+o(1)). Neither constant has moved since 2001.

## Progress so far

### Infrastructure (done)
- [`code/verify.py`](code/verify.py) — **standalone brute-force B₃ verifier** (integer and
  modular), independent of every solver/construction here. Everything below passes through it.
- [`code/brute_small.py`](code/brute_small.py) — verifier-driven ground-truth search
  (no clever logic), used to validate the fast solvers.

### Workstream A — exact values (running)
- [`code/fsolver.c`](code/fsolver.c) — exact f(N) via branch-and-bound (pair/triple-sum
  bitsets, Sidon invariant — lossless since B₃ ⇒ Sidon — and f-table suffix pruning).
  Matches ground truth for N ≤ 32; all witnesses re-verified by `verify.py`.
  **Currently past N ≈ 119**: f(119) = 7, ratio f/N^{1/3} ≈ 1.42 and slowly declining.
  Output: [`data/fN.csv`](data/fN.csv). OEIS A387704 cross-check pending (literature agents).
- [`code/gsolver.c`](code/gsolver.c) — exact g(m) for B₃ sets in ℤ/mℤ (honest pruning only).
  Matches ground truth for m ≤ 30. **Currently past m ≈ 199** (g = 6).
  Jump points (first m admitting size k): k=2→m=4, 3→13, 4→30, 5→65, 6→117.
  Output: [`data/gm.csv`](data/gm.csv).

### Workstream B — lower bound search (mostly done)
- [`code/bose_chowla.py`](code/bose_chowla.py) — Bose–Chowla h=3 over GF(q³) for prime
  powers q (hand-rolled field arithmetic; Rabin irreducibility test — the gcd step matters
  for composite extension degree, bug found & fixed at q=25). **Verified B₃ (mod q³−1 and
  as integers) for all q ∈ {2,…,64}** by the standalone verifier.
- **Empirical finding (augmentation):** the number of single elements of [1, q³−1] addable
  to the Bose–Chowla set collapses as q grows — 26 addable at q=4, 7 at q=8, 1 at q=11,
  **0 for every q ≥ 13 tested (13…64)**. Bose–Chowla sets become locally maximal; greedy
  augmentation ratios decline toward 1. Evidence *for* the conjecture from the disproof route.
  Data: [`data/bose_chowla.csv`](data/bose_chowla.csv), [`data/bose_chowla2.csv`](data/bose_chowla2.csv).

### Workstream C — upper bound program (centerpiece so far)
- **Modular theorem candidate** ([`code/moment_bound.py`](code/moment_bound.py), sympy-exact):
  if A ⊆ ℤ/mℤ is B₃ with |A| = k ≥ 2, then

  **m ≥ (k⁴ − 3k³ + 5k² − 4k)/(2k − 3), hence k³ ≤ 2m + 3k².**

  This improves the counting bound k³ ≤ 6m (constant 6 → 2 asymptotically, i.e.
  (6m)^{1/3} → (2m)^{1/3} ≈ 1.26 m^{1/3}, vs. the Bose–Chowla lower bound m^{1/3}).
  Chain: B₃ ⇒ Sidon (padding) ⇒ E₂ = 2k²−k exactly; B₃ ⇒ E₃ = 6k³−9k²+4k exactly;
  Parseval converts both to 4th/6th Fourier moments; Cauchy–Schwarz on nontrivial
  characters forces the 6th moment up, contradicting the E₃ cap unless k³ ≲ 2m.
  Validation so far: algebra machine-checked; E₂/E₃ exactness + Parseval confirmed
  numerically on Bose–Chowla sets; bound tight at k=2 (m=4, {0,1}); **consistent with all
  exact g(m) data** (actual extremal m exceeds the bound's requirement by 15–30%);
  within-method optimality proven (min E[y³] given E[y]=1, E[y²]=2, y≥0 is exactly 4,
  so constant 2 is the best this constraint set yields).
  **Pending:** novelty check against literature (agents fetched Green's paper + O'Bryant DS11,
  see `lit/`), adversarial proof-verification pass, full writeup.

### Literature pass (running)
Background multi-agent sweep has fetched primary sources into [`lit/`](lit/):
Green 2001 PDF, O'Bryant's DS11 dynamic survey, erdosproblems.com/241 page + forum/history,
OEIS A387704 raw. Verification memo to be written when the sweep completes.

## Still to do
1. Literature memo (incl. OEIS A387704 convention check, modular-bound novelty verdict).
2. Adversarial verification workflow on the modular theorem proof.
3. Green h=3 proof reconstruction; identify the slack step.
4. μ₃ = inf‖f∗f∗f‖∞ transfer derivation (exact factor) + numerics — expected negative
   result for beating 3.5N, to be documented with diagnosis.
5. Cilleruelo-identity / low-frequency-sector LP: ceiling of that scheme vs 3.5N.
6. Structured lower-bound family search (beyond augmentation).
7. `results.md` with THEOREM / COMPUTATION / HEURISTIC labels, final assessment.

## Layout

```
code/      verifier first, then solvers and constructions
data/      CSVs from the exact solvers and sweeps
lit/       fetched primary sources + literature notes
writeup/   memo, theorem writeup, results (in progress)
```
