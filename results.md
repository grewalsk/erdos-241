# results.md — Erdős Problem #241 (maximum size of B₃ sets)

Every item is tagged **THEOREM** (complete proof, checkable), **COMPUTATION** (code +
independent verification), or **HEURISTIC** (flagged; never a premise for a THEOREM).
Convention throughout: B₃ = every integer has ≤1 representation a+b+c, a≤b≤c, multisets,
repetition allowed; f(N)=max|A|, A⊆{1,…,N}. Matches OEIS A387704 and the Lean
formalization `Erdos241.f N 3` (see `writeup/literature_memo.md`).

Standalone verifier (ground truth for everything): [`code/verify.py`](code/verify.py).

---

## Headline outcomes

1. **COMPUTATION A1** — exact f(N) for N ≤ 159; matches OEIS A387704 on all 151 common
   terms (n=0..150) with zero mismatches, and **extends the b-file to n=159** (a(151..159)=8).
2. **THEOREM candidate (modular)** — for any finite abelian group G, |G|=m, a B₃ set has
   k³ ≤ 2m + 3k² (constant 6→2 vs the counting bound). Self-contained proof; adversarially
   verified. Possibly folklore (novelty not claimed). [`writeup/modular_theorem.md`]
3. **COMPUTATION C1 (certified)** — a rigorous improvement of Green's *actual* upper-bound
   constant: f(N) ≤ 1.51587·N^{1/3}(1+o(1)) vs Green's 1.51825, via certified optimization of
   his kernel. Small (~0.16%), idea-free, does not touch the headline "7/2"; but it is a
   rigorous improvement of the best known constant. [`writeup/upper_bound_program.md` §C3]
4. **COMPUTATION B1** — Bose–Chowla sets verified B₃ for all prime powers q ≤ 64; the count of
   single addable elements **collapses to 0 for every q ≥ 13**: BC sets are locally maximal.

---

## Workstream A — exact values and landscape

**COMPUTATION A1 (f(N)).** [`code/fsolver.c`] branch-and-bound with pair/triple-sum bitsets,
the lossless Sidon invariant (B₃⇒Sidon), and f-table suffix pruning. Cross-checked against the
verifier-driven brute force [`code/brute_small.py`] for all N ≤ 32 (identical), and **every
witness re-verified by `verify.is_b3`**. Reached N=159.
- Data: [`data/fN.csv`](data/fN.csv) = (N, f(N), ratio, seconds, nodes, witness).
- **OEIS A387704 cross-check:** all 151 terms n=0..150 agree; a(151..159)=8 extends the b-file.
- f(N) thresholds (smallest N with a size-k set) and the ratio at each (the per-plateau max):

  | k | min N | f/N^{1/3} | witness |
  |---|---|---|---|
  | 3 | 5 | 1.754 | {1,2,5} |
  | 4 | 12 | 1.747 | {1,2,8,12} |
  | 5 | 24 | 1.733 | {1,2,16,19,24} |
  | 6 | 46 | 1.675 | {1,3,12,27,43,46} |
  | 7 | 83 | 1.605 | {1,2,8,51,60,79,83} |
  | 8 | 130 | 1.579 | {1,3,6,35,75,108,121,130} |

  The threshold ratios **decline monotonically (1.75 → 1.58)** — consistent with f(N)~N^{1/3}
  but the descent is glacial and the ratio is still ≈1.58 ≫ 1 at N≈130. Bose–Chowla anchors
  N=q³−1 (7,26,63,124,342,…) give ratio→1 from above; between anchors the realized ratio sits
  well above 1. **No sign of the ratio plateauing above 1, and no sign of fast convergence.**

**COMPUTATION A2 (g(m), modular).** [`code/gsolver.c`], honest pruning only (counting bound;
*not* our spectral bound, to keep the data an independent test). Cross-checked vs brute force
for m ≤ 30. Reached m=200. Data: [`data/gm.csv`](data/gm.csv).
- g(m) thresholds: k=2→m=4, 3→13, 4→30, 5→65, 6→117.
- **Key observation:** g(m) hugs (2m)^{1/3}, NOT m^{1/3}: at m=117, g=6 gives g/(2m)^{1/3}=0.97
  but g/m^{1/3}=1.23. The densest *small* cyclic B₃ sets approach the **modular-theorem upper
  constant** (2m)^{1/3}, not the Bose–Chowla lower one. Whether limsup g(m)³/m is 1 (BC optimal
  in groups) or up to 2 is the modular shadow of #241; data to m=200 cannot decide (decreasing
  after k=4, but only 6 points). **HEURISTIC.**
- All g(m) data are consistent with the modular THEOREM candidate (no violations): the smallest
  m admitting size k always satisfies m ≥ (k⁴−3k³+5k²−4k)/(2k−3) with 15–30% slack.

---

## Workstream B — lower-bound search

**COMPUTATION B1 (Bose–Chowla + augmentation).** [`code/bose_chowla.py`], GF(q³) discrete-log
construction (Rabin irreducibility test, fixed for composite extension degree — bug caught at
q=25 via the standalone verifier). **All q ∈ {2,…,64} verified B₃ both mod q³−1 and over ℤ.**
Data: [`data/bose_chowla_all.csv`](data/bose_chowla_all.csv).
- Baseline ratio |A|/N^{1/3} at N=q³−1: → 1 from above (1.046 at q=2, 1.0007 at q=8, 1.00001 at q=64).
- **Addable-element collapse (FFT-exact count of single integers in [1,q³−1] addable keeping B₃):**

  | q | 2 | 3 | 4 | 5 | 7 | 8 | 9 | 11 | 13 | ≥13 |
  |---|---|---|---|---|---|---|---|---|---|---|
  | addable | 3 | 12 | 26 | 18 | 14 | 7 | 7 | 1 | 0 | **0** |

  For every q ≥ 13 tested, **zero** elements are addable: Bose–Chowla sets are locally maximal in
  their own interval. Greedy augmentation ratios decline toward 1. This is evidence **for** the
  conjecture from the disproof route — the natural dense family has no slack to exploit.

**COMPUTATION B2 (structured family search).** Workflow `wf_caeb452f-833`, five families, all
sets verified by the standalone checker. Synthesis: [`lit/lower_bound_search.md`]. **No family
sustains ratio 1+c over ℤ.** Findings:
- **Greedy** = OEIS A051912; ratio decays as ≈N^{−0.106} (so |A|~N^{0.228}, not N^{1/3}),
  crossing below 1 near N≈8000. The thin baseline, not a dense family.
- **Unions of dilated/translated BC:** pure translates A₀∪(A₀+T) are **never** B₃ for |A₀|≥3
  (explicit collision a_i+a_j+(a_k+T) = a_i+a_k+(a_j+T)); dilated unions need D growing with q,
  forcing ratio ≈2/D^{1/3}→0. Best (q=4) ratio 1.17 is an isolated small-N artifact.
- **Depth local search within [1,q³−1]:** a depth-1 *swap* (remove one BC element, refill)
  **does beat size q at small q** — e.g. size 14 in [1,2196] at q=13, where single-element
  augmentation is already locally maximal — but the surplus over q shrinks +2→+1→**0 for q≥16**.
  Refines B1: BC sets are single-element-maximal from q≥13 but only fully swap-maximal from q≥16;
  ratio still → 1.
- **Modular lift:** modular sets sustain ratio ≈2^{1/3}=1.26 in the torus (free wraparound), but
  lifting to ℤ and paying for max(A) erases the bonus — they tie or lose to f(N) at equal N.
- **Digit/algebraic:** cubes/squares are B₃ but too sparse (ratio→0); Golomb rulers are B₂ not B₃;
  **B₃ sets do not compose** — A+s·B fails B₃ for every spacing s (the triple-sum decouples the
  a- and b-sums, losing the pairing). This is the mechanistic reason the only dense construction
  must be algebraic (Bose–Chowla encodes the pairing into a single field element θ+c).

**Verdict (HEURISTIC, strongly supported):** no construction beats (1+o(1))N^{1/3}; the Bose–Chowla
constant 1 looks genuinely hard to improve. Its (asymptotic) local maximality, the non-composition
of B₃ sets, and the non-transfer of the modular density bonus are three independent pieces of
concrete evidence that the lower bound is tight in shape — i.e. evidence *for* the conjecture.

---

## Workstream C — upper-bound program

**Reconstruction (COMPUTATION C0).** Green's record (7/2)^{1/3}N^{1/3} confirmed from the PRIMARY
source (his Theorem 17) and O'Bryant DS11. His constant comes ENTIRELY from a kernel p he calls
"a simple function that gives a good bound," giving 4/(1+γ)=3.4997 (rounded to 7/2). `code/green_gamma.py`
reproduces his γ to 5 digits. The implicit optimization is min J(p)=Σ_{r≥1}|p̃(πr)|^{4/3}, ∫p=2.

**COMPUTATION C1 (certified kernel improvement).** [`code/green_optimize.py`, `code/green_certify2.py`]
- Exact step-function J via FFT + Hurwitz-zeta (no truncation), validated against Green's quartic.
- Convex minimization (analytic gradient) + Richardson: inf J ≈ 2.3751 (Green 2.4095).
- Certified explicit smooth positive kernel (degree-12, p(0)=p(1)=0): J ∈ [2.37996, 2.37996]
  with a rigorous analytic tail (<9×10⁻⁹), independently cross-checked by direct quadrature.
- ⇒ **f(N) ≤ 1.51587·N^{1/3}(1+o(1))**, below Green's 1.51825 (cubed: 3.4832 vs 3.4997).
- Honest scope: ~0.16% in the coefficient; no new idea (completes Green's own deferred
  optimization); rests on Green's published reduction; finite sum is float64 (margin 0.002 ≫
  error). A genuine improvement of the *known constant*; the asymptotic question is untouched.

**THEOREM candidate C2 (modular).** [`writeup/modular_theorem.md`] For any finite abelian group
G, |G|=m, and B₃ set A with k=|A|≥2:
> m ≥ (k⁴−3k³+5k²−4k)/(2k−3),  hence  **k³ ≤ 2m + 3k²**,  i.e. k ≤ (2m)^{1/3}(1+o(1)).
Proof: B₃⇒Sidon ⇒ E₂=2k²−k exactly; B₃ ⇒ E₃=6k³−9k²+4k exactly; Parseval; Cauchy–Schwarz on the
nontrivial spectrum. Improves the counting constant 6→2. Tight at k=2 ({0,1}⊆ℤ/4ℤ). Within-method
optimal (the moment problem min E[y³] s.t. E[y]=1,E[y²]=2,y≥0 equals 4). Verified: sympy-exact
algebra (`code/moment_bound.py`); E₂/E₃/Parseval confirmed on Bose–Chowla sets; consistent with
all exact g(m) data; **adversarial workflow** `wf_671d547a-fad` — 3/6 independent skeptics returned
"sound", 0 fatal/gap, 2 cosmetics (fixed); remaining lenses (counterexample search, non-cyclic
Fourier) covered by overlap. **Novelty not claimed** (likely folklore; h=2 shadow is classical).

**COMPUTATION C3 (μ₃ sup-norm relaxation — negative result).** [`code/mu3.py`]
μ₃=inf‖f∗f∗f‖∞ ≤ 3/4 (uniform), optimized ≤ 0.588. The provable transfer k³ ≤ (C/μ₃)N gives
constant ≥ (6/0.75)^{1/3}=2.0 — **far above Green's 1.519.** The sup-norm route is too lossy for
h=3 (the L^{4/3}-Fourier functional, i.e. Green's, is the working one). Confirms the brief's caution.

**C4 (Cilleruelo identity / LP ceiling).** Scoped, not completed. In groups the moment chain gives
the clean k³≤2m but the integer transfer loses interval structure (returns constant 6 < Green).
A sub-Green integer constant likely needs a new identity, not a re-optimization. Next step: the
LP-over-provable-constraints experiment (see final assessment).

---

## Candidate theorem (formal statement) and Lean target

See [`writeup/modular_theorem.md`] for the modular THEOREM (the cleanest self-contained result),
and [`writeup/upper_bound_program.md` §C3] for the certified integer-constant improvement. The
Lean target is `google-deepmind/formal-conjectures/FormalConjectures/ErdosProblems/241.lean`; its
`erdos_241.variants.upper_bound` encodes the constant (7/2)^{1/3}. A formal version of C1 would
replace (7/2)^{1/3} with our certified ≤1.51587 (needs interval-arithmetic on the kernel sum).

## Final assessment

See [`writeup/assessment.md`].
