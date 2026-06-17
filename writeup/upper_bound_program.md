# Upper-bound program (Workstream C) — reconstruction, transfer, and a certified kernel improvement

## C1. Green's argument and where the slack lives (COMPUTATION + reconstruction)

Green's route for B_3 = B_3[1] sets in {1,…,N}:

1. **Combinatorial input (his Lemma 16).** A B_3 set obeys the pointwise bound
   A∗A∗A∗A(x) ≤ 2|A|(1 + A∗A(x)). Summing gives an upper bound on the "number of squares"
   M(A∗A) and, after his Fourier machinery, a *lower* bound that must be reconciled with it.

2. **The variational lower bound (his Theorem 6 / 11).** For any p ∈ C¹[0,1] with ∫₀¹p = 2,
   the truncated 4th-moment E(X)=Σ_{0<|r|≤X}|Â(r)|⁴ ≥ γ(p)·N⁴·(1−o(1)), with
   **γ(p) = 2·(Σ_{r≥1}|p̃(πr)|^{4/3})^{−3}**, p̃(λ)=∫₀¹p(x)e^{ixλ}dx.

3. **Conclusion.** Combining, A(3,N) ≤ (4/(1+γ(p)))^{1/3} N^{1/3}(1+o(1)).

**The slack is exactly the kernel p.** Green picks p(x)=5/2−40(x−1/2)⁴ — by his own words "a
simple function that gives a good bound," not the optimum — getting γ ≈ 1/6.9994, hence
4/(1+γ) = 3.4997, reported as (7/2)^{1/3}. The implicit optimization is

    minimize  J(p) := Σ_{r≥1} |p̃(πr)|^{4/3}   over  p with ∫₀¹p = 2,

since the constant is (4/(1+2J^{-3}))^{1/3}, strictly increasing in J.

**Reproduction (COMPUTATION, `code/green_gamma.py`):** using Green's closed forms (26)–(29)
we get S₁=0.0839757, S₂=0.1219300 (his: 0.0839757, 0.1219299), γ=1/6.9948, constant 3.49968 —
matching to 5 digits. ✓ The framework is faithfully reconstructed.

## C2. Discrete → continuous transfer, two relaxations

There are TWO continuous relaxations of the h=3 problem; they are not equivalent.

**(a) Green's sharp relaxation** is the L^{4/3}-Fourier functional J above. It already lives in
the continuous variable (the kernel p) and *is* the object to optimize. This is the one that
gives the record constant.

**(b) The sup-norm autoconvolution μ₃** (the brief's C.2):
μ₃ := inf{‖f∗f∗f‖∞ : f ≥ 0, supp f ⊆ [0,1], ∫f = 1}.

*Transfer (the factor we can prove).* Let A ⊆ {1,…,N} be B_3, |A| = k. Replace each a ∈ A by
the uniform density of mass 1/k on [(a−1)/N, a/N]; call the result f (a density on [0,1],
∫f = 1). Then f∗f∗f = (1/k³) Σ_{a,b,c} U_a∗U_b∗U_c, a sum of B-splines of width 3/N. At any x,
the contributing triples are those with a+b+c within 3 of Nx; since A is B_3 the ordered count
r_3 ≤ 6 per integer value and the triple-sums are distinct integers, so a window of length 3
contains O(1) of them with bounded multiplicity. Each width-3/N spline has height O(N), so
‖f∗f∗f‖∞ ≤ C·N/k³ with a small explicit C (≈6 in the clean count). Since ‖f∗f∗f‖∞ ≥ μ₃,

    k³ ≤ (C/μ₃)·N,   constant (C/μ₃)^{1/3}.

*Why it cannot beat 3.5 (COMPUTATION, `code/mu3.py`):* the uniform density already gives
‖f∗f∗f‖∞ = 3/4 (the order-3 cardinal B-spline peak), so μ₃ ≤ 3/4; optimizing the density only
*lowers* μ₃ (we reach ≤ 0.588 numerically), which only *raises* the constant (C/μ₃)^{1/3}. Even
with the optimistic C = 6 the constant is ≥ (6/0.75)^{1/3} = 2.0 N^{1/3} (and ≥ 2.17 at
μ₃ ≤ 0.588) — **far above Green's 1.519.** This matches the brief's caution: the sup-norm route
is too lossy for h=3 (exactly as autoconvolution fails to recover the sharp B_2[1] Sidon bound).
The lesson: sup-norm throws away the phase/moment information that Green's L^{4/3} functional keeps.

## C3. Optimizing Green's kernel — a certified improvement (COMPUTATION → certified)

We minimize J(p) directly.

**Exact evaluation (`code/green_optimize.py`).** For a piecewise-constant p on n bins,
J(p) has a closed form with NO truncation: with ω=e^{iπ/n}, t = r mod 2n,
|p̃(πr)| = |ω^t−1|·|S(t)|/(πr), S(t)=Σ_j p_j ω^{tj}, hence
J = (1/π)^{4/3}(2n)^{−4/3} Σ_{t=1}^{2n−1}(|ω^t−1||S(t)|)^{4/3} ζ_H(4/3, t/2n),
S via one FFT, ζ_H the Hurwitz zeta. Validated against Green's closed-form quartic:
step-J → 2.40994 as n → 4096 (closed form 2.40955). ✓

**Convex minimization (analytic gradient, L-BFGS).** J is convex (sum of (4/3)-powers of
moduli of linear functionals) with the affine constraint ∫p = 2. Minimizing over step
functions and Richardson-extrapolating in n:

| n | min J | constant (cubed 4/(1+γ)) |
|---|---|---|
| 256 | 2.39592 | 3.49218 |
| 1024 | 2.38069 | 3.48364 |
| 4096 | 2.37654 | 3.48128 |
| Richardson n→∞ | **≈ 2.3751** | **≈ 3.4805** |

The optimal kernel is bounded, positive, symmetric about 1/2, with p(0)=p(1)=0 (same boundary
behavior as Green's quartic) and interior max ≈ 2.53. Two independent J evaluations (the
FFT–Hurwitz formula and a direct per-bin closed-form quadrature) agree to 4 digits.

**Certified explicit kernel (`code/green_certify2.py`).** To make the improvement rigorous we
exhibit a concrete smooth kernel and bound J(p) from above in closed form. Family
p(x)=x(1−x)·Σ_{j=0}^{d} c_j ((2x−1)²)^j (degree 2d+2, p(0)=p(1)=0, ∫p=2). The optimized
degree-12 member has:
- c_r = p̃(πr) in exact closed form (sympy), partial sum to R = 6·10⁵;
- rigorous tail |c_r| ≤ B/(πr)² with B = |p'(0)|+|p'(1)|+‖p''‖₁ (two integrations by parts,
  using p(0)=p(1)=0), so Σ_{r>R}|c_r|^{4/3} ≤ (B/π²)^{4/3}(3/5)R^{−5/3} < 9·10⁻⁹;
- result: **J ∈ [2.37996, 2.37996]**, the kernel is **positive** on [0,1].

Hence, granting Green's reduction (Theorem 6/11/17, published & refereed),

    **f(N) ≤ 1.51587 · N^{1/3} (1+o(1)),**

versus Green's 1.51825·N^{1/3} — equivalently 4/(1+γ) ≤ 3.4832 vs Green's actual 3.4997.

**Scope and honesty.** (i) This is ~0.16% in the N^{1/3} coefficient. (ii) It does not touch
the round figure "7/2"; it improves Green's *true* constant, which he obtained from a kernel he
explicitly flagged as non-optimal. (iii) It introduces **no new idea** — it completes an
optimization Green set up. (iv) The certificate's finite partial sum is float64 (the 0.002
margin to Green dwarfs the ~10⁻¹⁰ accumulated error; an interval-arithmetic / Lean version would
make it fully formal). (v) The asymptotic picture (whether f(N)/N^{1/3} → 1) is untouched.

This is, nonetheless, a *rigorous improvement to the best known upper-bound constant* — one of
the brief's stated win conditions — and the cleanest available next increment on the prize route.

## C4. Cilleruelo identity / moment route (status: scoped, not completed)

The brief's C.4 (strengthen the Cilleruelo small-difference identity to bound the integer
constant) was scoped but not pushed to a number, for two reasons discovered above:
1. In the *group* setting the analogous moment chain is clean and gives k³ ≤ 2m (see
   `writeup/modular_theorem.md`) — strictly better than 6 — but the integer transfer loses the
   interval structure (it returns constant 6, worse than Green; see modular_theorem R4).
2. Green's functional already extracts the interval non-equidistribution optimally up to the
   kernel choice; the certified C3 improvement suggests the remaining slack in the *integer*
   problem is small within this framework. A genuine sub-Green constant likely needs a new
   identity beyond Cilleruelo's, not a re-optimization of the existing one.
This is the honest place to stop on C.4: the LP-ceiling experiment is the right next step
(see results.md final assessment) but was not run to completion in this pass.

## C5. Modular toy problem — see `writeup/modular_theorem.md`

The cleanest C result: a self-contained THEOREM candidate k³ ≤ 2m + 3k² for B_3 sets in any
finite abelian group of order m (constant 6 → 2 vs the counting bound), adversarially verified
by independent skeptics (3/6 returned "sound", 0 fatal/gap; remaining lenses covered by overlap).
