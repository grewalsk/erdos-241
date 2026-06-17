# Literature verification memo — Erdős Problem #241

Compiled 2026-06-12. Sources fetched into `lit/`; every load-bearing claim below is
tagged by the strongest source actually read (PRIMARY = original paper read directly;
AUTHORITATIVE = erdosproblems.com / OEIS / Lean repo; SURVEY = O'Bryant DS11).

## Problem statement and status (AUTHORITATIVE)

From **erdosproblems.com/241** (verbatim, `lit/raw_forum_241.html`):
> "Originally asked to Erdős by Bose. Bose and Chowla [BoCh62] provided a construction
> proving one half of this, namely (1+o(1))N^{1/3} ≤ f(N). The best upper bound known to
> date is due to Green [Gr01], f(N) ≤ ((7/2)^{1/3}+o(1))N^{1/3} (note that (7/2)^{1/3}≈1.519).
> ... Is it true that f(N) ∼ N^{1/3}?"

Status: **OPEN.** Tags: additive combinatorics, Sidon sets. Guy's *Unsolved Problems in
Number Theory* problem C11.

## Convention (AUTHORITATIVE — fixes our whole project)

**Lean** `google-deepmind/formal-conjectures/FormalConjectures/ErdosProblems/241.lean`
defines `f N r` via
```
∀ m₁ m₂ : Multiset ℕ, m₁.card = r → m₂.card = r →
  (∀ x ∈ m₁, x ∈ A) → (∀ x ∈ m₂, x ∈ A) → m₁.sum = m₂.sum → m₁ = m₂
```
over `A ⊆ Icc 1 N`, target `(fun N ↦ f N 3) ~[atTop] (fun N ↦ N^(1/3))`.
This is exactly our **multiset / repetition-allowed** convention on the interval
{1,…,N}. The lower- and upper-bound variants are stated as separate `solved` theorems with
constants 1 and (7/2)^{1/3}. **Any theorem candidate we phrase should target `f N 3` and the
multiset predicate above.**

**OEIS A387704** (AUTHORITATIVE, `lit/a387704.txt`): "Size of the maximal subset S of
{1,2,...,n} such that for all a, b, c in S not necessarily distinct, a+b+c is unique up to
permutation." Offset 0,3; example: n=5 → {1,2,5}, "adding 3 would create the duplicate sums
1+1+3=1+2+2." This is **identical** to our convention (repetition allowed; multisets).
Author Sharvil Kesarwani, Dec 2025; b-file to n=150. Cross-refs A143824 (2-sums = Sidon),
A385931 (B_3 with **distinct** 3-sums — the other convention).
→ Our exact solver reproduces all 151 terms n=0..150 with **zero** mismatches and extends
the b-file to n=159 (see results.md, COMPUTATION A1).

## Lower bound (PRIMARY construction verified)

**Bose–Chowla 1962** ("Theorems in the additive theory of numbers", Comment. Math. Helv.):
for prime power q, a B_3 set of size q inside [1, q^3−1], giving f(N) ≥ (1+o(1))N^{1/3}.
Construction (from O'Bryant DS11 §3.3, `lit/obryant-ds11.txt` lines 206–234):
> "Bose_h(q,θ,k) := {a ∈ [q^h−1] : θ^a − kθ ∈ F_q}. Bose_h(q,θ,1) is a B_h (mod q^h−1) set."
For h=3, k=1: a with θ^a = θ + c, c ∈ F_q. **Our `bose_chowla.py` implements exactly this and
the standalone verifier confirms B_3 (mod q^3−1 and over ℤ) for all q ∈ {2,…,64}.**
The constant **1 has not been improved** (no source found doing so; consistent with
erdosproblems.com calling only the upper side "the best bound known").

## Upper bounds (chronology; record PRIMARY-verified)

| bound on f(N) | constant ×N^{1/3} | source | verification |
|---|---|---|---|
| trivial C(k+2,3) ≤ 3N | 18^{1/3} ≈ 2.621 | folklore | — |
| (4N)^{1/3} | 1.5874 | Chen; Graham (indep.) | SURVEY |
| (3.996N)^{1/3} | ≈1.587 | Graham via Li | SURVEY |
| (4/(1+16/(π+2)^4))^{1/3} | ≈1.5756 | Cilleruelo 2001 | SURVEY |
| **(7/2)^{1/3}** | **≈1.5183** | **Green 2001** | **PRIMARY** |

**Green 2001** ("The number of squares and B_h[g] sets", Acta Arith. 100, 365–390),
read directly (`lit/green-bhg.txt`):
- **Theorem 17 (verbatim):** "A(3,N) ≤ (7/2)^{1/3} N^{1/3}(1+o(1))." This is the record.
  Confirmed independently by O'Bryant DS11: "σ₃ ≤ (7/2)^{1/3} < 1.519" (two locations).
- **Proof skeleton (h=3):** Lemma 16: a B_3 set obeys A∗A∗A∗A(x) ≤ 2|A|(1+(A∗A)(x)); fed
  into the "number of squares" machinery M(f)=Σ_x(f∗f)(x)² = (1/(2N+v))Σ_r|f̂(r)|⁴.
  Theorem 6 gives, for any p ∈ C¹[0,1] with ∫₀¹p = 2,
  E(X)=Σ_{0<|r|≤X}|f̂(r)|⁴ ≥ γ(p)N⁴(1−o(1)), where **γ(p) = 2(Σ_{r≥1}|p̃(πr)|^{4/3})^{−3}**,
  p̃(λ)=∫₀¹p(x)e^{ixλ}dx. The B_3 bound is A(3,N) ≤ (4/(1+γ(p)))^{1/3}N^{1/3}(1+o(1)).
- **Where the slack is (the key finding):** the constant comes ENTIRELY from the choice of
  kernel p. Green writes (p.13, verbatim): "we have not been able to give a best possible
  choice in closed form. **A simple function that gives a good bound is** p(x)=5/2−40(x−1/2)⁴,"
  which yields γ(p) ≈ 1/6.9994, i.e. 4/(1+γ) = 3.4997, rounded up to 7/2 = 3.5.
  → He explicitly did **not** claim kernel-optimality. **We reproduce his γ to 5 digits and
  then minimize J(p)=Σ|p̃(πr)|^{4/3} over kernels** (results.md, COMPUTATION C1): inf J ≈ 2.3751
  (Green's quartic: 2.4095), giving a **certified** constant ≤ 1.51587 N^{1/3} < 1.51825 — a
  small but rigorous improvement of Green's *actual* value (the headline "3.5" is unaffected).

## Post-2001 (no g=1, h=3 improvement found)

- **Timmons 2016** ("An improved upper bound for B_h[g] sets", or the B_3[g] paper): improves
  B_3[**g**] only for **g ≥ 2**; explicitly NOT the g=1 problem here. (Brief's caution confirmed;
  no source contradicts it.)
- erdosproblems.com still lists **Green** as "the best upper bound known to date" (2025/26
  page), i.e. no improvement to (7/2)^{1/3} for g=1 has been recorded. The Lean repo's
  `upper_bound` variant also encodes (7/2)^{1/3} as current.
- **Triple autoconvolution** inf‖f∗f∗f‖∞ (the h=3 analogue of Martin–O'Bryant autoconvolution
  work): we found no published value; we compute it (results.md, COMPUTATION C3) and show this
  *sup-norm* relaxation is too lossy (constant ≥ 2.0), distinct from Green's working *L^{4/3}*
  functional.

## Modular case — novelty of our bound (PARTIAL — see open items)

The brief's in-house bound for B_3 sets in ℤ/mℤ — k³ ≤ 2m + 3k² via E_2 = 2k²−k (Sidon),
E_3 = 6k³−9k²+4k (B_3), Parseval, and a 6th-moment Cauchy–Schwarz — has the **h=2 shadow**
being the classical Sidon bound k² ≤ m+O(√m). The h=2 version is folklore (e.g. it underlies
the (7/2)/B_h machinery and the standard Singer/planar-difference-set extremal counts). We did
**not** locate the exact h=3 statement k³ ≤ (2+o(1))m in the literature in this pass; the most
likely homes (Cilleruelo–Ruzsa–Trujillo "Upper and lower bounds for finite B_h[g] sequences";
B_h-sets-in-finite-abelian-groups literature; O'Bryant DS11 modular section) should be checked
before any novelty claim. Our position: **treat the modular bound as a clean, self-contained
THEOREM candidate, explicitly flagged as "possibly folklore," not as a novelty claim.**

### Open items (literature, honest)
1. Exact provenance of the modular k³ ≤ (2+o(1))m bound (likely folklore; confirm or attribute).
2. Cilleruelo (2001) "New upper bounds for finite B_h sequences" — read the h=3 derivation
   directly (we have the constant from the brief/survey, not yet the primary proof).
3. Whether anyone has numerically optimized Green's γ(p) before (our improvement assumes not;
   Green's own "no closed form / simple function" remark suggests the optimization was left open).
