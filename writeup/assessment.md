# Final honest assessment — Erdős #241

## What was actually achieved

- **Exact landscape (A).** f(N) computed to N=159, reproducing OEIS A387704 exactly and
  extending it 9 terms; g(m) (modular) to m=200. Both via solvers cross-checked against a
  verifier-driven ground truth, with every witness independently re-verified.
- **A clean modular THEOREM candidate (C).** k³ ≤ 2m+3k² for B₃ sets in any finite abelian
  group — counting constant 6 driven down to 2 — with a fully elementary, adversarially
  vetted proof. Almost certainly the right modular constant order; possibly folklore.
- **A certified improvement of Green's constant (C).** f(N) ≤ 1.51587·N^{1/3} vs Green's
  1.51825, by certifying an optimized kernel in Green's own variational bound. Rigorous but
  small and idea-free — it finishes an optimization Green explicitly left open.
- **Lower-bound evidence (B).** Bose–Chowla sets are locally maximal (0 addable elements for
  all q≥13); no structured family beat ratio 1+c. Evidence *for* the conjecture.
- **A documented dead end (C).** The sup-norm triple-autoconvolution relaxation μ₃ caps the
  constant at ≥2.0 — too lossy; the L^{4/3}-Fourier (Green) functional is the one that works.

## Which direction is most promising

**Ranked by expected value:**

1. **Most promising for a real result: the modular problem (constant 1 vs 2 in groups).**
   We have k³ ≤ (2+o(1))m from above and m^{1/3} from Bose–Chowla below, and the exact g(m)
   data shows small extremal sets hugging the *upper* (2m)^{1/3} curve. The open modular
   question — is limsup g(m)³/m equal to 1 or 2? — is sharper, cleaner (perfect Fourier
   orthogonality), and independently publishable. **Concrete next experiment:** push the exact
   g(m) solver to m ≈ 2000 (feasible: sets have size ~(2m)^{1/3}≈16, and a dedicated
   meet-in-the-middle / clique solver should reach it), and fit g(m)³/m. If it trends to 1,
   that's strong evidence the modular truth is the BC constant and the spectral 2 is loose,
   pointing at exactly which Cauchy–Schwarz slack to close; if it plateaus above 1, that is a
   genuine structural surprise. Pair with: improve the spectral bound below 2 using the
   *two-level* structure of the extremal spectrum (Remark R3 shows the moment problem is tight
   only for a two-atom law y∈{0,2k}; real B₃ spectra are not two-level, so a 4th-moment or
   phase input should beat 2).

2. **Lower-effort, modest payoff: further certified kernel optimization (integer constant).**
   Our certified 1.51587 can be tightened toward the extrapolated inf J ≈ 2.3751 (constant
   ≈1.51544) with a higher-degree certified kernel, and formalized (interval arithmetic) to
   replace (7/2)^{1/3} in the Lean file. Real but incremental; no new mathematics.

3. **Highest ceiling, lowest probability: a new integer identity beyond Cilleruelo.** Beating
   Green by an Ω(1) margin needs information the kernel route discards. The group result shows
   the moment machinery is powerful when interval structure is removed; the missing idea is an
   *interval-aware* analogue of the B₃⇒(exact E₂,E₃) rigidity. **Concrete next experiment:**
   the C.4 LP — maximize the implied constant over all linear inequalities we can *prove*
   relating d_A(h), window-count variance (Cilleruelo's identity), and the triple-sum
   multiplicity cap — to learn the best constant this family of arguments can possibly yield
   *before* investing in a proof. If the LP optimum is ≥3.5, the route is exhausted; if it's
   below, it tells us which inequality to sharpen.

## On the main conjecture

Nothing here moves f(N)/N^{1/3} → 1. The evidence assembled is *consistent* with the
conjecture (declining exact ratios, locally-maximal BC sets, no dense outliers) but the
convergence — if real — is extraordinarily slow (ratio still ≈1.58 at N≈130), which is exactly
why the limit's existence remains open. The realistic near-term wins are the modular constant
and the certified integer-constant increment; the conjecture itself is untouched and the
"missing reduction" (an h=3 analogue of the Lindström B₄→Sidon transfer) remains absent.

## Reproducibility

All code in `code/` (verifier first); data in `data/`; sources in `lit/`; writeups in
`writeup/`. Run `python3 code/verify.py` (self-test), `python3 code/moment_bound.py` (modular
theorem checks), `python3 code/green_certify2.py` (certified kernel). Solvers: `cc -O2 fsolver.c`,
`gsolver.c`. Repo: github.com/grewalsk/erdos-241.
