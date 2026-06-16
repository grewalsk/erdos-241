# THEOREM CANDIDATE: a spectral upper bound for B₃ sets in finite abelian groups

Status: candidate — algebra machine-verified (code/moment_bound.py), numerically
consistent with all exact g(m) data (m ≤ 199 so far), novelty check in progress.
This file contains the complete proof; every step is elementary.

## Conventions

G is a finite abelian group, written additively, |G| = m. A ⊆ G with |A| = k.
A is a **B₃ set** if every g ∈ G has at most one representation g = a₁+a₂+a₃ with
a₁,a₂,a₃ ∈ A, counted up to reordering (multiset convention; repetition allowed).
A is a **Sidon (B₂) set** if the same holds for pair sums a₁+a₂.

Ordered representation functions: r₂(s) = #{(a,b) ∈ A² : a+b = s},
r₃(s) = #{(a,b,c) ∈ A³ : a+b+c = s}.

## Theorem 1

Let A ⊆ G be a B₃ set with k = |A| ≥ 2. Then

    m ≥ (k⁴ − 3k³ + 5k² − 4k) / (2k − 3).                                   (1)

**Corollary 1.** k³ ≤ 2m + 3k² for every B₃ set with k ≥ 2.

**Corollary 2.** Combining Corollary 1 with the counting bound k ≤ (6m)^{1/3}+O(1) to
control the k² term, k³ ≤ 2m + 3·(6m)^{2/3}, hence

    k ≤ (2m)^{1/3} (1 + O(m^{−1/3})).

**Comparison.** The counting bound (all C(k+2,3) multiset sums distinct in G) gives only
k³ + 3k² + 2k ≤ 6m, asymptotically k ≤ (6m)^{1/3} ≈ 1.817·m^{1/3}. Theorem 1 gives
k ≤ (2m)^{1/3}(1+o(1)) ≈ 1.260·m^{1/3}. Bose–Chowla provides B₃ sets in ℤ/(q³−1)ℤ of
size q = (m+1)^{1/3}, so the truth for cyclic groups lies in [m^{1/3}, (2m)^{1/3}(1+o(1))].

---

## Proof

### Step 1: every B₃ set with k ≥ 1 is Sidon.

Suppose a+b = c+d with a,b,c,d ∈ A and {a,b} ≠ {c,d} as multisets. Pick any x ∈ A
(possible, k ≥ 1). Then a+b+x = c+d+x, and the multisets {a,b,x}, {c,d,x} are distinct:
removing one copy of x from each leaves {a,b} ≠ {c,d}. (This cancellation is valid for
multisets regardless of whether x equals any of a,b,c,d.) Two distinct multiset
representations of the same group element contradict B₃. ∎

### Step 2: exact additive energies.

Define E₂ = Σ_s r₂(s)² = #{(a,b,c,d) ∈ A⁴ : a+b = c+d} and
E₃ = Σ_s r₃(s)² = #{(a₁,…,a₆) ∈ A⁶ : a₁+a₂+a₃ = a₄+a₅+a₆}.

Because A is Sidon (Step 1), each s ∈ G is the sum of at most one multiset {a,b}.
A multiset {a,b} with a ≠ b contributes r₂ = 2 (two orderings); {a,a} contributes r₂ = 1.
There are C(k,2) of the former and k of the latter, all with distinct sums, so

    E₂ = 4·C(k,2) + 1·k = 2k² − k.                                          (2)

Because A is B₃, each s is the sum of at most one multiset {a,b,c}. Multisets of three
distinct elements contribute r₃ = 6; type {a,a,b}, a ≠ b, contribute r₃ = 3; type {a,a,a}
contribute r₃ = 1. Counts: C(k,3), k(k−1), k respectively, all sums distinct
(distinctness ACROSS types is exactly the B₃ condition — e.g. a sum of type {a,a,b}
coinciding with a sum of type {c,c,c} would be two distinct multiset representations). So

    E₃ = 36·C(k,3) + 9·k(k−1) + 1·k = 6k³ − 9k² + 4k.                       (3)

Both formulas hold in any abelian group: torsion can only merge sums of distinct
multisets, which B₃ forbids. ∎

### Step 3: Fourier identities.

Let Ĝ be the character group (|Ĝ| = m), and Â(χ) = Σ_{a∈A} χ(a). For j ≥ 1, expanding
|Â(χ)|^{2j} as a 2j-fold sum and using character orthogonality
(Σ_χ χ(g) = m·1_{g=0}):

    Σ_{χ∈Ĝ} |Â(χ)|^{2j} = m · E_j,   E_j := #{(a₁..a_{2j}) : a₁+…+a_j = a_{j+1}+…+a_{2j}}.

With E₁ = k (pairs a₁ = a₂), E₂, E₃ as above, and Â(χ₀) = k at the trivial character:

    S₁ := Σ_{χ≠χ₀} |Â(χ)|²  = mk − k²
    S₂ := Σ_{χ≠χ₀} |Â(χ)|⁴  = m(2k² − k) − k⁴                              (4)
    S₃ := Σ_{χ≠χ₀} |Â(χ)|⁶  = m(6k³ − 9k² + 4k) − k⁶.

All three are sums of nonnegative reals y_χ := |Â(χ)|², y_χ², y_χ³.

### Step 4: Cauchy–Schwarz on the nontrivial spectrum.

Writing y² = y^{1/2} · y^{3/2} and applying Cauchy–Schwarz over {χ ≠ χ₀}:

    S₂² ≤ S₁ · S₃.                                                          (5)

### Step 5: algebra.

Substituting (4) into (5) and expanding (machine-verified in code/moment_bound.py;
also a routine hand computation):

    P(m,k) := S₁S₃ − S₂² = m·(α·m + β),
    α = k²(k−1)(2k−3),
    β = −k³(k⁴ − 4k³ + 8k² − 9k + 4) = −k³(k−1)(k³ − 3k² + 5k − 4).

For k ≥ 2 we have α > 0 and m > 0, so P ≥ 0 forces

    m ≥ −β/α = k(k³ − 3k² + 5k − 4)/(2k − 3) = (k⁴ − 3k³ + 5k² − 4k)/(2k − 3),

which is (1). ∎

### Corollary 1 proof.

2(k⁴ − 3k³ + 5k² − 4k) − (k³ − 3k²)(2k − 3) = 3k³ + k² − 8k = k(3k² + k − 8) > 0 for
k ≥ 2, so (1) implies the strict bound m > (k³ − 3k²)/2, a fortiori m ≥ (k³ − 3k²)/2,
i.e. k³ ≤ 2m + 3k². ∎

---

## Remarks

**R1 (tightness at k = 2).** A = {0,1} ⊆ ℤ/4ℤ has triple sums 0,1,2,3 — distinct — and
(1) demands m ≥ 4. So the inequality is achieved.

**R2 (where the strength comes from).** A flat nontrivial spectrum (y_χ ≡ S₁/(m−1) ≈ k)
would have 4th moment ≈ k²m. But Sidon-ness forces the 4th moment to be ≈ 2k²m — twice
the flat value — because E₂ cannot drop below its trivial-solution count 2k² − k. The
spectrum must therefore carry substantial mass at level y ≳ 2k. The B₃ condition caps the
6th moment at ≈ 6k³m, limiting how much such mass is affordable: by (5),
6k³m ≳ S₃ ≥ S₂²/S₁ ≈ (2k²m)²/(km) = 4k³m + k⁶/m-term, and the leftover budget 2k³m is
what bounds k⁶ ≤ 2k³m, i.e. k³ ≤ 2m. The h = 2 analogue of this chain is the classical
proof of k² ≤ m + O(√m)-type bounds for Sidon sets in groups; the point here is that for
h = 3 the *exactness of the lower-order energy* (E₂) is itself a consequence of B₃, which
is what the integer-case literature (Chen, Graham, Cilleruelo, Green) exploits in
interval-specific ways — the group setting isolates it cleanly.

**R3 (optimality within the method).** The proof uses only the first three power sums of
the nontrivial spectrum. For y ≥ 0 and any c > 0, y(y−c)² ≥ 0 gives S₃ ≥ 2cS₂ − c²S₁,
optimized at c = S₂/S₁ — which reproduces (5) exactly. The moment problem
min{E[y³] : E[y] = 1, E[y²] = 2, y ≥ 0} has value exactly 4, attained by the two-level
law y ∈ {0, 2} with equal mass. Hence **no argument using only (S₁, S₂, S₃) can improve
the constant 2**. Improvement requires new information: phase relations among Â(χ),
higher moments (E₄ is NOT controlled by B₃), or structural input about two-level spectra.

**R4 (integer corollary, for calibration).** A B₃ set A ⊆ {1,…,N} embeds into ℤ/Mℤ with
M = 3N preserving B₃: its triple sums are distinct integers in [3, 3N], an interval
containing 3N−2 integers any two of which differ by at most 3N−3 < 3N = M, so reduction
mod M is injective on them and no collision is created. (Indeed any M ≥ 3N−2 works.)
Theorem 1 then gives k³ ≤ 2M + 3k² = 6N + 3k², i.e. asymptotic constant 6 — below trivial
18 but ABOVE Chen–Graham's 4 and far above Green's 3.5. So as an *integer* bound this is
not competitive; its value is entirely in the clean modular setting, where no interval
exists and the interval-exploiting arguments of Green/Cilleruelo do not apply. The spectral
input here (exactness of E₂ from B₃⇒Sidon, exactness of E₃ from B₃, the 6th-moment
Cauchy–Schwarz) is what survives the loss of interval structure.

**R5 (empirical landscape, exact g(m) ≤ 199).** First m admitting a size-k B₃ set in
ℤ/mℤ: k=2→4, 3→13, 4→30, 5→65, 6→117. The bound (1) requires m ≥ 4, 11, 26, 51, 90
respectively: consistent, with 15–30% slack. The density ratio k³/m at those record m:
2.00, 2.08, 2.13, 1.92, 1.85 — small-k transients above 2 are permitted by the finite
form (1) (which only forces k³/m ≤ 2.46 at k=4); the asymptotic constraint is k³/m ≤ 2+o(1),
and the data trend after k=4 is decreasing. Whether lim sup k³/m = 1 (Bose–Chowla optimal
in groups) or something in (1, 2] is the modular shadow of Erdős #241.
