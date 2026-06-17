# Union of dilated/translated Bose-Chowla copies — B_3 ratio study

**Question.** Take a Bose-Chowla (BC) B_3 set `A0` of size `q` in `[1, q^3-1]`.
Form `A = union_i (d_i * A0 + t_i)` for a small number of copies. Can a union of
2-3 scaled BC blocks sustain a ratio `|A| / N^{1/3}` (with `N = max A`) above
`1 + c` for fixed `c > 0` as `N -> infinity`?

**Verification.** Every candidate set was checked with the standalone brute-force
verifier `code/verify.py` (`is_b3`, multiset convention) — imported, never
reimplemented. BC sets came from `code/bose_chowla.py`. Scoring used the honest
`N = max(A)` (not the modulus `q^3-1`).

## Answer: NO. The family does not sustain ratio > 1. The only wins are tiny-N artifacts.

### 1. Pure translates never work
`A0 union (A0 + T)` is **never** B_3 for `|A0| >= 3`, for any offset `T`. A full
scan of `T` over `[1, 6*max(A0)]` (step 1) for `q in {4,5,7,8,9,11}` found zero
working offsets. Three independent offsets (3-copy) also found nothing.

The obstruction is explicit. For `A0` from `q=5` and `T=500`, the verifier's first
collision is
```
(1, 1, 513) and (1, 13, 501)   both sum to 515
```
i.e. `a_i + a_j + (a_k + T)` collides with `a_i + a_k + (a_j + T)` because both
equal `a_i + a_j + a_k + T`. Mixing "two-from-block-0 + one-from-block-1" against a
reshuffled split is unavoidable once the block has more than one pairwise sum.
Translation alone cannot separate these — every translated copy sum-interacts with
the original at the same additive offset.

### 2. Dilated copies can work at tiny q, but the ratio collapses
`A = A0 union (D*A0 + T)`. Dilation makes the second copy multiplicatively
independent, which CAN kill the cross-block collisions — but only if `D` is large
enough, and the required `D` grows with `q`. Exhaustive `(D, T)` search
(`D <= 30`, `T` full range):

| q | \|A0\| | max(A0) | min viable D | best union ratio | single BC ratio |
|---|--------|---------|--------------|------------------|-----------------|
| 4 | 4  | 60  | 4    | **1.1696** | 1.0217 |
| 5 | 5  | 104 | 9    | 0.9461     | 1.0632 |
| 7 | 7  | 297 | 14   | 0.8277     | 1.0492 |
| 8 | 8  | 422 | none (<=30) | — | 1.0666 |
| 9 | 9  | 723 | none (<=22) | — | 1.0028 |

The best union ratio **decays monotonically**: 1.17 -> 0.95 -> 0.83 -> no valid set.
It crosses below 1 already at `q=5` and below the *single-BC* baseline immediately
after `q=4`. For `q >= 8` no dilation `<= 30` yields a B_3 union at all.

### 3. Why it must decay (the mechanism)
To clear all cross-block triple collisions, the dilation `D` must grow with `q`
(denser elements -> more collision opportunities -> larger separation required;
empirically min-D = 4, 9, 14, ... ). But then `N = max(A) ~ D * max(A0) ~ D * q^3`,
while `|A|` only doubles to `2q`. So
```
ratio = 2q / (D q^3)^{1/3} = 2 / D^{1/3} * (q / q) ... -> 2 / D^{1/3}
```
and since `D -> infinity` with `q`, the ratio `-> 0`. A constant number of dilated
copies buys a constant factor in `|A|` but pays a `D^{1/3}` -> infinity penalty in
the span. There is no fixed `c > 0`.

### 4. The one "win" is a finite-size effect
The single genuine improvement is at `q = 4`:
```
A0 = {1, 6, 29, 60},  D = 4, T = 80
A  = {1, 6, 29, 60, 84, 104, 196, 320}   (verified B_3 by verify.is_b3)
|A| = 8,  N = 320,  ratio = 1.1696
```
This is real (independently re-verified) but isolated: it is the smallest case,
and the construction's advantage evaporates by `q = 5`. It is a textbook small-N
artifact, not a sustained family. (Compare: greedy augmentation of a single BC set
gives ratio ~1.57 at q=2 decaying to ~1.09 at q=11 — same story, different family.)

## Conclusion
Unions of 2-3 dilated/translated BC blocks **do not** beat the Bose-Chowla
`1 + o(1)` ratio in any sustained way. Pure translates are never B_3; dilated
unions require a dilation that grows with `q`, forcing the span up faster than the
size, so the ratio decays to 0. The only set that exceeds the baseline (q=4,
ratio 1.17) is a finite-size accident. No fixed `c > 0` is achievable by this
family.

Reproduce: `/tmp/union_bc.py`, `/tmp/union_search.py`, `/tmp/obstruction.py`,
`/tmp/dilate_search.py`, `/tmp/dilate_fine.py`, `/tmp/scaling.py`,
`/tmp/final_verify.py`.
