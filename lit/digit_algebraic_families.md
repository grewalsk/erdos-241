# Digit-constrained & small-algebraic B_3 families (Erdős #241)

**Verdict: every family in this group either (a) is not B_3 at all beyond
trivial sizes, or (b) is B_3 but has ratio `|A|/N^{1/3}` decaying to 0. None
sustains ratio `1+c` for fixed `c>0`. The only ">1" numbers seen are tiny-N
finite-size artifacts, and the highest ones come from sets that FAIL the B_3
check (so they do not count).**

All sets verified with the standalone `verify.is_b3` (imported, not
reimplemented); multiset convention. Scripts in `/tmp/fam_*.py`.
`N = max(A)`, `ratio = |A| / N^{1/3}`. Bose–Chowla baseline → 1+o(1).

## 1. Cubes `a_i = i^3`

The set `{1^3, …, n^3}` has `ratio ≡ 1.000` *by construction* (`N=n^3`,
`|A|=n`) — but it is **not B_3** for `n ≥ 6` (three-cube near-Fermat
collisions, e.g. relations among `1,6,8,9`). Verified:

| set        | \|A\| | N      | ratio | B_3   |
|------------|-------|--------|-------|-------|
| `1..5`³    | 5     | 125    | 1.000 | True  |
| `1..8`³    | 8     | 512    | 1.000 | **False** |
| `1..50`³   | 50    | 125000 | 1.000 | **False** |

Greedy over cubes (keep only cubes that preserve B_3) IS B_3 but collapses,
because admissible cubes are sparse and `N` grows like `i^3`:

| candidates | \|A\| | N           | ratio |
|------------|-------|-------------|-------|
| `i<50`     | 22    | 110,592     | 0.458 |
| `i<100`    | 33    | 970,299     | 0.333 |
| `i<200`    | 52    | 7,414,875   | 0.267 |
| `i<800`    | 117   | 469,097,433 | 0.151 |

Ratio → 0. Cubes are the wrong density (they thin out as `i^3`).

## 2. Squares `a_i = i^2`

Squares are not even B_2 in a useful density here. Greedy-over-squares B_3
sets decay too: `i<50 → 0.921`, `i<100 → 0.782`, `i<400 → 0.554`. Monotone to 0.

## 3. Perfect / optimal rulers (Golomb = Sidon = B_2)

Optimal Golomb rulers are **B_2 by definition but not B_3**. *Every* tested
order (4–13) fails the B_3 check. Their high "ratios" (2.0–2.7) are
meaningless because `B_3 = False`. Concretely `B_2 ⊋ B_3`: a perfect ruler
controls pair sums/differences, not triple sums.

| order | marks N | ratio | B_2  | B_3   |
|-------|---------|-------|------|-------|
| 4     | 7       | 2.09  | True | **False** |
| 10    | 56      | 2.61  | True | **False** |
| 13    | 107     | 2.74  | True | **False** |

Rulers / Sidon constructions do not transfer to B_3.

## 4. Digit / base-b product (`generalized Sidon-from-digits`)

Idea: take a B_3 alphabet `S ⊆ {0,…,m-1}`, base `b`, and form all numbers
whose base-`b` digits lie in `S`. **This product set is NOT B_3** (digit
permutations collide under 3-fold sums). For `S={1,2,5}`, `b=16`:
`ndig=1` B_3 True; `ndig=2,3` B_3 **False**. The ">2" ratios are spurious.

## 5. Direct sum / packing `A + s·B` of two B_3 sets

A natural "B_3-by-construction" hope: pick spacing `s` large enough that the
low and high blocks decouple, expecting `A + s·B` to be B_3 when `A,B` are.
**This is false at every spacing** — including `s = 10^6`, far past any carry
threshold. Reason (the load-bearing finding): the triple sum is
`(Σ a_i) + s·(Σ b_i)`; for large `s` the `a`-sum and `b`-sum decouple, but
the **pairing** between which `a` came with which `b` is lost. Explicit
collision for `{1,2,5}+s·{1,2,5}` at *any* `s`:

```
(1+1s, 1+1s, 2+2s)  and  (1+1s, 2+1s, 1+2s)   both sum to 4 + 4s
   pairs {(1,1),(1,1),(2,2)}      pairs {(1,1),(2,1),(1,2)}
```

Both give `(Σa, Σb) = (4,4)` as different multisets of pairs. So **no
spacing** rescues it. Searched `s ∈ [1,120)` for `{1,2,5}` and `{1,2,5,14}`:
zero B_3-preserving spacings. Consequence: B_3 sets **do not compose** the way
B_h sets are usually concatenated; digit/direct-sum packing cannot build large
B_3 sets at all, let alone dense ones. (This is exactly why the genuine
construction, Bose–Chowla, is *algebraic* — it encodes the pairing into a
single field element `θ+c`, which a digit decomposition cannot.)

## 6. Why the only ">1" survivor is a finite-size blip

The one genuinely-B_3 family in this neighborhood that breaks 1 is the greedy
B_3 sequence (A051912, covered in `greedy_b3_family.md`): ratio peaks at
**1.754 at the trivial `{1,2,5}` (k=3, N=5)**, then decays monotonically
through 1.0 (≈ k=15, N≈3100) to 0.52 by k=64. Verified prefix
`{1,2,5,14,33,72,125,219,376,573}` (k=10) has ratio 1.204, on the decay slope.
The peak is driven entirely by the tiny leading terms; it is not sustained.

## Bottom line

Digit-restricted, ruler-based, and small-algebraic (cube/square/power)
families give **no** sustained excess over ratio 1. The structural reason is
sharp: B_3 is a *triple-sum* condition that does not factor through digit
blocks or through Sidon/ruler (pair) structure. The only constructions that
hold ratio ≈ 1 at large N are the algebraic Bose–Chowla family — and the
question of beating `1+c` is not resolved upward by any family here.
