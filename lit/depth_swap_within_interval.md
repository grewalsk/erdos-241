# Depth local search on Bose–Chowla B_3 sets within the fixed interval [1, q^3-1]

**Question (Erdős #241 angle).** The Bose–Chowla (BC) construction gives a B_3
set of size exactly `q` inside `[1, n]`, `n = q^3 - 1`, with ratio
`|A| / n^{1/3} -> 1` from above. Single-element augmentation of BC is known to
collapse to **0 addable** elements for `q >= 13`. The family explored here:
**DEPTH local search** — remove `j` elements of the BC set and greedily re-add
`j+1` (or more) elsewhere, staying inside the SAME interval `[1, q^3-1]`, trying
to net-increase size beyond `q`. If we ever reach size `q+1` inside `[1,q^3-1]`,
that beats ratio 1.

All sets reported were verified with the standalone `verify.is_b3` (imported, not
reimplemented). The depth search uses an exact set-based `addable(S,N)`
(antitone, equivalent to the FFT mask in `bose_chowla.py`).

## What the search does

For each `q`:
1. **Greedy baseline.** Greedily augment the full BC set inside `[1,n]` (both
   ascending and descending candidate order). This is the j=0 case.
2. **Depth-j swap.** For every `j`-subset `R` of the BC set, remove `R`, then
   greedily re-fill the remainder inside `[1,n]`. Record the largest B_3 set
   obtained. Tried up to `j = 4` at the critical `q = 13`, and up to `j = 3` at
   `q = 16, 17, 19`.

## Results

| q  | n=q^3-1 | BC size | best size found | gain over q | ratio = |A|/n^{1/3} |
|----|---------|---------|-----------------|-------------|---------------------|
| 4  | 63      | 4       | 6               | **+2**      | 1.508               |
| 5  | 124     | 5       | 7               | **+2**      | 1.404               |
| 7  | 342     | 7       | 9               | **+2**      | 1.287               |
| 8  | 511     | 8       | 10              | **+2**      | 1.251               |
| 9  | 728     | 9       | 11              | **+2**      | 1.223               |
| 11 | 1330    | 11      | 12              | **+1**      | 1.091               |
| 13 | 2196    | 13      | 14              | **+1**      | 1.077               |
| 16 | 4095    | 16      | 16              | **0**       | 1.000               |
| 17 | 4912    | 17      | 17              | **0**       | 1.000               |
| 19 | 6858    | 19      | 19              | **0**       | 1.000               |
| 23 | 12166   | 23      | 23              | **0**       | 1.000               |
| 25 | 15624   | 25      | 25              | **0**       | 1.000               |
| 27 | 19682   | 27      | 27              | **0**       | 1.000               |

(`ratio` uses `n = q^3 - 1`. Using `N = max(A)` instead gives slightly larger
numbers for small q, e.g. 1.092 at q=13, because the achieved max element sits
below `n`; this does not change any conclusion. The 1.034/1.046 figures at
q=25/27 for `max(A)`-normalization are an artifact of BC's largest element
landing well below `q^3-1` — the set is still size exactly `q`, gain 0.)

### The headline positive at q=13 (where single augmentation gives 0)

At `q = 13` single-element augmentation of BC is locally maximal (0 addable —
matches the known collapse). The depth-1 swap (remove one BC element, then
greedily re-fill) breaks past it and reaches **size 14 > q = 13**:

```
A = [1, 99, 105, 129, 194, 351, 683, 844, 1189, 1413, 2073, 2081, 2100, 2107]
```

`verify.is_b3(A) == True`, size 14, max = 2107 < 2196. So inside `[1, 2196]`
one can indeed exceed `q`. Going deeper (j up to 4) does **not** get past 14.

## Honest verdict: finite-size effect, decays to 1

The depth swap genuinely does what the question asks **for small q**: it exceeds
size `q` inside `[1, q^3-1]`, even at `q = 13` where the simpler single-element
move is dead. But the effect is a **finite-size artifact, not a sustained
improvement**:

- The net gain over `q` **shrinks** as `q` grows: **+2** (q=4–9) -> **+1**
  (q=11,13) -> **0** (q >= 16). Adding more depth does not rescue it: at
  q = 16, 17, 19 even removing 3 elements and re-filling cannot beat `q`.
- Consequently the ratio **decays monotonically to 1** (1.51, 1.40, 1.29, 1.25,
  1.22, 1.09, 1.08, then 1.000, 1.000, ...), i.e. it converges to the same
  `1 + o(1)` as plain BC — there is **no fixed `c > 0`** that it sustains.

The mechanism is intuitive: at fixed `q`, the BC set occupies a vanishing
density `q / q^3 = q^{-2}` of `[1, q^3-1]`, so for small `q` there is slack —
the interval is "loose" enough that a B_3 set of size a couple above `q` fits,
and local search finds it. As `q` grows the BC set becomes essentially optimal
for the interval (it nearly saturates the cube-root bound), the slack closes,
and remove-and-refill can no longer net any element. The few extra elements
findable at tiny `q` are a constant-order surplus that is swamped by `q^{1/3}`
in the denominator.

**Conclusion.** This depth-local-search family does **not** sustainably exceed
ratio `1 + c` for any fixed `c > 0`. It only beats 1 at small `q` (up to q=13)
and decays to exactly 1 by q=16. It is a finite-size effect, consistent with the
Bose–Chowla baseline being asymptotically `1 + o(1)`.

## Reproduction

Scripts in `/tmp` (`depth_search.py`, `depth2.py`, `depth3.py`, `depth_all.py`,
`depth_big.py`, `depth2_big.py`), run with
`PYTHONPATH=.../erdos-241/code python3 ...`. Every candidate set is checked with
`verify.is_b3`; BC sets generated by `bose_chowla.bose_chowla`.
