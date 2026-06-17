# Greedy B_3 family (Erdős #241): the thin baseline

**Verdict: greedy FALLS BELOW Bose-Chowla. It beats ratio 1 only at tiny N
(finite-size effect) and then decays to 0 as a slow power law. It cannot
sustain ratio 1+c for any fixed c>0.**

## Construction

Pure greedy: `S = {1}`; repeatedly add the smallest integer keeping `S` a
B_3 set, up to a ceiling `N`. Score `ratio = |S| / Nmax^{1/3}`, `Nmax = max(S)`.
Implemented in `/tmp/greedy_b3.py`; admissibility maintained incrementally via
pair-sum and triple-sum sets, and **every** produced set is verified with
`verify.is_b3` (imported, not reimplemented).

The sequence produced is the classical greedy B_3 sequence
**OEIS A051912**: `1, 2, 5, 14, 33, 72, 125, 219, 376, 573, 745, 1209, 1557,
2442, 3098, 4048, 5298, 6704, 7839, 10987, 12332, 15465, ...` — the B_3 analogue
of Mian–Chowla (which is the greedy B_2/Sidon sequence).

## Numbers

Requested grid (`N`, `|A|`, ratio at `Nmax = max(S)`), all verified B_3:

| N    | Nmax | \|A\| | ratio  |
|------|------|-------|--------|
| 50   | 33   | 5     | 1.5588 |
| 100  | 72   | 6     | 1.4422 |
| 200  | 125  | 7     | 1.4000 |
| 500  | 376  | 9     | 1.2469 |
| 1000 | 745  | 11    | 1.2134 |
| 2000 | 1557 | 13    | 1.1216 |
| 5000 | 4048 | 16    | 1.0039 |

Monotone decay. Continuing past the requested grid, it crosses below 1 and
keeps dropping:

| N (ceiling) | Nmax    | \|A\| | ratio  |
|-------------|---------|-------|--------|
| 10,000      | 7,839   | 19    | 0.9565 |
| 30,000      | 28,974  | 25    | 0.8140 |
| 100,000     | 88,545  | 33    | 0.7404 |
| 300,000     | 269,154 | 42    | 0.6505 |
| 1,000,000   | 942,493 | 55    | 0.5610 |

## Asymptotics

Log-log fit over `Nmax in [10^3, 10^6]`:

    ratio ~ Nmax^{-0.106}   (slow power-law decay to 0)
    |A|   ~ Nmax^{0.228}     (NOT Nmax^{1/3})
    a_k   ~ k^{4.4}          (greedy terms grow faster than k^3)

So `|A| / Nmax^{1/3} -> 0`. The "beats 1" at N <= 5000 is a pure finite-size
effect from the small leading terms {1,2,5,14}; it is not sustained.

## Comparison to Bose–Chowla at matched scale

Bose–Chowla -> 1 from above and stays essentially pinned at 1:

| q   | N        | \|A\| | ratio    |
|-----|----------|-------|----------|
| 11  | 1,330    | 11    | 1.00025  |
| 16  | 4,095    | 16    | 1.00008  |
| 64  | 262,143  | 64    | 1.000003 |
| 101 | 1,030,300| 101   | 1.000001 |

At N ≈ 10^6: Bose–Chowla has |A| = 101 (ratio ≈ 1.000), greedy has |A| = 55
(ratio ≈ 0.561). Greedy gives **roughly half** as many elements per cube-root
window and the gap widens with N.

## Takeaway for the question

Greedy is the natural baseline and it is the WRONG direction: it is thin, not
dense. It does not beat Bose–Chowla asymptotically, let alone sustain `1+c`.
Any family that hopes to exceed `1+c` for fixed `c>0` must be *structured*
(algebraic / projective-plane type), not greedy. The greedy curve is a useful
floor to plot structured candidates against: a structured family is only
interesting if its ratio stays flat (or rising) where greedy is already
visibly collapsing (by N ≈ 10^4 greedy is below 1 and falling like N^{-0.1}).
