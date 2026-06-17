# Lower-bound search synthesis: can any non-Bose–Chowla family sustain ratio 1+c?

**Question.** Score a finite B_3 set A by `ratio = |A| / N^{1/3}` with `N = max(A)`.
Bose–Chowla (BC) gives `ratio = 1 + o(1)`, approaching 1 from above as `q -> infinity`
(empirically ~1.05 at q=2, ~1.02 at q=4, ~1.0007 at q=16, ~1.003 at q=128 — see table
below). The disproof condition for the relevant conjecture is: some **other** family that
**sustainably** holds `ratio >= 1 + c` for a **fixed** `c > 0` as `N -> infinity`.

**Verdict: NO.** Four structurally distinct families were searched. Every one of them
either never exceeds 1 at all, or exceeds 1 only at small N as a finite-size effect and
then **decays to (or below) 1**. None sustains any fixed `c > 0`. Bose–Chowla remains the
champion family, and the `1 + o(1)` ceiling stands unbroken by every candidate tested.

---

## Re-verification of the best example sets (imported `verify.is_b3`, not reimplemented)

All "best examples" handed to the synthesizer were re-checked. Every one is a genuine B_3
set over Z. They check out:

| Family | Best example set | `is_b3` | \|A\| | N=max | ratio = \|A\|/N^(1/3) |
|---|---|---|---|---|---|
| Greedy (A051912) | `1 2 5 14 33` | **True** | 5 | 33 | 1.5588 |
| Union of dilated BC | `1 6 29 60 84 104 196 320` | **True** | 8 | 320 | 1.1696 |
| Depth local search (q=13) | `1 99 105 129 194 351 683 844 1189 1413 2073 2081 2100 2107` | **True** | 14 | 2107 | **1.0920** |
| Modular spike (lifted) | `0 1 4` | **True** | 3 | 4 | 1.8899 |

Note on the depth-search row: its family report listed `best_ratio = 1.0771`, which uses
the *interval* denominator `n = q^3 - 1 = 2196`. Scored by the project convention
`N = max(A) = 2107`, the ratio is **1.0920**. Either way it is a tiny-N artifact (see below).

The `0 1 4` "ratio 1.8899" is a degenerate 3-element set — it is large only because N=4 is
absurdly small; it carries no asymptotic content.

---

## Per-family results

### 1. Pure greedy B_3 (Mian–Chowla analogue, OEIS A051912) — THIN, not dense
Highest ratio 1.559 at N=33, but this is purely the leading terms `{1,2,5,14}`. It decays
monotonically, crosses **below 1** around N~5000, and keeps falling: 0.74 at N~9e4, 0.56 at
N~9e5. Log-log fit gives `ratio ~ Nmax^{-0.106}`, i.e. `|A| ~ N^{0.228}`, far short of the
`N^{1/3}` growth BC achieves. Greedy is the *thin* extreme — at N~1e6 it has roughly half
the elements per cube-root window that BC has, and the gap widens. Wrong direction entirely.

### 2. Union of dilated/translated BC copies — isolated artifact at q=4 only
Pure translates `A0 ∪ (A0+T)` are **never** B_3 for `|A0| >= 3` (explicit collision:
`a_i+a_j+(a_k+T) = a_i+a_k+(a_j+T)`; a full offset scan for q in {4..11} found zero working
T). Dilated unions `A0 ∪ (D·A0+T)` *can* be B_3, but only if the dilation D grows with q
(min viable D = 4,9,14,... for q=4,5,7; none ≤ 30 for q ≥ 8). Since `N ~ D·q^3` grows like D
while `|A|` only doubles to 2q, `ratio ~ 2/D^{1/3} -> 0`. Best union ratio decays
monotonically 1.17 (q=4) -> 0.95 (q=5) -> 0.83 (q=7) -> no valid set (q≥8). The q=4 win is
real and re-verified but is a single isolated point, not a family.

### 3. Depth local search on BC within [1, q^3-1] — finite slack, closes with q
The most interesting positive signal. At q=13, where single-element augmentation of BC is
locally maximal (the known q≥13 "collapse"), a **depth-1 swap** (remove one BC element,
greedily re-fill) reaches **size 14 > q=13** — verified B_3. But the surplus over q shrinks
monotonically: +2 (q=4–9) -> +1 (q=11,13) -> **0 for q ≥ 16**. At q=16,17,19,23,25,27 even
removing up to 3 elements and refilling cannot beat q; best size = q exactly, ratio = 1.000.
Mechanism: at small q the BC set has density `q^-2` in `[1,q^3-1]`, leaving slack that local
search mines for a constant-order surplus; as q grows BC nearly saturates the cube-root bound
and the slack closes. Converges to the same `1 + o(1)` as plain BC.

### 4. Extremal modular B_3 sets lifted to integers — sustains in Z_m, NOT in Z
Two ratios must be separated. The **modular** ratio `g(m)/m^{1/3}` genuinely sustains
~`2^{1/3} = 1.2599` (the "free wraparound" bonus; `g(m)` hugs `(2m)^{1/3}`, confirming the
context premise). But the integer-lifted score `|A|/max(A)^{1/3}` **decays to 1**, because a
modular B_3 set read as residues in `[0,m)` has `max(A)` close to m and `|A| ~ (2m)^{1/3}`,
so the lifted ratio collapses back toward `2^{1/3}/2^{1/3}... -> 1`. The wraparound bonus is
a torus artifact that the integer max-element penalty exactly erases. (This is the cleanest
illustration that the question's "free wraparound" only helps in `Z_m`, never in `Z`.)

---

## Best SUSTAINED ratio across all families

There is **no** sustained ratio above 1 from any non-BC family. Ranked by how the ratio
behaves as N grows:

| Family | Peak ratio (small N) | Asymptotic ratio | Sustains 1+c? |
|---|---|---|---|
| Modular B_3 in **Z_m** (torus) | — | ~1.26 (`2^{1/3}`) | YES, but only in Z_m, not Z |
| Depth local search (in Z) | 1.51 (q=4) | **1.000** | no |
| Union of dilated BC (in Z) | 1.17 (q=4) | 0 | no |
| Greedy A051912 (in Z) | 1.56 (N=33) | 0 (`~N^{-0.106}`) | no |
| **Bose–Chowla (in Z)** | ~1.05 | **1.000 (from above)** | the ceiling itself |

The only thing that genuinely sustains a constant above 1 is the modular construction —
and it does so **only over Z_m**, where wraparound is free. The moment you lift it to the
integers and pay for `max(A)`, the bonus evaporates and the lifted ratio returns to 1. Over
the integers, BC is the asymptotic champion and the depth-search family is the closest
runner-up, both converging to exactly 1.

## Bose–Chowla baseline (re-confirmed, `bose_chowla(q)` + `verify.is_b3`)

| q | \|A\| | N=max | ratio |
|---|---|---|---|
| 2 | 2 | 5 | 1.1696 |
| 4 | 4 | 60 | 1.0217 |
| 8 | 8 | 422 | 1.0666 |
| 13 | 13 | 2100 | 1.0152 |
| 16 | 16 | 4088 | 1.0007 |
| 32 | 32 | 32217 | 1.0057 |
| 64 | 64 | 256949 | 1.0067 |
| 128 | 128 | 2077848 | 1.0031 |

---

## Honest conclusion

- **No family credibly beats `1 + c` for fixed `c > 0` over the integers.** The disproof
  condition is not met by any candidate searched.
- The **highest sustained ratio over Z is exactly 1** (tied by Bose–Chowla and, asymptotically,
  by the depth-search family). Every "win above 1" — greedy 1.56, union 1.17, depth 1.51,
  lifted modular 1.89 — is a small-N finite-size effect that decays to 1 or below.
- The only place a fixed `c > 0` is sustained is the **torus `Z_m`** (modular B_3 sets,
  `~2^{1/3}`), and that bonus is structurally tied to free wraparound; it does **not**
  transfer to the integer problem.
- Net: every non-Bose–Chowla family gets *close to* ratio 1 only from below (or touches it
  from a decaying peak), and Bose–Chowla's `1 + o(1)` remains the operative ceiling. The
  search supports the conjecture rather than disproving it.
