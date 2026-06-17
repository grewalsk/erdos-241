# Modular B_3 sets in Z_m (m not q^3-1): do dense modular sets lift to dense integer sets?

**Family explored:** exact extremal B_3 sets in Z_m for moduli m NOT of the form
q^3-1, recomputed by brute force (gsolver, validated against `verify.is_b3` over
Z_m) for m up to 200. Question: when such a set is read as an integer set in
[0, m), is it B_3 over the **integers**, and does its integer score
`ratio = |A| / N^{1/3}` (N = max(A)) sustainably beat 1?

**Scoring & verification.** Every set below was checked with the standalone
`verify.is_b3` (imported, not reimplemented), both mod m and over Z. The
Bose-Chowla baseline ratio -> 1+ from above.

## Headline answer: NO sustained beat. The modular density bonus does NOT transfer to Z.

There are two distinct ratios, and conflating them is the trap:

| ratio | what it measures | behavior as the size grows |
|---|---|---|
| `g(m)/m^{1/3}` (modular) | density per modulus-scale | **genuinely sustains ~ 2^{1/3} = 1.2599** (the "free wraparound" bonus; g(m) hugs (2m)^{1/3}) |
| `|A|/max(A)^{1/3}` (integer-lifted) | the actual EP-241 integer score | **decays toward 1**, tracking the integer envelope f(N)/N^{1/3} -> 1 |

So the context premise is confirmed in the modular world (g(m) is denser per
m^{1/3} than integer sets, by a factor ~2^{1/3}) — but that bonus **evaporates on
lift to Z**.

## Why the bonus does not transfer (the mechanism)

A B_3 set that achieves the high modular density `g(m)/m^{1/3} ~ 1.26` does so by
**packing its elements into a SHORT interval** and letting modular wraparound
separate the triple-sums that would otherwise collide over Z. Concretely, the
densest extremal witnesses have `span = max(A)` much smaller than m:

```
m=30  g=4 span=19 span/m=0.63 ratio_m=1.287
m=13  g=3 span= 4 span/m=0.31 ratio_m=1.276
m=32  g=4 span=13 span/m=0.41 ratio_m=1.260
m=65  g=5 span=24 span/m=0.37 ratio_m=1.244
m=128 g=6 span=119 span/m=0.93 ratio_m=1.191  (closest to full span -> LOWER ratio_m)
```

The high-`ratio_m` sets have small `span/m`. When you read such a set over Z, N =
max(A) = span is *small*, so `|A|/N^{1/3}` is just the ordinary small-N integer
ratio — a finite-size effect, not a sustained gain. The only way the modular bonus
could carry over is if the densest sets used (nearly) full span (`span ~ m`); but
those are exactly the sets with the LOWEST modular ratio (see m=128 above). The
two objectives are in tension, and the integer objective always loses.

## Every extremal modular witness is trivially B_3 over Z — and that's not a win

For all 198 exact extremal witnesses up to m=200, the set (which contains 0 and has
`C(g+2,3) <= m`) is **also** B_3 over the integers. This is automatic: g(m) is so
small that the elements are tiny and no wraparound is actually load-bearing for the
*recorded* extremal witness. Reading it over Z just gives an ordinary (and in fact
**sub-optimal**) integer B_3 set. Head-to-head against the true integer optimum
f(N) at the same N = max(A):

- when the modular witness happens to be integer-optimal, it **ties** f(N);
- otherwise it is **strictly worse** ("no" in the table below);
- it **never** beats the integer optimum (it can't — f is the optimum).

```
 m  |A|  N   mod_ratio  f(N)  fN_ratio   verdict
30   4  19    1.4990     4    1.4990     ties (this IS an integer-optimal set)
32   4  13    1.7012     4    1.7012     ties
65   5  24    1.7334     5    1.7334     ties
69   5  53    1.3311     6    1.5973     WORSE than integer optimum
73   5  46    1.3955     6    1.6746     WORSE
```

I also enumerated ALL maximal modular B_3 sets (not just the gsolver witness) for
m = 4,13,30,32,65 and took the most generous integer scoring (N = span = max-min
after shifting min to 0). Best integer ratios: m=13 -> 1.890 (span 4), m=32 ->
1.747 (span 12, set {0,1,5,12}), m=65 -> 1.758 (span 23, set {0,1,15,18,23}). All
B_3 over Z, all just small-N integer sets matching the f(N) ceiling at their tiny
N. None give a sustained-`>1+c` family.

## Asymptotic envelope (the honest decay)

Peak integer-lifted ratio per modulus window (decaying):

```
m in [ 13, 40]: peak |A|/max(A)^{1/3} = 1.890  (m=13,  N=4,  |A|=3)
m in [ 41, 80]: peak                  = 1.733  (m=65,  N=24, |A|=5)
m in [ 81,130]: peak                  = 1.575  (m=85,  N=32, |A|=5)
m in [131,200]: peak                  = 1.587  (m=131, N=54, |A|=6)
```

These large numbers are entirely a small-N artifact: max(A) is 4, 24, 32, 54 —
tiny. They sit on the f(N)/N^{1/3} curve (f(50)/50^{1/3}=1.63, f(100)=1.51,
f(159)=1.48 ...), which itself -> 1. There is no modulus window in which the
integer-lifted ratio stays above 1 + c for fixed c as N -> infinity.

## Verdict for the structured output

- The modular density premise is TRUE (g(m)/m^{1/3} sustains ~1.26 = 2^{1/3}).
- It does NOT transfer to Z: reading extremal modular sets over Z gives ordinary,
  often sub-optimal, integer B_3 sets whose score decays to the f(N)/N^{1/3} -> 1
  envelope.
- Apparent integer ratios of 1.5-1.9 are finite-size effects (N = max(A) is tiny),
  exactly the mirage to be honest about. This family does NOT sustainably beat 1+c.
