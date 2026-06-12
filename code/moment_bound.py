"""Exact derivation + validation of the modular moment-chain bound.

THEOREM CANDIDATE (modular):
  If A is a B_3 set (multiset convention) in a finite abelian group G of
  order m, with k = |A| >= 2, then
      m >= (k^4 - 3k^3 + 5k^2 - 4k) / (2k - 3),
  hence k^3 <= 2m + 3k^2 <= 2m + 3*(6m)^{2/3}.

Chain:
  (1) B_3 => Sidon (padding argument)
  (2) E_2 := #{(a,b,c,d): a+b=c+d}        = 2k^2 - k        (exact, Sidon)
      E_3 := #{6-tuples: a1+a2+a3=a4+a5+a6} = 6k^3 - 9k^2 + 4k (exact, B_3)
  (3) Parseval: sum_chi |Ahat(chi)|^4 = m E_2 ; sum_chi |Ahat(chi)|^6 = m E_3 ;
      sum_chi |Ahat(chi)|^2 = mk ; Ahat(triv) = k.
  (4) Cauchy-Schwarz over nontrivial chi: (S2')^2 <= S1' * S3'
      where Sj' = sum_{chi nontrivial} |Ahat(chi)|^{2j}.

This script:
  A. symbolically expands (4) into the polynomial inequality and verifies the
     claimed closed forms (sympy, exact);
  B. proves the corollary chain k^3 <= 2m + 3k^2 (k>=2) exactly;
  C. numerically validates E_2/E_3 exactness + Parseval on Bose-Chowla sets
     (integer convolution counts, no floats in the identities);
  D. tabulates the bound vs the counting bound vs exact g(m) data;
  E. verifies the moment-problem optimality remark: min E[y^3] s.t.
     E[y]=1, E[y^2]=2, y>=0 equals 4 (so constant 2 is the best this
     constraint set can give).
"""

import sys
import numpy as np
import sympy as sp

k, m = sp.symbols('k m', positive=True)

print("=" * 70)
print("A. Symbolic expansion of the Cauchy-Schwarz inequality")
print("=" * 70)

E2 = 2 * k**2 - k
E3 = 6 * k**3 - 9 * k**2 + 4 * k
S1 = m * k - k**2          # sum over nontrivial chi of |Ahat|^2
S2p = m * E2 - k**4        # ... of |Ahat|^4
S3p = m * E3 - k**6        # ... of |Ahat|^6

P = sp.expand(S1 * S3p - S2p**2)   # C-S requires P >= 0
print("P(m,k) = S1*S3' - S2'^2 =", P)

# claim: P = m * (alpha*m + beta) with
alpha = k**2 * (2 * k - 3) * (k - 1)
beta = -(k**7 - 4 * k**6 + 8 * k**5 - 9 * k**4 + 4 * k**3)
check = sp.expand(P - m * (alpha * m + beta))
print("P - m*(alpha*m + beta) =", check, " (must be 0)")
assert check == 0

mstar = sp.cancel(-beta / alpha)
mstar_simpl = sp.cancel(sp.factor(-beta) / sp.factor(alpha))
print("P >= 0  <=>  m >= -beta/alpha =", sp.factor(-beta), "/", sp.factor(alpha))
print("          =", mstar_simpl)
target = (k**4 - 3 * k**3 + 5 * k**2 - 4 * k) / (2 * k - 3)
assert sp.simplify(mstar - target) == 0, "closed form mismatch"
print("Verified: m >= (k^4 - 3k^3 + 5k^2 - 4k)/(2k - 3)   for k >= 2 (alpha>0).")

print()
print("=" * 70)
print("B. Corollary  k^3 <= 2m + 3k^2  (k >= 2), and k <= (2m)^(1/3) + 2")
print("=" * 70)
# (k^4-3k^3+5k^2-4k)/(2k-3) >= (k^3-3k^2)/2  <=>  2(k^4-3k^3+5k^2-4k) - (k^3-3k^2)(2k-3) >= 0
diff = sp.expand(2 * (k**4 - 3 * k**3 + 5 * k**2 - 4 * k) - (k**3 - 3 * k**2) * (2 * k - 3))
print("2*mstar*(2k-3)/... difference polynomial:", diff)
# diff = 3k^3 + k^2 - 8k = k(3k^2 + k - 8) > 0 for k >= 2
assert sp.simplify(diff - (3 * k**3 + k**2 - 8 * k)) == 0
print("= k(3k^2+k-8) > 0 for k>=2  =>  m >= (k^3 - 3k^2)/2  =>  k^3 <= 2m + 3k^2. OK")

# sharper integer corollary k <= ceil((2m)^(1/3)) + 1 style: test numerically below.

print()
print("=" * 70)
print("C. Numerical validation of E_2 / E_3 exactness + Parseval (Bose-Chowla)")
print("=" * 70)
sys.path.insert(0, '.')
from bose_chowla import bose_chowla  # noqa: E402
from verify import is_b3             # noqa: E402

for q in [3, 4, 5, 7, 8, 9]:
    A, n = bose_chowla(q)
    mm = n
    kk = len(A)
    assert is_b3(A, modulus=mm)
    # exact integer energies via convolution counting
    r2 = np.zeros(mm, dtype=np.int64)
    for a in A:
        for b in A:
            r2[(a + b) % mm] += 1
    E2_actual = int((r2.astype(object) ** 2).sum())
    r3 = np.zeros(mm, dtype=np.int64)
    for a in A:
        for b in A:
            for c in A:
                r3[(a + b + c) % mm] += 1
    E3_actual = int((r3.astype(object) ** 2).sum())
    e2_pred = 2 * kk**2 - kk
    e3_pred = 6 * kk**3 - 9 * kk**2 + 4 * kk
    # Parseval check (floats, but only as a redundant cross-check)
    fa = np.fft.fft(np.bincount(np.array(A) % mm, minlength=mm))
    y = np.abs(fa) ** 2
    par4 = y.sum() ** 0 * (y**2).sum()  # sum |Ahat|^4
    par6 = (y**3).sum()
    ok2 = E2_actual == e2_pred
    ok3 = E3_actual == e3_pred
    okp4 = abs(par4 - mm * E2_actual) < 1e-3 * mm * E2_actual
    okp6 = abs(par6 - mm * E3_actual) < 1e-3 * mm * E3_actual
    # Cauchy-Schwarz slack on nontrivial modes
    ynz = np.sort(y)[:-1]  # drop the largest = k^2 at chi=0... safer: remove index 0
    ynz = np.delete(y, 0)
    s1, s2p, s3p = ynz.sum(), (ynz**2).sum(), (ynz**3).sum()
    cs_ratio = s1 * s3p / s2p**2
    print(f"q={q:2d} m={mm:5d} k={kk:2d}  E2 exact:{ok2}  E3 exact:{ok3}  "
          f"Parseval4:{okp4} Parseval6:{okp6}  CS-slack S1*S3'/S2'^2 = {cs_ratio:.4f}")
    assert ok2 and ok3 and okp4 and okp6

print()
print("=" * 70)
print("D. Bound table: smallest m admitting a size-k B_3 set in Z_m")
print("=" * 70)
print(" k   counting C(k+2,3)<=m   moment bound m>=ceil(...)   ratio k^3/m_bound")
for kk in range(2, 26):
    mc = sp.binomial(kk + 2, 3)
    mb = sp.ceiling(sp.Rational(kk**4 - 3 * kk**3 + 5 * kk**2 - 4 * kk, 2 * kk - 3))
    print(f"{kk:3d}   {int(mc):8d}              {int(mb):8d}              {kk**3 / int(mb):.4f}")

print()
print("=" * 70)
print("E. Moment-problem optimality: min E[y^3] s.t. E[y]=1, E[y^2]=2, y>=0")
print("=" * 70)
# Two-point candidate {0, a} with mass p at a: pa=1, pa^2=2 => a=2, p=1/2, E[y^3]=4.
# Lower bound: for all y>=0: y^3 - lambda1*y - lambda2*y^2 - lambda0 >= 0 with
# lambda0 + lambda1*1 + lambda2*2 = 4.  Take y^3 >= 3y^2 - ... find multipliers:
# y^3 - (3y^2 - 2y - ...)?  Use: y^3 + 4y >= 4y^2 - 0? <=> y(y-2)^2 >= 0. TRUE.
y = sp.symbols('y', nonnegative=True)
ident = sp.factor(y**3 - (4 * y**2 - 4 * y))
print("y^3 - 4y^2 + 4y =", ident, ">= 0 for y >= 0")
print("=> E[y^3] >= 4E[y^2] - 4E[y] = 8 - 4 = 4, matched by y ~ {0:1/2, 2:1/2}.")
print("So within {E[y]=1, E[y^2]=2, y>=0} the constant 2 in k^3<=2m is optimal.")
