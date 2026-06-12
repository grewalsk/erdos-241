"""Bose-Chowla B_3 sets, h = 3, for prime powers q.

Construction (Bose-Chowla 1962): let theta generate GF(q^3)^*.  Put
    A = { i in [1, q^3-1] : theta^i = theta + c for some c in GF(q) }.
Then |A| = q and A is a B_3 set modulo q^3 - 1 (hence also as integers).

Proof idea (for the writeup): if a1+a2+a3 = b1+b2+b3 (mod q^3-1) then
prod(theta+c_aj) = prod(theta+c_bj) in GF(q^3); both sides are values at
theta of monic cubics over GF(q); their difference is a poly of degree <= 2
vanishing at theta, which has degree 3 over GF(q); so the cubics are equal,
so {c_a} = {c_b} as multisets, so {a_j} = {b_j}.

Implementation: hand-rolled GF(p)[x]/(f) arithmetic, f irreducible of degree
3e where q = p^e; the subfield GF(q) is {0} union the powers of
g^((p^{3e}-1)/(q-1)) for a generator g; theta = g.  Discrete logs by a single
walk through all powers of g.  EVERY output is checked by the standalone
brute-force verifier (verify.is_b3), both mod q^3-1 and as integers.

Also: FFT-based exact augmentation analysis -- how many x in [1, q^3-1]
can be added to the Bose-Chowla set keeping B_3 (as integers), and greedy
augmentation counts.
"""

import sys
import numpy as np
from sympy import factorint, isprime
from verify import is_b3

sys.setrecursionlimit(10000)


# ---------- GF(p)[x]/(f) arithmetic, coefficients as tuples (low degree first)

def poly_mul_mod(a, b, f, p):
    d = len(f) - 1  # f monic, degree d
    res = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                res[i + j] = (res[i + j] + ai * bj) % p
    # reduce mod f (monic)
    for i in range(len(res) - 1, d - 1, -1):
        c = res[i]
        if c:
            res[i] = 0
            for j in range(d):
                res[i - d + j] = (res[i - d + j] - c * f[j]) % p
    return tuple(res[:d] + [0] * (d - len(res)))[:d] if len(res) < d else tuple(res[:d])


def poly_pow_mod(a, e, f, p):
    d = len(f) - 1
    result = tuple([1] + [0] * (d - 1))
    base = a
    while e:
        if e & 1:
            result = poly_mul_mod(result, base, f, p)
        e >>= 1
        if e:
            base = poly_mul_mod(base, base, f, p)
    return result


def _poly_norm(a, p):
    a = [c % p for c in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def _poly_divmod(a, b, p):
    a = _poly_norm(a, p)
    b = _poly_norm(b, p)
    binv = pow(b[-1], p - 2, p)
    q = [0] * max(1, len(a) - len(b) + 1)
    r = list(a)
    while len(r) >= len(b) and any(r):
        c = (r[-1] * binv) % p
        s = len(r) - len(b)
        q[s] = c
        for i, bi in enumerate(b):
            r[s + i] = (r[s + i] - c * bi) % p
        r = _poly_norm(r, p)
        if not r:
            break
    return q, r


def _poly_gcd(a, b, p):
    a = _poly_norm(a, p)
    b = _poly_norm(b, p)
    while b:
        _, r = _poly_divmod(a, b, p)
        a, b = b, r
    return a


def is_irreducible(f, p):
    """Rabin's irreducibility test for monic f of degree d over F_p:
    x^(p^d) == x (mod f), and gcd(x^(p^(d/r)) - x, f) = 1 for primes r | d.
    (The gcd step is essential for composite d: without it, a product of
    irreducible factors of degrees e.g. {1,2,3} passes for d = 6.)"""
    d = len(f) - 1
    x = tuple([0, 1] + [0] * (d - 2)) if d >= 2 else (0,)
    xp = poly_pow_mod(x, p ** d, f, p)
    if xp != x:
        return False
    for r in set(factorint(d)):
        xr = poly_pow_mod(x, p ** (d // r), f, p)
        diff = [(xi - yi) % p for xi, yi in zip(xr, x)]
        g = _poly_gcd(list(f), diff, p)
        if len(g) - 1 != 0:  # gcd not constant => reducible
            return False
    return True


def find_irreducible(d, p, rng):
    while True:
        coeffs = [rng.randrange(p) for _ in range(d)] + [1]
        if coeffs[0] == 0:
            coeffs[0] = 1 + rng.randrange(p - 1) if p > 1 else 1
        if is_irreducible(coeffs, p):
            return coeffs


def bose_chowla(q, seed=0):
    """Return (A, n) with A a Bose-Chowla B_3 set of size q in [1, n], n = q^3-1."""
    import random
    rng = random.Random(seed)
    fac = factorint(q)
    assert len(fac) == 1, f"{q} not a prime power"
    p, e = next(iter(fac.items()))
    d = 3 * e
    n = q ** 3 - 1            # = p^(3e) - 1
    f = find_irreducible(d, p, rng)
    one = tuple([1] + [0] * (d - 1))
    # find generator g of GF(p^d)^*
    prim_divs = list(factorint(n))
    while True:
        g = tuple(rng.randrange(p) for _ in range(d))
        if all(c == 0 for c in g):
            continue
        if all(poly_pow_mod(g, n // r, f, p) != one for r in prim_divs):
            break
    # subfield GF(q): {0} + powers of g^(n/(q-1))
    sub = {tuple([0] * d)}
    if q > 2 or True:
        h = poly_pow_mod(g, n // (q - 1), f, p) if q > 1 else one
        cur = one
        for _ in range(q - 1):
            sub.add(cur)
            cur = poly_mul_mod(cur, h, f, p)
    assert len(sub) == q, (len(sub), q)
    # targets: {g + c : c in GF(q)}  ->  collect discrete logs by one walk
    targets = set()
    for c in sub:
        t = tuple((gi + ci) % p for gi, ci in zip(g, c))
        targets.add(t)
    assert len(targets) == q
    A = []
    cur = g
    for i in range(1, n + 1):
        if cur in targets:
            A.append(i)
        cur = poly_mul_mod(cur, g, f, p)
    assert len(A) == q, f"got {len(A)} logs, expected {q}"
    return sorted(A), n


# ---------- exact augmentation analysis (integers in [1, n]) ----------

def addable_mask(A, N):
    """Boolean mask over [0..N]: addable[x] = True iff A + {x} is still B_3
    (as integers in [1,N]), for x not in A.  Exact, via integer FFT
    cross-correlations + exact small checks.  A must be B_3 (hence Sidon)."""
    A = sorted(A)
    k = len(A)
    L = 3 * N + 4
    aind = np.zeros(L, dtype=np.float64)
    for a in A:
        aind[a] = 1.0
    s2 = np.zeros(L, dtype=np.float64)
    s3 = np.zeros(L, dtype=np.float64)
    pair_sums = set()
    for i in range(k):
        for j in range(i, k):
            pair_sums.add(A[i] + A[j])
    triple_sums = set()
    for i in range(k):
        for j in range(i, k):
            for l in range(j, k):
                triple_sums.add(A[i] + A[j] + A[l])
    for s in pair_sums:
        s2[s] = 1.0
    for s in triple_sums:
        s3[s] = 1.0
    # cond1: exists a in A with x+a in S2  <=> corr1[x] > 0
    # cond2: exists p in S2 with x+p in S3 <=> corr2[x] > 0
    F = 1 << int(np.ceil(np.log2(2 * L)))
    fa = np.fft.rfft(aind, F)
    f2 = np.fft.rfft(s2, F)
    f3 = np.fft.rfft(s3, F)
    corr1 = np.fft.irfft(np.conj(fa) * f2, F)[:L]   # corr1[x] = sum_a s2[x+a]... careful
    corr2 = np.fft.irfft(np.conj(f2) * f3, F)[:L]
    # np convention: irfft(conj(F(a))*F(b))[x] = sum_t a[t] b[t+x]  (circular);
    # our arrays are zero-padded to F >= 2L so no wraparound for 0<=x<L.
    bad = np.zeros(N + 1, dtype=bool)
    xs = np.arange(N + 1)
    bad |= corr1[: N + 1] > 0.5                     # x+a in S2
    bad |= corr2[: N + 1] > 0.5                     # x+p in S3 for p in S2
    twox = 2 * xs
    bad |= s2[twox] > 0.5                            # 2x in S2
    bad |= s3[3 * xs] > 0.5                          # 3x in S3
    # 2x + a in S3 for some a: small loop (k terms)
    for a in A:
        bad |= s3[twox + a] > 0.5
    for a in A:
        bad[a] = True                                # already in the set
    bad[0] = True
    return ~bad


def greedy_augment(A, N, order="asc", cand=None):
    """Greedily add addable elements (re-deriving addability exactly each time
    with set arithmetic); returns augmented set.

    cand: optional candidate list. Passing the initially-addable elements is
    lossless for greedy: adding elements only adds pair/triple sums, so an x
    not addable to A is not addable to any superset of A (antitone)."""
    A = sorted(A)
    pair_sums = set()
    triple_sums = set()
    k = len(A)
    for i in range(k):
        for j in range(i, k):
            pair_sums.add(A[i] + A[j])
    for i in range(k):
        for j in range(i, k):
            for l in range(j, k):
                triple_sums.add(A[i] + A[j] + A[l])

    def ok(x, S):
        for a in S:
            if x + a in pair_sums:
                return False
        if 2 * x in pair_sums:
            return False
        for p in pair_sums:
            if x + p in triple_sums:
                return False
        for a in S:
            if 2 * x + a in triple_sums:
                return False
        if 3 * x in triple_sums:
            return False
        return True

    def add(x, S):
        newt = {x + p for p in pair_sums} | {2 * x + a for a in S} | {3 * x}
        newp = {x + a for a in S} | {2 * x}
        triple_sums.update(newt)
        pair_sums.update(newp)
        S.append(x)

    S = list(A)
    if cand is None:
        cand = [x for x in range(1, N + 1) if x not in set(A)]
    if order == "desc":
        cand = list(cand)[::-1]
    added = []
    for x in cand:
        if ok(x, S):
            add(x, S)
            added.append(x)
    return sorted(S), added


if __name__ == "__main__":
    qs = [int(x) for x in sys.argv[1:]] or [2, 3, 4, 5, 7, 8, 9, 11, 13, 16, 19, 23, 25, 27]
    print("q,n=q^3-1,|A|,ratio=|A|/n^(1/3),verified_mod,verified_int,addable_count,greedy_final,greedy_ratio")
    for q in qs:
        A, n = bose_chowla(q)
        okm = is_b3(A, modulus=n)
        oki = is_b3(A)
        assert okm and oki, f"VERIFIER FAILED for q={q}"
        mask = addable_mask(A, n)
        nadd = int(mask.sum())
        candidates = [int(x) for x in np.where(mask)[0]]
        S, added = greedy_augment(A, n, cand=candidates)
        assert is_b3(S), f"augmented set fails verifier for q={q}!"
        ratio = len(A) / n ** (1 / 3)
        gratio = len(S) / n ** (1 / 3)
        print(f"{q},{n},{len(A)},{ratio:.6f},{okm},{oki},{nadd},{len(S)},{gratio:.6f}", flush=True)
