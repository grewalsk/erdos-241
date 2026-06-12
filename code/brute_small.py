"""Ground-truth exact values of f(N) and g(m) for small arguments.

Deliberately dumb: the ONLY legality test is the standalone verifier
(verify.is_b3) applied to the whole candidate set at every extension step.
No incremental bitset logic, no shared code with the fast C solver.
Used to validate the C solver on small instances.

Pruning used (both provably safe):
  - remaining-elements bound:  |A| + (#candidates left) <= best  => cut
  - counting bound: C(k+2,3) <= 3N-2 for f, C(k+2,3) <= m for g
"""

import sys
from math import comb
from verify import is_b3


def f_exact(N):
    """Exact f(N), plus one witness, by verifier-driven DFS."""
    best = {"size": 0, "set": []}

    def ext(A, start):
        if len(A) > best["size"]:
            best["size"] = len(A)
            best["set"] = A[:]
        for x in range(start, N + 1):
            if len(A) + 1 + (N - x) <= best["size"]:
                break  # even taking everything from x on can't beat best
            cand = A + [x]
            if is_b3(cand):
                ext(cand, x + 1)

    ext([], 1)
    return best["size"], best["set"]


def g_exact(m):
    """Exact g(m) for B_3 sets in Z_m, plus witness.

    WLOG 0 in A for any nonempty set (translation preserves B_3 in Z_m:
    sums shift by 3t, a bijection on sums).  g(m) >= 1 always.
    """
    # max k allowed by counting: C(k+2,3) <= m
    kmax = 0
    while comb(kmax + 3, 3) <= m:
        kmax += 1
    best = {"size": 1, "set": [0]}

    def ext(A, start):
        if len(A) > best["size"]:
            best["size"] = len(A)
            best["set"] = A[:]
        if len(A) >= kmax + 1:
            return
        for x in range(start, m):
            if len(A) + 1 + (m - 1 - x) <= best["size"]:
                break
            cand = A + [x]
            if is_b3(cand, modulus=m):
                ext(cand, x + 1)

    ext([0], 1)
    return best["size"], best["set"]


if __name__ == "__main__":
    mode = sys.argv[1]  # "f" or "g"
    lo, hi = int(sys.argv[2]), int(sys.argv[3])
    for n in range(lo, hi + 1):
        if mode == "f":
            k, w = f_exact(n)
        else:
            k, w = g_exact(n)
        print(f"{n},{k},\"{' '.join(map(str, w))}\"", flush=True)
