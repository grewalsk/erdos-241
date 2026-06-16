"""Green's variational constant gamma(p), reproduction + convex optimization.

Green (Acta Arith. 100 (2001), 365-390) proves, for B_3 sets A in {1..N}:
    |A| <= (4/(1+gamma))^{1/3} N^{1/3} (1+o(1)),
for ANY real p in C^1[0,1] with int_0^1 p = 2, where
    gamma(p) = 2 * ( sum_{r>=1} |c_r(p)|^{4/3} )^{-3},
    c_r(p)   = int_0^1 p(x) e^{-i pi r x} dx.
His choice p(x) = 5/2 - 40(x-1/2)^4 gives gamma ~ 1/6.9994, i.e. the constant 7/2.

J(p) := sum_{r>=1} |c_r(p)|^{4/3} is CONVEX in p (sum of convex functions of
linear functionals), the constraint is affine, so the global infimum is
computable.  This script:
  1. reproduces Green's S1, S2, gamma exactly from (26)-(29) of the paper;
  2. evaluates J for several closed-form candidate kernels;
  3. solves the discretized convex program min J(p) numerically;
  4. reports the implied B_3 constant 4/(1+gamma*).

Frequency split: for r = 2s even, c_r = standard circle Fourier coefficient
p_hat(s); for r = 2s+1 odd, c_r = (p(x)e^{-i pi x})_hat(s).
"""

import numpy as np
from scipy.optimize import minimize

# ---------------------------------------------------------------- reproduction
def green_reproduction():
    # (26): |p~(pi r)| = 40/pi^2 |1/r^2 - 24/(pi^2 r^4)|   (r even, r != 0)
    #                   = 40/pi^2 * (6/pi) |1/r^3 - 8/(pi^2 r^5)|   (r odd)
    R = 4_000_000
    r = np.arange(1, R + 1, dtype=np.float64)
    even = r[1::2]  # 2,4,6,...
    odd = r[0::2]   # 1,3,5,...
    S1 = np.sum(np.abs(1 / even**2 - 24 / (np.pi**2 * even**4)) ** (4 / 3))
    S2 = np.sum(np.abs(1 / odd**3 - 8 / (np.pi**2 * odd**5)) ** (4 / 3))
    # tails: |1/r^2 - ...|^{4/3} <= r^{-8/3}; integral comparison
    tail1 = 3 / 5 * (R ** (-5 / 3))   # sum_{r>R} r^{-8/3} <= int_R^inf = (3/5) R^{-5/3}
    tail2 = 1 / 3 * (R ** (-3))       # sum_{r>R} r^{-4} <= (1/3) R^{-3}
    gamma = 2 * (np.pi**2 / 40) ** 4 / (S1 + (6 / np.pi) ** (4 / 3) * S2) ** 3
    J = (40 / np.pi**2) ** (4 / 3) * (S1 + (6 / np.pi) ** (4 / 3) * S2)
    print(f"S1 = {S1:.7f}   (Green: 0.0839757)  [tail < {tail1:.2e}]")
    print(f"S2 = {S2:.7f}   (Green: 0.1219299)  [tail < {tail2:.2e}]")
    print(f"gamma(p_Green) = {gamma:.7f} = 1/{1/gamma:.4f}   (Green: 1/6.9994)")
    print(f"J(p_Green) = {J:.6f}")
    print(f"B_3 constant 4/(1+gamma) = {4/(1+gamma):.6f}  (Green: 3.5)")
    return J, gamma


# ------------------------------------------------------- generic J evaluation
def J_of_p(pvals, R=20000):
    """J(p) for p given as values on the midpoint grid x_j=(j+1/2)/n.
    Midpoint quadrature for c_r; adequate for exploration (certification of
    any winner is done separately in exact arithmetic)."""
    n = len(pvals)
    x = (np.arange(n) + 0.5) / n
    # c_r for r = 1..R via FFT-like direct evaluation in blocks
    rs = np.arange(1, R + 1)
    # phase matrix would be R x n; do it blockwise
    total = 0.0
    B = 2000
    for lo in range(0, R, B):
        rblk = rs[lo:lo + B]
        ph = np.exp(-1j * np.pi * np.outer(rblk, x))
        c = ph @ pvals / n
        total += np.sum(np.abs(c) ** (4 / 3))
    return total


def gamma_of_J(J):
    return 2.0 / J**3


# ------------------------------------------------------------- convex program
def optimize_p(n=400, R=6000, iters=3, x0=None):
    """min J(p) over p on n-point grid with mean exactly 2 (projected)."""
    x = (np.arange(n) + 0.5) / n
    rs = np.arange(1, R + 1)
    PH = np.exp(-1j * np.pi * np.outer(rs, x))  # R x n

    def obj_grad(p):
        p = p - (p.mean() - 2.0)               # project onto mean-2 affine set
        c = PH @ p / n                          # R complex
        a = np.abs(c)
        f = np.sum(a ** (4 / 3))
        # d|c|^{4/3}/dc = (4/3)|c|^{1/3} * c/|c| ; chain rule through real p
        w = np.where(a > 1e-14, a ** (-2 / 3), 0.0) * (4 / 3)
        g = (np.conj(PH).T @ (w * c)).real / n
        g = g - g.mean()                        # gradient within the affine set
        return f, g

    if x0 is None:
        x0 = 2.5 - 40 * (x - 0.5) ** 4 + 0.0   # start near Green's kernel shape
        x0 = x0 - x0.mean() + 2.0
    res = minimize(obj_grad, x0, jac=True, method="L-BFGS-B",
                   options={"maxiter": 2000, "ftol": 1e-14, "gtol": 1e-12})
    p = res.x - (res.x.mean() - 2.0)
    return p, res.fun


if __name__ == "__main__":
    print("=" * 70)
    print("1. Reproduction of Green's numbers")
    print("=" * 70)
    Jg, gg = green_reproduction()

    print()
    print("=" * 70)
    print("2. Sanity: J for closed-form kernels (grid evaluation, R=20000)")
    print("=" * 70)
    n = 4000
    x = (np.arange(n) + 0.5) / n
    for name, vals in [
        ("constant 2", np.full(n, 2.0)),
        ("Green quartic 5/2-40(x-1/2)^4", 2.5 - 40 * (x - 0.5) ** 4),
        ("parabola 3-6(x-1/2)^2*2", 3 - 12 * (x - 0.5) ** 2),
    ]:
        J = J_of_p(vals, R=20000)
        g = gamma_of_J(J)
        print(f"{name:34s} J = {J:.5f}  gamma = 1/{1/g:.3f}  const = {4/(1+g):.5f}")

    print()
    print("=" * 70)
    print("3. Convex optimization (discretized), increasing resolution")
    print("=" * 70)
    best = None
    for n_, R_ in [(200, 4000), (400, 8000), (800, 12000)]:
        p, J = optimize_p(n=n_, R=R_)
        g = gamma_of_J(J)
        print(f"n={n_:4d} R={R_:6d}:  J* = {J:.6f}  gamma* = 1/{1/g:.4f}  "
              f"B3 const 4/(1+gamma*) = {4/(1+g):.6f}")
        best = (p, J, n_)
    p, J, n_ = best
    np.save("/Users/kabirgrewal/projects/erdos-241/data/popt.npy", p)
    print("optimal p saved to data/popt.npy; p range:",
          float(p.min()), "..", float(p.max()))
