"""The triple-autoconvolution sup-norm relaxation mu_3 (Workstream C.2/C.3).

mu_3 := inf { ||f*f*f||_infty : f >= 0, supported on [0,1], int f = 1 }.

TRANSFER (proved in writeup/upper_bound_program.md):  for a B_3 set A in {1..N},
|A| = k, the normalized smeared measure f has ||f*f*f||_infty <= 6N/k^3 * (1+o(1)),
while ||f*f*f||_infty >= mu_3 by definition.  Hence
        k^3 <= (6/mu_3) N (1 + o(1)),    constant (6/mu_3).
The point of this script: estimate mu_3 to decide whether (6/mu_3)^{1/3} can beat
Green's 1.519.  (Spoiler/known caution from the brief: the sup-norm relaxation is
LOSSY; Green's sharp route uses an L^{4/3} Fourier functional, not the sup norm.)

We bound mu_3 two ways:
  * upper bound: evaluate ||f*f*f||_infty for candidate f (uniform, and an
    optimized step density) -> gives an over-estimate of the achievable constant's
    denominator, i.e. a *lower* bound on the constant (6/mu_3)^{1/3} we could hope.
    Wait: larger mu_3 => smaller constant. So to show the route CANNOT beat 3.5 we
    need an UPPER bound on mu_3 (=> lower bound on the constant).  ||f*f*f||_inf for
    ANY specific f is an UPPER bound on mu_3 (inf <= any value).  Uniform f gives
    mu_3 <= 3/4.  Optimizing f only DECREASES ||f*f*f||_inf, lowering this upper
    bound on mu_3 further, which only RAISES the resulting constant.  So already
    mu_3 <= 3/4 gives constant >= (6/(3/4))^{1/3} = 8^{1/3} = 2.0 — already proving
    the sup-norm route cannot reach Green's 1.519 (since 2.0 > 1.519).
"""
import numpy as np
from scipy.optimize import minimize


def conv3_max(f, n):
    """||f*f*f||_inf for step density f on n bins of width 1/n (sum f/n = 1)."""
    g = np.convolve(np.convolve(f, f), f) / (n * n)   # density values at resolution 1/n
    return g.max()


def uniform_mu3():
    # f = 1 on [0,1]; f*f*f is the order-3 B-spline; peak = 3/4 at center.
    for n in [100, 1000, 10000]:
        f = np.ones(n)              # value 1 on each bin (density 1)
        print(f"  uniform, n={n}: ||f*f*f||_inf = {conv3_max(f, n):.6f}  (analytic 0.75)")


def optimize_mu3(n=120):
    """Minimize ||f*f*f||_inf over step densities (only lowers the upper bound on mu_3)."""
    def obj(z):
        f = np.abs(z)
        s = f.sum() / n
        f = f / s
        return conv3_max(f, n)
    best = None
    for seed in range(4):
        rng = np.random.RandomState(seed)
        z0 = np.ones(n) + 0.1 * rng.randn(n)
        res = minimize(obj, z0, method="Nelder-Mead",
                       options={"maxiter": 60000, "xatol": 1e-6, "fatol": 1e-9, "adaptive": True})
        if best is None or res.fun < best:
            best = res.fun
    return best


if __name__ == "__main__":
    print("mu_3 = inf ||f*f*f||_inf over densities on [0,1]")
    print("=" * 60)
    print("Specific densities give UPPER bounds on mu_3:")
    uniform_mu3()
    mlo = optimize_mu3(n=120)
    print(f"  optimized step density (n=120): ||f*f*f||_inf = {mlo:.5f}  => mu_3 <= {mlo:.5f}")
    print()
    print("Resulting upper-bound constant for B_3 via this route: (6/mu_3)^(1/3)")
    for mu in [0.75, mlo]:
        print(f"  mu_3 <= {mu:.5f}  =>  constant >= (6/{mu:.5f})^(1/3) = {(6/mu)**(1/3):.4f} N^(1/3)")
    print()
    print("CONCLUSION: even the most favorable (smallest) mu_3 makes the constant")
    print(">= 2.0 N^(1/3), FAR above Green's 1.519.  The sup-norm autoconvolution")
    print("relaxation is too lossy for h=3; the working route is Green's L^{4/3}")
    print("Fourier functional (see green_optimize.py / green_certify2.py).")
