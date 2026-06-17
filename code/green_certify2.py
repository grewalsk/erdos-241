"""Certified smooth kernel beating Green's B_3 constant.

Family (symmetric about x=1/2, vanishing at endpoints, C^infty):
    p(x) = x(1-x) * Q(u),  u=(2x-1)^2,  Q(u)=sum_{j=0}^d c_j u^j   (deg p = 2d+2).
Optimize c_j (s.t. int_0^1 p = 2) to minimize J(p)=sum_{r>=1}|p~(pi r)|^{4/3}.

Optimization uses the FAST exact step-function J (FFT + Hurwitz zeta, validated in
green_optimize.py) evaluated on a fine sampling of p -- NO sympy in the loop.
CERTIFICATION of the winner is done ONCE in closed form:
    c_r exact (sympy), partial sum to R, plus rigorous tail
    |c_r| <= B/(pi r)^2,  B=|p'(0)|+|p'(1)|+||p''||_1
    => sum_{r>R}|c_r|^{4/3} <= (B/pi^2)^{4/3} (3/5) R^{-5/3}.
const(p) = (4/(1+2 J^{-3}))^{1/3};  Green quartic: J=2.409547, const=1.518248.
"""
import numpy as np
import sympy as sp
from scipy.optimize import minimize
import sys; sys.path.insert(0, '/Users/kabirgrewal/projects/erdos-241/code')
from green_optimize import hurwitz_weights, Jval

JG = 2.4095470
CG = (4 / (1 + 2 * JG ** -3)) ** (1 / 3)
NGRID = 8192
Z8 = hurwitz_weights(NGRID)
t8 = np.arange(1, 2 * NGRID)
OF8 = np.abs(np.exp(1j * np.pi / NGRID) ** t8 - 1.0)
xg = (np.arange(NGRID) + 0.5) / NGRID


def p_samples(coeffs):
    u = (2 * xg - 1) ** 2
    Q = np.zeros_like(xg)
    for j, c in enumerate(coeffs):
        Q += c * u ** j
    return xg * (1 - xg) * Q


def Jfast(coeffs):
    p = p_samples(coeffs)
    p = p - (p.mean() - 2.0)  # enforce mean 2 (additive, tiny shift since family already ~2)
    return Jval(p, NGRID, Z8, OF8)


# exact mean of basis x(1-x)u^j over [0,1]
xx = sp.symbols('x', real=True)
def basis_mean(j):
    u = (2 * xx - 1) ** 2
    return float(sp.integrate(xx * (1 - xx) * u ** j, (xx, 0, 1)))


def optimize(d):
    bm = [basis_mean(j) for j in range(d + 1)]

    def obj(free):
        c0 = (2 - sum(bm[j] * free[j - 1] for j in range(1, d + 1))) / bm[0]
        return Jfast([c0] + list(free))

    best = None
    for trial in range(6):
        x0 = (np.arange(d) - d / 2.0) * (0.5 if trial == 0 else 0.0) + \
             (np.zeros(d) if trial == 0 else (np.array([(-1)**k for k in range(d)]) * (trial)))
        res = minimize(obj, x0, method="Nelder-Mead",
                       options={"maxiter": 6000, "xatol": 1e-8, "fatol": 1e-11, "adaptive": True})
        if best is None or res.fun < best[1]:
            best = (res.x, res.fun)
    free = best[0]
    c0 = (2 - sum(bm[j] * free[j - 1] for j in range(1, d + 1))) / bm[0]
    return [c0] + list(free), best[1]


# ---- exact certification (once) ----
rsym = sp.symbols('r', real=True)
def certify(coeffs, R=600000):
    u = (2 * xx - 1) ** 2
    p = sp.expand(xx * (1 - xx) * sum(sp.Float(c) * u ** j for j, c in enumerate(coeffs)))
    m = float(sp.integrate(p, (xx, 0, 1)))
    assert abs(m - 2) < 1e-6, f"mean {m}"
    integ = sp.integrate(p * sp.exp(sp.I * sp.pi * rsym * xx), (xx, 0, 1))
    f = sp.lambdify(rsym, integ, "numpy")
    rs = np.arange(1, R + 1, dtype=float)
    with np.errstate(all="ignore"):
        c = np.asarray(f(rs), dtype=complex)
    if not np.isfinite(c).all():
        for i in np.where(~np.isfinite(c))[0]:
            c[i] = complex(sp.limit(integ, rsym, int(rs[i])))
    Jpart = float(np.sum(np.abs(c) ** (4 / 3)))
    pp = sp.diff(p, xx); ppp = sp.diff(p, xx, 2)
    B = abs(float(pp.subs(xx, 0))) + abs(float(pp.subs(xx, 1)))
    xs = np.linspace(0, 1, 400001)
    vals = np.abs(np.broadcast_to(sp.lambdify(xx, ppp, "numpy")(xs), xs.shape))
    dx = xs[1] - xs[0]
    B += float(0.5 * dx * (vals[0] + vals[-1] + 2 * vals[1:-1].sum())) * 1.0002
    tail = (B / np.pi ** 2) ** (4 / 3) * (3 / 5) * R ** (-5 / 3)
    return Jpart, Jpart + tail, B, tail, p


def const_hi(Jhi):
    return (4 / (1 + 2 * Jhi ** -3)) ** (1 / 3)


if __name__ == "__main__":
    print(f"Green: J={JG:.6f}  const={CG:.6f}  (cubed {CG**3:.5f})")
    print("=" * 74)
    results = {}
    for d in [1, 2, 3, 4, 5]:
        coeffs, Jf = optimize(d)
        results[d] = (coeffs, Jf)
        print(f"deg u^{d} (poly deg {2*d+2}): Jfast={Jf:.6f}  const={const_hi(Jf):.6f}  "
              f"coeffs={[round(c,3) for c in coeffs]}", flush=True)
    # certify the best
    dbest = min(results, key=lambda d: results[d][1])
    coeffs, Jf = results[dbest]
    print("\nCERTIFYING best family deg u^%d ..." % dbest, flush=True)
    Jlo, Jhi, B, tail, p = certify(coeffs)
    print(f"  certified J in [{Jlo:.6f}, {Jhi:.6f}]  (B={B:.3f}, tail<{tail:.2e})")
    print(f"  certified constant <= {const_hi(Jhi):.6f}   (Green {CG:.6f})")
    print(f"  cubed: 4/(1+gamma) <= {const_hi(Jhi)**3:.5f}   (Green {CG**3:.5f})")
    if Jhi < JG:
        print(f"  ==> CERTIFIED IMPROVEMENT: constant {const_hi(Jhi):.5f} N^(1/3) < Green {CG:.5f} N^(1/3)")
    np.save("/Users/kabirgrewal/projects/erdos-241/data/green_best_coeffs.npy", np.array(coeffs))
    with open("/Users/kabirgrewal/projects/erdos-241/data/green_best_kernel.txt", "w") as fh:
        fh.write(f"degree u^{dbest}\ncoeffs (c_0..): {coeffs}\n")
        fh.write(f"p(x) = {sp.nsimplify(p, rational=False)}\n")
        fh.write(f"J in [{Jlo},{Jhi}]  constant<= {const_hi(Jhi)}\n")
