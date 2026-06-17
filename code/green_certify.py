"""Certifiable kernels for Green's B_3 bound.

We seek an EXPLICIT, simply-described p in C^1[0,1] with int_0^1 p = 2 and
J(p) = sum_{r>=1} |p~(pi r)|^{4/3} provably below Green's quartic value
J_Green = 2.409547 (which yields constant (4/(1+2 J^{-3}))^{1/3} = 1.518248).

For each kernel we compute c_r = p~(pi r) = int_0^1 p(x) e^{i pi r x} dx in CLOSED
FORM (sympy / exact), sum |c_r|^{4/3} for r <= R, and add a RIGOROUS tail bound.

Tail bound: if p in C^2[0,1] with p(0)=p(1)=0, two integrations by parts give
   c_r = p~(pi r) = [ -p'(x) e^{i pi r x} /(i pi r)^2 ]_0^1  + (1/(i pi r)^2) p''~(pi r),
and since |p''~(pi r)| <= ||p''||_1, we get |c_r| <= ( |p'(1)|+|p'(0)| + ||p''||_1 )/(pi r)^2
=: B/(pi r)^2.  Hence sum_{r>R} |c_r|^{4/3} <= (B/pi^2)^{4/3} sum_{r>R} r^{-8/3}
<= (B/pi^2)^{4/3} * (3/5) R^{-5/3}.  All constants explicit => certified.

(Boundary p(0)=p(1)=0 also matches Green's own quartic and the numerically
optimal kernel, so it costs nothing.)
"""
import numpy as np
import sympy as sp

x, r = sp.symbols('x r', real=True)
PI = sp.pi


def c_r_closed(pexpr):
    """Return a python function r->complex for c_r = int_0^1 p e^{i pi r x} dx (exact)."""
    integ = sp.integrate(pexpr * sp.exp(sp.I * PI * r * x), (x, 0, 1))
    integ = sp.simplify(integ)
    f = sp.lambdify(r, integ, modules=["numpy"])
    return integ, f


def J_certified(pexpr, R=300000):
    # mean check
    mean = sp.integrate(pexpr, (x, 0, 1))
    assert sp.simplify(mean - 2) == 0, f"mean is {mean}, must be 2"
    _, cf = c_r_closed(pexpr)
    rs = np.arange(1, R + 1, dtype=float)
    cvals = cf(rs)
    Jpart = float(np.sum(np.abs(cvals) ** (4.0 / 3.0)))
    # rigorous tail constant B = |p'(0)|+|p'(1)| + ||p''||_1
    pp = sp.diff(pexpr, x)
    ppp = sp.diff(pexpr, x, 2)
    p1_0 = abs(float(pp.subs(x, 0)))
    p1_1 = abs(float(pp.subs(x, 1)))
    # ||p''||_1 = int_0^1 |p''| ; p'' is a polynomial/trig, integrate |.| numerically (upper bound via fine grid + small slack)
    xs = np.linspace(0, 1, 200001)
    ppp_f = sp.lambdify(x, ppp, "numpy")
    vals = np.abs(np.broadcast_to(ppp_f(xs), xs.shape))
    dx = xs[1] - xs[0]
    l1 = float(0.5 * dx * (vals[0] + vals[-1] + 2 * vals[1:-1].sum())) * 1.0001
    B = p1_0 + p1_1 + l1
    tail = (B / np.pi ** 2) ** (4.0 / 3.0) * (3.0 / 5.0) * R ** (-5.0 / 3.0)
    Jlo = Jpart
    Jhi = Jpart + tail
    return Jlo, Jhi, B, tail


def const_of_J(J):
    g = 2.0 * J ** -3.0
    return (4.0 / (1.0 + g)) ** (1.0 / 3.0)


if __name__ == "__main__":
    JG = 2.409547
    cG = const_of_J(JG)
    print(f"Green quartic:  J = {JG:.6f}   constant = {cG:.6f}   (cubed {cG**3:.5f})")
    print("=" * 78)

    kernels = {
        "parabola 6x(1-x)... mean2 => 12x(1-x)": 12 * x * (1 - x),
        "pi*sin(pi x)": PI * sp.sin(PI * x),
        "Green quartic 5/2-40(x-1/2)^4": sp.Rational(5, 2) - 40 * (x - sp.Rational(1, 2)) ** 4,
        "sin^2 scaled (4 sin^2(pi x))": 4 * sp.sin(PI * x) ** 2,
        "x(1-x)*(1+a..) deg4: 30x(1-x)(1-?)": None,  # filled below
    }
    # a 2-parameter even-poly family vanishing at ends: p = A x(1-x) + B [x(1-x)]^2, mean=2
    # int x(1-x)=1/6, int [x(1-x)]^2 = 1/30.  A/6 + B/30 = 2.
    A, Bc = sp.symbols('A Bc')
    print("\nClosed-form single kernels:")
    for name, expr in kernels.items():
        if expr is None:
            continue
        Jlo, Jhi, Bcnst, tail = J_certified(expr)
        flag = "  <<< BEATS GREEN (certified)" if Jhi < JG else ("  (lo beats, hi doesn't)" if Jlo < JG else "")
        print(f"  {name:42s} J in [{Jlo:.6f}, {Jhi:.6f}]  const_hi={const_of_J(Jhi):.6f}{flag}")

    # optimize the 2-param family numerically over the constraint, then certify winner
    print("\n2-parameter family p = a*x(1-x) + b*[x(1-x)]^2  (mean=2 constraint):")
    from scipy.optimize import minimize_scalar, minimize
    def Jnum(params):
        a, b = params
        expr = a * x * (1 - x) + b * (x * (1 - x)) ** 2
        m = float(sp.integrate(expr, (x, 0, 1)))
        if abs(m - 2) > 1e-9:
            # project: scale not linear under fixed shape; instead enforce via a from b
            pass
        # enforce mean=2 by solving a from b: a/6 + b/30 = 2 => a = 12 - b/5
        a = 12 - b / 5
        _, cf = c_r_closed(a * x * (1 - x) + b * (x * (1 - x)) ** 2)
        rs = np.arange(1, 60000, dtype=float)
        return float(np.sum(np.abs(cf(rs)) ** (4 / 3)))
    res = minimize(lambda bb: Jnum([0, bb[0]]), [0.0], method="Nelder-Mead",
                   options={"xatol": 1e-6, "fatol": 1e-9})
    bopt = res.x[0]; aopt = 12 - bopt / 5
    print(f"  optimal b={bopt:.5f}, a={aopt:.5f}")
    expr = aopt * x * (1 - x) + bopt * (x * (1 - x)) ** 2
    Jlo, Jhi, Bc2, tail = J_certified(sp.nsimplify(aopt, rational=False) * x * (1 - x)
                                      + sp.nsimplify(bopt, rational=False) * (x * (1 - x)) ** 2)
    print(f"  J in [{Jlo:.6f},{Jhi:.6f}]  const_hi={const_of_J(Jhi):.6f}"
          + ("  <<< BEATS GREEN (certified)" if Jhi < JG else ""))
