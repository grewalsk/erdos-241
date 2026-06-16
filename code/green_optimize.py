"""Exact minimization of Green's kernel functional  J(p) = sum_{r>=1} |p~(pi r)|^{4/3}.

Green (Acta Arith. 100, 2001): for any real p in C^1[0,1] with int_0^1 p = 2,
every B_3 set A in {1..N} obeys |A| <= (4/(1+gamma(p)))^{1/3} N^{1/3}(1+o(1)),
where gamma(p) = 2 J(p)^{-3} and p~(lambda) = int_0^1 p(x) e^{i x lambda} dx.
His quartic p gives J = 2.4095, gamma = 1/6.9948, constant 3.4997 (rounded to 3.5).

=> The current world-record constant is governed by  inf_p J(p).
   constant(p) = (4/(1+2 J^{-3}))^{1/3}.  Smaller J  <=>  smaller (better) constant.
   J such that constant = 3.5 exactly: gamma = 1/7 => J = 14^{1/3} = 2.41014.

EXACT computation for piecewise-constant p on n equal bins (value p_j on
[j/n,(j+1)/n)).  With omega = e^{i pi / n}:
   p~(pi r) = (omega^r - 1)/(i pi r) * S(r),   S(r) = sum_j p_j omega^{r j}.
S has period 2n in r, and omega^r depends only on t = r mod 2n, so
   |p~(pi r)| = |omega^t - 1| * |S(t)| / (pi r),    t = r mod 2n.
Therefore, with no truncation error,
   J = (1/pi)^{4/3} sum_{t=1}^{2n-1} ( |omega^t - 1| |S(t)| )^{4/3} * Z(t),
   Z(t) = sum_{q>=0} (2 n q + t)^{-4/3} = (2n)^{-4/3} * zeta_Hurwitz(4/3, t/(2n)),
and the t=0 term vanishes (|omega^0-1| = 0).  S(t) for t=0..2n-1 is one FFT of
the zero-padded coefficient vector.  This is exact for step functions; the true
optimal p is smooth, so step functions only UNDER-resolve (give an upper bound on
inf J that improves as n grows) -- safe for the question "can we beat 3.5?".

We (i) validate against Green's closed-form numbers, (ii) convex-minimize J.
Convexity: J is a sum of (4/3)-powers of moduli of linear functionals of p, hence
convex; the constraint mean(p)=2 is affine; the global min is found.
"""
import numpy as np
from scipy.special import zeta as _zeta   # Hurwitz: zeta(s, a)
from scipy.optimize import minimize


def hurwitz_weights(n):
    """Z(t) = sum_{q>=0} (2 n q + t)^{-4/3}, t = 1..2n-1, exact via Hurwitz zeta."""
    N2 = 2 * n
    t = np.arange(1, N2)
    return (N2 ** (-4.0 / 3.0)) * _zeta(4.0 / 3.0, t / N2)


def Jval(p, n, Z=None, omega_fac=None):
    """Exact J(p) for step function p (length n)."""
    N2 = 2 * n
    if Z is None:
        Z = hurwitz_weights(n)
    if omega_fac is None:
        t = np.arange(1, N2)
        omega = np.exp(1j * np.pi / n)
        omega_fac = np.abs(omega ** t - 1.0)
    S = np.fft.fft(p, N2)            # S[t] = sum_j p_j omega^{t j}, omega = e^{i pi/n}?
    # np.fft uses e^{-2pi i t j / N2} = e^{-i pi t j / n} = omega^{-t j}.
    # We need sum p_j omega^{+t j}; that's conj(FFT) for real p. Use magnitude only:
    St = np.abs(S[1:N2])
    base = (omega_fac * St / np.pi) ** (4.0 / 3.0)
    return float(np.sum(base * Z))


def constant_of_J(J):
    """Returns (cuberoot_constant, gamma).  cuberoot_constant multiplies N^{1/3};
    its cube (4/(1+gamma)) is the '3.5' figure."""
    gamma = 2.0 * J ** (-3.0)
    return (4.0 / (1.0 + gamma)) ** (1.0 / 3.0), gamma


def J_and_grad(p, n, Z, omega_fac, omega_pow):
    """Exact J(p) and gradient dJ/dp for step function p (length n), with the
    mean-2 projection's null direction removed from the gradient.

    J = (1/pi)^{4/3} sum_{t=1}^{2n-1} (a_t |S_t|)^{4/3} Z_t,  a_t=|omega^t-1|,
    S_t = sum_j p_j omega^{t j} (omega=e^{i pi/n}); here we use np.fft which gives
    F_t = sum_j p_j e^{-2pi i t j/(2n)} = sum_j p_j omega^{-t j} = conj(S_t) (real p),
    so |S_t| = |F_t|.  d|S_t|^{4/3}/dp_j = (4/3)|S_t|^{-2/3} Re(conj(S_t) omega^{t j}).
    """
    N2 = 2 * n
    F = np.fft.fft(p, N2)                 # F_t = conj(S_t)
    St = np.abs(F[1:N2])
    coef = (omega_fac / np.pi) ** (4.0 / 3.0) * Z      # length 2n-1, multiplies |S_t|^{4/3}
    J = float(np.sum(coef * St ** (4.0 / 3.0)))
    # gradient: g_j = sum_t coef_t * (4/3) |S_t|^{-2/3} Re(conj(S_t) omega^{t j})
    # conj(S_t) = F_t.  Build weighted spectrum W_t = coef_t*(4/3)|S_t|^{-2/3} * F_t,
    # then g_j = Re sum_t W_t omega^{t j} = Re sum_t W_t e^{i pi t j/n}.
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(St > 1e-15, coef * (4.0 / 3.0) * St ** (-2.0 / 3.0), 0.0)
    W = np.zeros(N2, dtype=complex)
    W[1:N2] = scale * F[1:N2]
    # g_j = Re sum_{t=0}^{2n-1} W_t e^{+i pi t j/n} = Re( IFFT-like ).  e^{+i pi t j/n}
    # = e^{+2pi i t j/(2n)} = conj of np.fft.fft basis; so sum_t W_t e^{+...} = N2 * ifft(W).
    g_full = (N2 * np.fft.ifft(W)).real[:n]
    g = g_full - g_full.mean()           # remove component along the mean-constraint normal
    return J, g


def optimize_grad(n=512, Z=None, omega_fac=None, starts=None, maxiter=5000):
    from scipy.optimize import minimize as _min
    N2 = 2 * n
    if Z is None:
        Z = hurwitz_weights(n)
    t = np.arange(1, N2)
    if omega_fac is None:
        omega_fac = np.abs(np.exp(1j * np.pi / n) ** t - 1.0)
    x = (np.arange(n) + 0.5) / n

    def proj(p):
        return p - (p.mean() - 2.0)

    def fg(p):
        pp = proj(p)
        J, g = J_and_grad(pp, n, Z, omega_fac, None)
        return J, g

    if starts is None:
        starts = [
            ("const2", np.full(n, 2.0)),
            ("green", 2.5 - 40 * (x - 0.5) ** 4),
            ("deg6", 2.0 + 6 * ((x - 0.5) ** 2 * 4 - 1.0 / 3) - 200 * (x - 0.5) ** 6),
            ("cos1", 2.0 + 1.0 * np.cos(2 * np.pi * x)),
        ]
    best = None
    for name, p0 in starts:
        p0 = proj(np.asarray(p0, float))
        res = _min(fg, p0, jac=True, method="L-BFGS-B",
                   options={"maxiter": maxiter, "ftol": 1e-15, "gtol": 1e-12})
        if best is None or res.fun < best[1]:
            best = (name, res.fun, proj(res.x))
    return best


# ----- exact J for a function given by samples of a SMOOTH p (closed-form-ish):
def J_smooth_quartic(R=2_000_000):
    """Green's quartic via his closed forms (26)-(29): independent check."""
    r = np.arange(1, R + 1, dtype=float)
    even = r[1::2]; odd = r[0::2]
    S1 = np.sum(np.abs(1/even**2 - 24/(np.pi**2*even**4))**(4/3))
    S2 = np.sum(np.abs(1/odd**3 - 8/(np.pi**2*odd**5))**(4/3))
    return (40/np.pi**2)**(4/3) * (S1 + (6/np.pi)**(4/3) * S2)


def optimize(n=256, restarts=None):
    Z = hurwitz_weights(n)
    N2 = 2 * n
    t = np.arange(1, N2)
    omega_fac = np.abs(np.exp(1j*np.pi/n) ** t - 1.0)
    x = (np.arange(n) + 0.5) / n

    def proj(p):
        return p - (p.mean() - 2.0)

    def obj(p):
        return Jval(proj(p), n, Z, omega_fac)

    starts = restarts or [
        ("const2", np.full(n, 2.0)),
        ("green", 2.5 - 40*(x-0.5)**4),
        ("cos", 2.0 + 1.5*np.cos(2*np.pi*x)),
        ("U-ramp", np.where(x < 0.5, 4.0*(1-x), 4.0*x*0+ 0.0) ),
    ]
    best = None
    for name, p0 in starts:
        p0 = proj(np.asarray(p0, float))
        res = minimize(obj, p0, method="Nelder-Mead",
                       options={"maxiter": 40000, "xatol": 1e-9, "fatol": 1e-12,
                                "adaptive": True})
        J = res.fun
        if best is None or J < best[1]:
            best = (name, J, proj(res.x))
    return best


if __name__ == "__main__":
    print("=" * 72)
    print("Validation: exact step-function J vs Green's closed-form quartic J")
    print("=" * 72)
    Jq = J_smooth_quartic()
    cq, gq = constant_of_J(Jq)
    print(f"Green quartic (closed form):  J = {Jq:.6f}  gamma = 1/{1/gq:.4f}  const = {cq:.6f}")
    for n in [64, 256, 1024, 4096]:
        x = (np.arange(n)+0.5)/n
        pg = 2.5 - 40*(x-0.5)**4
        pg = pg - pg.mean() + 2.0
        pc = np.full(n, 2.0)
        Jg = Jval(pg, n); Jc = Jval(pc, n)
        cg,_ = constant_of_J(Jg); cc,_ = constant_of_J(Jc)
        print(f"  n={n:5d}: step-J(green quartic)={Jg:.6f} (const {cg:.5f}) | "
              f"step-J(const 2)={Jc:.6f} (const {cc:.5f})")

    print()
    print("=" * 72)
    print("Gradient-based convex minimization of J (true inf over step functions)")
    print("=" * 72)
    Jgreen_cf = Jq
    cgreen, ggreen = constant_of_J(Jgreen_cf)
    print(f"  Green's reference: J = {Jgreen_cf:.6f}, gamma = 1/{1/ggreen:.4f}, "
          f"cubed-constant 4/(1+g) = {4/(1+ggreen):.5f}")
    print(f"  To BEAT Green we need J strictly below {Jgreen_cf:.6f}.")
    for n in [128, 256, 512, 1024]:
        name, J, p = optimize_grad(n=n)
        c, g = constant_of_J(J)
        cubed = 4 / (1 + g)
        verdict = ("BEATS Green" if J < Jgreen_cf - 1e-4
                   else "matches Green" if abs(J - Jgreen_cf) <= 1e-3
                   else "above Green")
        print(f"  n={n:5d}: min J = {J:.6f}  gamma = 1/{1/g:.4f}  cubed = {cubed:.5f}  "
              f"[{verdict}] (best start: {name})")
        np.save(f"/Users/kabirgrewal/projects/erdos-241/data/popt_n{n}.npy", p)
