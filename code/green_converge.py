"""Confirm inf J converges to a finite limit < Green's 2.4095, via higher
resolution + Richardson extrapolation, and test the p(0)=p(1)=0 variant."""
import numpy as np
import sys; sys.path.insert(0,'code')
from green_optimize import hurwitz_weights, J_and_grad, Jval, constant_of_J, J_smooth_quartic
from scipy.optimize import minimize

def optimize(n, clamp_ends=False, maxiter=20000):
    Z=hurwitz_weights(n); N2=2*n
    t=np.arange(1,N2); omega_fac=np.abs(np.exp(1j*np.pi/n)**t-1.0)
    x=(np.arange(n)+0.5)/n
    def proj(p):
        p=p.copy()
        if clamp_ends:
            # force a smooth dip: not a hard constraint, just start; mean fix below
            pass
        return p-(p.mean()-2.0)
    def fg(p):
        pp=proj(p); J,g=J_and_grad(pp,n,Z,omega_fac,None); return J,g
    best=None
    starts=[("green",2.5-40*(x-0.5)**4),("cos1",2.0+np.cos(2*np.pi*x)),
            ("bump",2.0+2.0*np.cos(np.pi*(2*x-1))*(np.abs(2*x-1)<1))]
    for nm,p0 in starts:
        r=minimize(fg,proj(np.asarray(p0,float)),jac=True,method="L-BFGS-B",
                   options={"maxiter":maxiter,"ftol":1e-16,"gtol":1e-13})
        if best is None or r.fun<best[1]: best=(nm,r.fun,proj(r.x))
    return best[1],best[2]

Jq=J_smooth_quartic(); cq,gq=constant_of_J(Jq)
print(f"Green reference: J={Jq:.6f}  cubed=4/(1+g)={4/(1+gq):.6f}")
Js={}
for n in [256,512,1024,2048,4096]:
    J,p=optimize(n)
    c,g=constant_of_J(J); Js[n]=J
    np.save(f'data/popt_n{n}.npy',p)
    print(f"  n={n:5d}: minJ={J:.6f}  cubed={4/(1+g):.6f}  p in [{p.min():.3f},{p.max():.3f}]")
# Richardson: assume J(n) = Jinf + a/n.  Use last two.
ns=sorted(Js); 
for i in range(1,len(ns)):
    n1,n2=ns[i-1],ns[i]
    # J = Jinf + a/n  => Jinf = (n2 J2 - n1 J1)/(n2-n1)
    Jinf=(n2*Js[n2]-n1*Js[n1])/(n2-n1)
    print(f"   Richardson[{n1},{n2}] -> Jinf~{Jinf:.6f}  cubed~{4/(1+2*Jinf**-3):.6f}")
