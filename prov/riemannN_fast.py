#!/usr/bin/env python3
"""Fast independent shooting for CCM lam=sqrt2 operator eigenvalues, multiprocessed.
See riemannN_delta_Sscan.py for the derivation. Exact transform u=sqrt(x^2-2)*xi."""
import numpy as np, sys, time
from scipy.integrate import solve_ivp, cumulative_trapezoid
from multiprocessing import Pool

FP = 8*np.pi**2
sp = np.sqrt(2.0)
X = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0
TR = sys.argv[2] if len(sys.argv) > 2 else 'exact'

def Q(x, a):
    d = x*x - 2.0
    return FP + (2*FP - a)/d + 2.0/d**2

def rhs(x, y, a):
    return [y[1], -(2*x*y[1] + (FP*x*x - a)*y[0])/(x*x - 2.0)]

def Bcos(gamma):
    a = -(gamma/2.0)**2
    dxi0 = -(2*FP - a)/(2*sp)
    x0 = sp + 1e-6
    npts = 1600
    xs = np.linspace(x0, X, npts)
    s = solve_ivp(rhs, [x0, X], [1.0, dxi0], args=(a,), method='DOP853',
                  rtol=1e-10, atol=1e-12, t_eval=xs, dense_output=False)
    xi = s.y[0]
    u = (np.sqrt(xs*xs-2.0) if TR=='exact' else xs)*xi
    n = int(npts*0.85)
    xt, ut = xs[n:], u[n:]
    k = np.sqrt(Q(xt, a))
    psit = cumulative_trapezoid(k, xt, initial=0.0)
    amp = k**-0.5
    M = np.column_stack([amp*np.sin(psit), amp*np.cos(psit)])
    (A, B), *_ = np.linalg.lstsq(M, ut, rcond=None)
    return B

def bisect(a, b, fa, fb, it=45):
    for _ in range(it):
        m = 0.5*(a+b); fm = Bcos(m)
        if fa*fm <= 0: b, fb = m, fm
        else: a, fa = m, fm
    return 0.5*(a+b)

if __name__ == '__main__':
    t0 = time.time()
    gs = np.linspace(13.5, 96.0, 330)
    with Pool(8) as p:
        Bs = np.array(p.map(Bcos, gs))
    roots = []
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1] < 0:
            roots.append(bisect(gs[i], gs[i+1], Bs[i], Bs[i+1]))
    roots = np.array(roots)
    np.save(f'/tmp/ops_X{int(X)}_{TR}.npy', roots)
    print(f"# X={X} tr={TR} #eig={len(roots)}  ({time.time()-t0:.1f}s)")
    for i, r in enumerate(roots):
        print(f"{i+1:3d} {r:12.6f}")
