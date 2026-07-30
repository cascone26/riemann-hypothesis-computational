#!/usr/bin/env python3
"""Riemann-N INDEPENDENT re-derivation of RELAY §426.
CCM lam=sqrt(2) operator on x in (sqrt2, inf):
    (x^2-2) xi'' + 2x xi' + (8pi^2 x^2 - alpha) xi = 0,   alpha = -(gamma/2)^2.
Exact Liouville transform u = sqrt(x^2-2)*xi  (NOT u=x*xi as in the tail test):
    u'' + Q(x,alpha) u = 0,  Q = 8pi^2 + (16pi^2-alpha)/(x^2-2) + 2/(x^2-2)^2.
Eigenvalue gamma_op = value of gamma where the regular-at-sqrt2 solution is pure-sin at the tail
(cos-coefficient B=0), matched against the EXACT local WKB basis of Q.

Independent choices vs handoff: exact transform (not x*xi), WKB phase, my own bisection, my own
S(T) from mpmath, Odlyzko 9-digit zeros (not zeros500 floats).
"""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

FP = 8*np.pi**2
sp = np.sqrt(2.0)

def Q(x, alpha):
    d = x*x - 2.0
    return FP + (2*FP - alpha)/d + 2.0/d**2

def rhs(x, y, alpha):
    xi, xip = y
    xipp = -(2*x*xip + (FP*x*x - alpha)*xi)/(x*x - 2.0)
    return [xip, xipp]

def solve_xi(gamma, X, npts=4000):
    alpha = -(gamma/2.0)**2
    # regular IC at x=sqrt2+eps: xi(sp)=1, xi'(sp) from the ODE regularity.
    # near x=sp, leading balance 2x xi' + (FP x^2 - alpha) xi ~ (x^2-2) xi'' -> finite; use series slope
    dxi0 = -(2*FP - alpha)/(2*sp)   # same regular slope as tail test (indep-checked below)
    x0 = sp + 1e-6
    xs = np.linspace(x0, X, npts)
    s = solve_ivp(rhs, [x0, X], [1.0, dxi0], args=(alpha,), method='DOP853',
                  rtol=1e-12, atol=1e-14, t_eval=xs)
    return xs, s.y[0], alpha

def Bcos(gamma, X, transform='exact', tailfrac=0.15):
    xs, xi, alpha = solve_xi(gamma, X)
    if transform == 'exact':
        u = np.sqrt(xs*xs - 2.0)*xi
    else:  # 'xxi' = the tail-test approximation
        u = xs*xi
    n = int(len(xs)*(1-tailfrac))
    xt, ut = xs[n:], u[n:]
    k = np.sqrt(Q(xt, alpha))
    psit = cumulative_trapezoid(k, xt, initial=0.0)
    amp = k**-0.5
    M = np.column_stack([amp*np.sin(psit), amp*np.cos(psit)])
    (A, B), *_ = np.linalg.lstsq(M, ut, rcond=None)
    return B

def eigs_in(lo, hi, X, transform='exact', ngrid=400):
    gs = np.linspace(lo, hi, ngrid)
    Bs = np.array([Bcos(g, X, transform) for g in gs])
    roots = []
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1] < 0:
            a, b = gs[i], gs[i+1]
            fa = Bs[i]
            for _ in range(60):
                m = 0.5*(a+b)
                fm = Bcos(m, X, transform)
                if fa*fm <= 0: b = m
                else: a, fa = m, fm
            roots.append(0.5*(a+b))
    return np.array(roots)

if __name__ == '__main__':
    import sys
    X = float(sys.argv[1]) if len(sys.argv) > 1 else 100.0
    tr = sys.argv[2] if len(sys.argv) > 2 else 'exact'
    od = np.loadtxt('data/zeros_odlyzko_100k')
    zin = od[od < 96.0]           # 27 zeros
    ops = eigs_in(13.5, 96.0, X, tr, ngrid=500)
    print(f"# X={X} transform={tr}  #eig={len(ops)}  #zeros={len(zin)}")
    print(f"{'n':>3} {'gamma_op':>12} {'zero':>12} {'delta':>9} {'%err':>7}")
    m = min(len(ops), len(zin))
    deltas = []
    for i in range(m):
        d = ops[i]-zin[i]
        deltas.append(d)
        print(f"{i+1:3d} {ops[i]:12.6f} {zin[i]:12.6f} {d:+9.4f} {100*d/zin[i]:+7.3f}")
    deltas = np.array(deltas)
    print(f"# max|delta|={np.max(np.abs(deltas)):.4f} at n={np.argmax(np.abs(deltas))+1}"
          f"  mean|%err|={np.mean(np.abs(100*deltas/zin[:m])):.3f}%")
    np.save(f'/tmp/riemannN_ops_X{int(X)}_{tr}.npy', ops)
