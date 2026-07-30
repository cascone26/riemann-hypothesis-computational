#!/usr/bin/env python3
"""Weyl: Dirac's pre-registered decisive test. Exact ODE for u=x*xi (§464 build, lam=sqrt2):
    u'' + bracket(x,alpha) u = 0,
    bracket(x,alpha) = 8pi^2 + (16pi^2-alpha)/(x^2-2) + 2/(x^2-2)^2      [EXACT -- derived by hand,
    no series truncation: bracket = g - f'/2 - f^2/4 for f=2x/(x^2-2), g=(8pi^2 x^2-alpha)/(x^2-2)].
Naive build matches u(x) at finite X against RAW {sin(w x), cos(w x)}, w=2pi*sqrt2 -- that basis solves
u''+8pi^2 u=0 exactly, i.e. it drops the (16pi^2-alpha)/(x^2-2)+2/(x^2-2)^2 correction entirely. This is
Dirac's "bare-sinusoid match". The decisive replacement: use the exact local WKB phase/amplitude of
bracket(x,alpha) itself (k(x)=sqrt(bracket), psi(x)=int k dx, amplitude ~ k^{-1/2}) as the reference basis,
carrying the true, alpha-DEPENDENT phase (this is Dirac's "u''+(w^2-alpha/x^2)u=0" comparison equation,
done here to full exactness in x rather than the leading 1/x^2 truncation, since the exact bracket() is
cheap and already in hand).
Test: does replacing bare-sin/cos with the WKB-corrected basis collapse the ~3% floor, and does the
residual now shrink monotonically with X (finite-window artifact) rather than sit pinned (intrinsic)?
"""
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

Z = np.load('zeros500.npy')
FP = 8*np.pi**2
w  = 2*np.pi*np.sqrt(2)
sp = np.sqrt(2.0)

def bracket(x, alpha):
    d = x*x - 2.0
    return FP + (2*FP - alpha)/d + 2.0/d**2

def rhs(x, y, alpha):
    xi, xip = y
    xipp = -(2*x*xip + (FP*x*x - alpha)*xi) / (x*x - 2)
    return [xip, xipp]

def solve_tail(gamma, X, npts=1600):
    alpha = -(gamma/2.0)**2
    dxi0 = -(2*FP - alpha) / (2*sp)
    xs = np.linspace(sp+1e-6, X, npts)
    s = solve_ivp(rhs, [sp+1e-6, X], [1.0, dxi0], args=(alpha,), method='DOP853',
                   rtol=1e-12, atol=1e-14, t_eval=xs)
    u = xs * s.y[0]
    return xs, u, alpha

def B_naive(gamma, X, tailfrac=0.15):
    xs, u, alpha = solve_tail(gamma, X)
    n = int(len(xs)*(1-tailfrac))
    xt, ut = xs[n:], u[n:]
    M = np.column_stack([np.sin(w*xt), np.cos(w*xt)])
    (A, B), *_ = np.linalg.lstsq(M, ut, rcond=None)
    return B

def B_wkb(gamma, X, tailfrac=0.15):
    # WKB phase must be integrated only over the tail window (far from the x=sp pole in bracket()
    # where naive quadrature garbles the phase) -- reference psi=0 at the tail window's own start.
    xs, u, alpha = solve_tail(gamma, X)
    n = int(len(xs)*(1-tailfrac))
    xt, ut = xs[n:], u[n:]
    k = np.sqrt(bracket(xt, alpha))
    psit = cumulative_trapezoid(k, xt, initial=0.0)         # psi=0 at xt[0], local to the tail only
    ampt = k**-0.5
    M = np.column_stack([ampt*np.sin(psit), ampt*np.cos(psit)])
    (A, B), *_ = np.linalg.lstsq(M, ut, rcond=None)
    return B

def roots_from(Bfun, X, gs):
    Bs = np.array([Bfun(g, X) for g in gs])
    roots = []
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1] < 0:
            a, b = gs[i], gs[i+1]
            for _ in range(30):
                m = (a+b)/2
                if Bfun(a, X)*Bfun(m, X) < 0: b = m
                else: a = m
            roots.append((a+b)/2)
    return np.array(roots)

gs = np.linspace(12, 32, 100)
print(f"{'X':>5} | {'method':8} | {'#roots':>6} | mean%err(first 6) | first 6 %err")
for X in [12.0, 20.0, 35.0, 60.0]:
    zin = Z[(Z>=12)&(Z<=32)]
    for name, fn in [('naive', B_naive), ('wkb', B_wkb)]:
        r = roots_from(fn, X, gs)
        nn = min(len(r), len(zin), 6)
        if nn == 0:
            print(f"{X:5.0f} | {name:8} | {len(r):6d} | (no roots in window)")
            continue
        errs = 100*(r[:nn]-zin[:nn])/zin[:nn]
        mean_abs = np.mean(np.abs(errs))
        print(f"{X:5.0f} | {name:8} | {len(r):6d} | {mean_abs:8.3f}%          | " +
              ", ".join(f"{e:+.2f}%" for e in errs))
print("""
reading: if 'wkb' mean%err << 'naive' at fixed X, and/or 'wkb' error keeps shrinking with X while 'naive'
plateaus near 3%, the floor was a bare-sinusoid phase-matching artifact (Dirac's hypothesis CONFIRMED,
Report-tier numpy/scipy -- not yet Arb-certified). If 'wkb' still floors near the same ~3%, the floor
survives the exact local-WKB correction and is NOT explained by the phase tail alone -- points at either
the x=sqrt(2) boundary condition (theta_reg) or something intrinsic to D.
""")
