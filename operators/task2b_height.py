#!/usr/bin/env python3
"""Height-robustness of the §466 NULL: does bare F_0 stay uncorrelated with S(T) up to gamma<=150?
And does |F_0| stay O(1) (bounded) while |S(T)| grows? (empirical support for the Task-3 O(1)-vs-unbounded
argument). Bare build, tail-corrected phase, X=100."""
import os, numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr
_ZP = 'zeros500.npy' if os.path.exists('zeros500.npy') else 'operators/zeros500.npy'
Z = np.load(_ZP)
FP = 8*np.pi**2; w = 2*np.pi*np.sqrt(2); sp = np.sqrt(2.0); eps = 1e-6
X = 100.0

def rhs(x, y, alpha):
    xi, xip = y
    return [xip, -(2*x*xip + (FP*x*x - alpha)*xi)/(x*x - 2)]
def Bval(gamma):
    alpha = -(gamma/2.0)**2; dxi0 = -(2*FP - alpha)/(2*sp)
    xs = np.linspace(X-6, X, 120)
    s = solve_ivp(rhs, [sp+eps, X], [1.0, dxi0], args=(alpha,), method='DOP853',
                  rtol=1e-8, atol=1e-10, t_eval=xs)
    u = xs*s.y[0]; ph = w*xs - (gamma/2.0)**2/(2*w*xs)
    M = np.column_stack([np.sin(ph), np.cos(ph)]); (A, B), *_ = np.linalg.lstsq(M, u, rcond=None); return B
def bare_roots(gmax):
    gs = np.arange(3.0, gmax, 0.10); Bs = np.array([Bval(g) for g in gs]); r=[]
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1] < 0:
            a,b=gs[i],gs[i+1]
            for _ in range(22):
                m=(a+b)/2
                if Bval(a)*Bval(m)<0: b=m
                else: a=m
            r.append((a+b)/2)
    return np.sort(np.array(r))
def Nsmooth(T): return (T/(2*np.pi))*np.log(T/(2*np.pi*np.e)) + 7.0/8.0

for GMAX in [95.0, 120.0, 150.0]:
    grid = np.linspace(14.0, GMAX, int((GMAX-14)*4))
    Ssm = Nsmooth(grid); S_true = np.searchsorted(Z, grid, side='right') - Ssm
    roots = bare_roots(GMAX)
    Fop = np.searchsorted(roots, grid, side='right') - Ssm
    r,_ = pearsonr(Fop, S_true); rms = np.sqrt(np.mean((Fop-S_true)**2))
    ms = 2*np.pi/np.log(GMAX/(2*np.pi))
    print(f"GMAX={GMAX:6.0f}: nroot={len(roots):3d} nzeros={int(np.sum(Z<GMAX)):3d}  r={r:7.4f}  "
          f"RMS(F-S)={rms:.3f}  max|F_0|={np.max(np.abs(Fop)):.2f}  max|S|={np.max(np.abs(S_true)):.2f}",
          flush=True)
