#!/usr/bin/env python3
"""Task 1 (roundtable §427 handoff): confirm the BARE build's S(T) correlation r stays ~0 across
X in {60,80,100,130} -> the 'wrong-shape' is not an X-artifact. Fast settings (tail-corrected phase,
rtol 1e-8); grades by pearsonr(F_op, S_true) which is invariant to the WKB tail offset.
Second instrument for §427 (roundtable vn_FvsS.py got r~0.095 at X=80)."""
import os, numpy as np
from scipy.integrate import solve_ivp
from scipy.stats import pearsonr
_ZP = 'zeros500.npy' if os.path.exists('zeros500.npy') else 'operators/zeros500.npy'
Z = np.load(_ZP)
FP = 8*np.pi**2; w = 2*np.pi*np.sqrt(2); sp = np.sqrt(2.0); eps = 1e-6
GMAX = 95.0

def rhs(x, y, alpha):
    xi, xip = y
    return [xip, -(2*x*xip + (FP*x*x - alpha)*xi)/(x*x - 2)]

def Bval(gamma, X):
    alpha = -(gamma/2.0)**2
    dxi0 = -(2*FP - alpha)/(2*sp)
    xs = np.linspace(X-6, X, 120)
    s = solve_ivp(rhs, [sp+eps, X], [1.0, dxi0], args=(alpha,), method='DOP853',
                  rtol=1e-8, atol=1e-10, t_eval=xs)
    u = xs*s.y[0]
    ph = w*xs - (gamma/2.0)**2/(2*w*xs)
    M = np.column_stack([np.sin(ph), np.cos(ph)])
    (A, B), *_ = np.linalg.lstsq(M, u, rcond=None)
    return B

def bare_roots(X):
    gs = np.arange(3.0, GMAX, 0.10)
    Bs = np.array([Bval(g, X) for g in gs])
    r = []
    for i in range(len(gs)-1):
        if Bs[i]*Bs[i+1] < 0:
            a, b = gs[i], gs[i+1]
            for _ in range(22):
                m = (a+b)/2
                if Bval(a, X)*Bval(m, X) < 0: b = m
                else: a = m
            r.append((a+b)/2)
    return np.array(r)

def Nsmooth(T): return (T/(2*np.pi))*np.log(T/(2*np.pi*np.e)) + 7.0/8.0
grid = np.linspace(14.0, GMAX, 300)
Ssm = Nsmooth(grid)
S_true = np.searchsorted(Z, grid, side='right') - Ssm
nz = int(np.sum(Z < GMAX)); ms = 2*np.pi/np.log(GMAX/(2*np.pi))
print(f"Task 1 bare-build S(T) correlation vs X.  n_true(<{GMAX})={nz}, spacing~{ms:.3f}")
print(f"{'X':>6} {'nroot':>6} {'r':>8} {'RMS':>7} {'RMS/spacing':>12}", flush=True)
for X in [60.0, 80.0, 100.0, 130.0]:
    roots = bare_roots(X)
    Nop = np.searchsorted(np.sort(roots), grid, side='right')
    Fop = Nop - Ssm
    r, _ = pearsonr(Fop, S_true)
    rms = np.sqrt(np.mean((Fop - S_true)**2))
    print(f"{X:6.0f} {len(roots):6d} {r:8.4f} {rms:7.3f} {rms/ms:12.3f}", flush=True)
