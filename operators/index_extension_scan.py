#!/usr/bin/env python3
"""Find the correct self-adjoint extension of the CCM prolate operator IN THE INDEX (matrix) picture,
by scanning a boundary parameter — the index-space analog of the §466 U(2) scan. The pentadiagonal
W = -J^2 + 4pi^2(4N+1) - 1/4 (lambda=sqrt2) decouples into even/odd index sublattices (each a tridiagonal
3-term recurrence). Hard truncation (Friedrichs) gives the WRONG gamma=17.1 (§469c); deficiency (2,2) means
each sublattice has a 1-parameter U(1) family of extensions set by a boundary term at the large-n (singular
x=inf) end. Scan that boundary parameter t and look for the value(s) where the low gamma=2*sqrt|neg eig|
hit the TRUE zeros. HARD GATE: must reproduce base zeros 14.13/21.02/25.01/30.42/32.94 before trusting any
semilocal run.
"""
import numpy as np

def bcoef_base(n):   # b_n = sqrt((n+1/2)(n+1)), the arithmetic-free base scaling Jacobi off-diagonal
    return np.sqrt((n+0.5)*(n+1.0))

def build_W_sublattice(bfun, K, parity, t_boundary=0.0, lam=np.sqrt(2.0)):
    """Tridiagonal even/odd-index block of W with a boundary self-energy t at the last (large-n) site."""
    ns = np.arange(parity, 2*K, 2)        # indices of this sublattice
    m = len(ns)
    b = np.array([bfun(n) for n in range(2*K+3)])
    d = 4*np.pi**2*(4*ns+1) - 0.25
    diag = d - (b[ns]**2 + b[ns+1]**2)                 # W_{n,n} = d_n - (b_n^2 + b_{n+1}^2)
    off  = -b[ns[:-1]+1]*b[ns[:-1]+2]                  # W_{n,n+2} = -b_{n+1} b_{n+2}
    T = np.diag(diag) + np.diag(off,1) + np.diag(off,-1)
    T[-1,-1] += t_boundary                             # boundary extension parameter at the singular end
    return T

def low_gammas(bfun, K, t_e, t_o, nlow=6):
    """Return the smallest-|neg| gammas (=2*sqrt|neg eig|) pooling even+odd sublattices."""
    gs = []
    for parity, t in [(0,t_e),(1,t_o)]:
        ev = np.sort(np.linalg.eigvalsh(build_W_sublattice(bfun, K, parity, t)))
        neg = ev[ev<0]
        if len(neg):
            gs += list(2*np.sqrt(-neg[-nlow:]))        # least-negative -> smallest gamma
    return np.sort(np.array(gs))

TRUE = np.array([14.1347,21.0220,25.0109,30.4249,32.9351,37.5862])
K = 120

# Baseline (Friedrichs, t=0) for reference
g0 = low_gammas(bcoef_base, K, 0.0, 0.0)
print(f"Friedrichs (t=0) low gammas: {np.round(g0[:6],3)}   (should be the wrong ~17.1/28/... )")
print(f"TRUE zeros:                  {np.round(TRUE,3)}\n")

# Scan boundary parameters. Grade each (t_e,t_o) by how well the low gammas match TRUE zeros.
def match_score(g, K=5):
    # for each of the first K true zeros, distance to nearest operator gamma
    if len(g)==0: return 1e9
    return np.mean([np.min(np.abs(g - z)) for z in TRUE[:K]])

ts = np.concatenate([ -np.logspace(3.5,0.5,25), [0.0], np.logspace(0.5,3.5,25) ])
best = (1e9,None,None)
print("Scanning boundary self-energy t_e (even) x t_o (odd)  [coarse]...")
for te in ts:
    for to in ts:
        g = low_gammas(bcoef_base, K, te, to)
        sc = match_score(g)
        if sc < best[0]:
            best = (sc, te, to, g)
print(f"BEST coarse: score={best[0]:.4f}  t_e={best[1]:.3g}  t_o={best[2]:.3g}")
print(f"  gammas: {np.round(best[3][:8],3)}")
print(f"  TRUE:   {np.round(TRUE,3)}")
gate = best[0] < 0.5
print(f"\nGATE (mean-nearest-dist to first 5 zeros < 0.5): {'PASS' if gate else 'FAIL'}")
if not gate:
    print("=> the boundary self-energy at the last site is NOT the right parametrization of the extension,")
    print("   OR the correct extension needs the full limit-circle phase (2-solution admixture), not a")
    print("   diagonal boundary term. Honest negative — do NOT apply to semilocal until the base gate passes.")
