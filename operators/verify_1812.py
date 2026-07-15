#!/usr/bin/env python3
"""Quick verify lsq=1812 with G=50000 — checking if positive eigenvalue is real."""
import time, numpy as np
from scipy.special import digamma as _digamma
import scipy.sparse, scipy.sparse.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi; LOGPI = float(np.log(np.pi))

def check(lsq, N, G):
    L = float(mp.log(lsq)); sh = float(np.sinh(L/4)**2)
    tab = []
    for k in range(2, lsq+1):
        n, p, f, ok = k, 0, 2, True
        while f*f<=n:
            if n%f==0:
                if p==0: p=f
                elif p!=f: ok=False; break
                while n%f==0: n//=f
            f+=1
        if n>1:
            if p==0: p=n
            elif p!=n: ok=False
        if ok and p>0: tab.append((k, float(mp.log(p))))
    idx = np.arange(-N, N+1, dtype=float); ns = idx
    R = TWO_PI*N/L + 80.0; grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]; w = np.full(G, h); w[0]=w[-1]=0.5*h
    Om = np.real(_digamma(0.25+0.5j*grid)) - LOGPI
    d = (TWO_PI*ns/L)[:,None] - grid[None,:]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L), (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/(d+1e-300*tiny))
    arch = ((E*(w*Om))@E.conj().T).real/(2.0*np.pi)
    B = (32.0*L*sh*(L**2-16.0*np.pi**2*np.outer(ns,ns))/
         np.outer(L**2+16.0*np.pi**2*ns**2, L**2+16.0*np.pi**2*ns**2))
    if tab:
        lkpp = np.array([float(mp.log(k)) for k,_ in tab])
        wts = np.array([lk*k**(-0.5) for k,lk in tab])
        ph = TWO_PI*np.outer(ns,lkpp)/L; S = np.sin(ph)@wts
        dn = ns[:,None]-ns[None,:]
        with np.errstate(divide='ignore',invalid='ignore'):
            P = (S[None,:]-S[:,None])/(np.pi*dn)
        P[dn==0]=0.0; np.fill_diagonal(P, np.cos(ph)@(2*wts*(1-lkpp/L)))
    else:
        P = np.zeros((len(idx),len(idx)))
    tau = (arch+B-P); tau=(tau+tau.T)/2.0
    tau_sp = scipy.sparse.csr_matrix(tau)
    t0 = time.time()
    try:
        ev, _ = scipy.sparse.linalg.eigsh(tau_sp, k=1, which='SA', tol=1e-8, maxiter=10000)
        lam = float(ev[0])
    except Exception:
        # Dense fallback
        ev2 = scipy.linalg.eigvalsh(tau)
        lam = float(ev2[0])
    dt = time.time()-t0
    return lam, dt

print("Verifying lsq=1812 positive eigenvalue claim")
for G in [30000, 50000]:
    lam, dt = check(1812, 400, G)
    sign = "POSITIVE!" if lam > 0 else "negative"
    print(f"  G={G:5d}: λ={lam:.6e}  ({dt:.0f}s)  [{sign}]")

# Also check neighbors
print("\nNeighbors at G=50000:")
for lsq in [1808, 1810, 1812, 1814, 1816]:
    lam, dt = check(lsq, 400, 50000)
    sign = "POSITIVE!" if lam > 0 else ""
    print(f"  lsq={lsq}: λ={lam:.6e}  ({dt:.0f}s) {sign}")
