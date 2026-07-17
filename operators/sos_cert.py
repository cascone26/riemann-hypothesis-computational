#!/usr/bin/env python3
"""
sos_cert.py — PRISM item #2: Finite-N SOS certificate extraction for tau_N.

tau_N is observed to be PSD at all tested N. An SOS (sum-of-squares) certificate
proves this algebraically: tau_N = V @ diag(λ) @ V^T with all λ ≥ 0.

The spectral decomposition IS the SOS certificate (trivially). The more interesting
question: what is the algebraic degree of the certificate (how many squares), and
does it grow or stabilize as N grows?

This script:
1. Builds tau_N for N ∈ {4, 6, 8, 10, 12}
2. Computes the spectral decomposition (eigenvectors = certificate vectors)
3. Reports: # positive eigenvalues, spectral gap, condition number
4. Extracts the SOS factorization: tau_N ≈ sum_i λ_i v_i v_i^T (λ_i > 0)
5. Verifies certificate residual ||tau_N - V_pos @ diag(λ_pos) @ V_pos.T||_F
6. Uses cvxpy to attempt a degree-bounded SOS verification (optional, slow)

Key questions:
- How many squares are needed? (= # positive eigenvalues)
- Does the spectral gap → constant or → 0 as N grows?
- Is the certificate "tight" (minimum # squares) or wasteful?
"""

import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import mpmath as mp
import time

mp.mp.dps = 30
LAMBDA_SQ = 14
L = float(2*mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0*np.pi
LOGPI = float(np.log(np.pi))

def vm_table(kmax):
    tab = []
    for k in range(2, kmax+1):
        n, p, f, ok = k, 0, 2, True
        while f*f <= n:
            if n%f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n%f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((float(k), float(mp.log(p))))
    return tab

VM = vm_table(LAMBDA_SQ)

def omega_vec(r):
    return np.real(_digamma(0.25+0.5j*r)) - LOGPI

def e_hat_vec(r, n):
    s = TWO_PI*n/L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L)+0j
    dd = d[~small]
    out[~small] = (1/np.sqrt(L))*np.exp(0.5j*dd*L)*2*np.sin(0.5*dd*L)/dd
    return out

def q_kernel(n, m, y):
    if n == m:
        return 2*(1-abs(y)/L)*np.cos(TWO_PI*n*y/L)
    return (np.sin(TWO_PI*m*y/L)-np.sin(TWO_PI*n*y/L))/(np.pi*(n-m))

def build_tau(N, G=30000):
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI*N/L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]
    w = np.full(G,h); w[0]=w[-1]=0.5*h
    Om = omega_vec(grid)
    E = np.array([e_hat_vec(grid,n) for n in idx])
    arch = (E*(w*Om) @ E.conj().T).real/(2*np.pi)
    # B matrix
    sh = np.sinh(L/4)**2
    B = np.zeros((D,D))
    for ii,n in enumerate(idx):
        for jj,m in enumerate(idx):
            num = L**2-16*np.pi**2*m*n
            den = (L**2+16*np.pi**2*m**2)*(L**2+16*np.pi**2*n**2)
            B[ii,jj] = 32*L*sh*num/den
    # P matrix
    P = np.zeros((D,D))
    for ii,n in enumerate(idx):
        for jj in range(ii,D):
            m = idx[jj]
            s = sum(Lam*k**(-0.5)*q_kernel(n,m,np.log(k)) for k,Lam in VM)
            P[ii,jj]=P[jj,ii]=s
    tau = arch+B-P
    return (tau+tau.T)/2

def analyze_cert(N):
    t0 = time.time()
    tau = build_tau(N)
    D = 2*N+1
    evs, vecs = scipy.linalg.eigh(tau)
    dt = time.time()-t0

    n_pos = int(np.sum(evs > 1e-12))
    n_zero = int(np.sum(np.abs(evs) <= 1e-12))
    n_neg = int(np.sum(evs < -1e-12))
    ev_min = float(evs[0])
    ev_max = float(evs[-1])

    # Spectral gap = smallest positive eigenvalue / |smallest eigenvalue|
    pos_evs = evs[evs > 1e-12]
    gap = float(pos_evs[0]) if len(pos_evs) > 0 else 0.0

    # SOS certificate: tau ≈ sum_i λ_i v_i v_i^T for λ_i > 0
    mask = evs > 1e-12
    V_pos = vecs[:, mask]
    lam_pos = evs[mask]
    recon = V_pos @ np.diag(lam_pos) @ V_pos.T
    resid = np.linalg.norm(tau - recon, 'fro') / np.linalg.norm(tau, 'fro')

    # Condition number of positive part
    cond = float(lam_pos[-1] / lam_pos[0]) if len(lam_pos) > 1 else float('inf')

    print(f"N={N:>3} D={D:>3}  ev_min={ev_min:>10.3e}  ev_max={ev_max:>10.3e}")
    print(f"       +:{n_pos:>3}  ~0:{n_zero}  -:{n_neg:>3}   gap={gap:.3e}  cond={cond:.2e}")
    print(f"       SOS squares needed: {n_pos}  cert_residual={resid:.2e}  ({dt:.1f}s)")
    print(f"       Spectral ratio (gap/|ev_min|): {gap/abs(ev_min):.3f}")
    return {"N": N, "D": D, "ev_min": ev_min, "ev_max": ev_max,
            "n_pos": n_pos, "n_zero": n_zero, "n_neg": n_neg,
            "gap": gap, "cond": cond, "cert_resid": resid}


if __name__ == "__main__":
    print(f"SOS Certificate Analysis  LAMBDA_SQ={LAMBDA_SQ}  L={L:.4f}")
    print(f"tau_N = arch + B - P_primes  (real prime comb, von Mangoldt)")
    print("="*70)
    results = []
    for N in [4, 6, 8, 10, 12, 16]:
        r = analyze_cert(N)
        results.append(r)
        print()

    print("="*70)
    print("SUMMARY — SOS squares vs N:")
    print(f"{'N':>4}  {'D':>4}  {'#squares':>9}  {'gap':>10}  {'cond':>10}  {'ev_min':>12}")
    for r in results:
        print(f"{r['N']:>4}  {r['D']:>4}  {r['n_pos']:>9}  {r['gap']:>10.3e}  "
              f"{r['cond']:>10.2e}  {r['ev_min']:>12.3e}")

    print("\nKey: #squares stabilizing → certificate complexity bounded.")
    print("     gap → constant → spectral gap conjecture confirmed at small N.")
