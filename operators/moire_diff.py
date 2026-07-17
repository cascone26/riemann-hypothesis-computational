#!/usr/bin/env python3
"""
moire_diff.py — PRISM EXTREME item #1: Moiré overlay diff.

Superpose tau_N at two different lsq values (lsq=14 vs lsq=49) cell-by-cell
and analyze the seam (difference matrix). Tony's Moiré principle: "the seam
dominates; emergence is an interface phenomenon." The interference pattern
may carry more signal than either operator alone.

Key questions:
1. What are the eigenvalues of Delta = tau_14 - tau_49 (the seam)?
2. Does Delta have a clean spectral structure (few dominant modes)?
3. What is the spectral overlap between tau_14 and tau_49 eigenvectors?
"""

import math
import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import time
import json

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)

def omega_vec(r):
    return np.real(_digamma(0.25+0.5j*r)) - LOGPI

def build_sieve(lsq):
    sieve = np.ones(lsq+1, dtype=bool)
    sieve[0] = sieve[1] = False
    for i in range(2, int(lsq**0.5)+1):
        if sieve[i]: sieve[i*i::i] = False
    lkpp, wts = [], []
    for p in range(2, lsq+1):
        if not sieve[p]: continue
        pk = p
        while pk <= lsq:
            lkpp.append(math.log(pk))
            wts.append(math.log(p)/math.sqrt(pk))
            pk *= p
    return np.array(lkpp), np.array(wts)

def build_tau(N, lsq, G=None, G_BATCH=2000, K_BATCH=5000):
    L = 2 * math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D = len(ns)
    if G is None:
        G = max(8000, 500*N)

    # Arch
    R = TWO_PI*N/L + 60.0
    grid = np.linspace(-R, R, G)
    h = grid[1]-grid[0]
    w = np.full(G, h); w[0]=w[-1]=0.5*h
    Om = omega_vec(grid)
    arch = np.zeros((D,D))
    for g0 in range(0, G, G_BATCH):
        g1 = min(g0+G_BATCH, G)
        sl = slice(g0, g1)
        c = TWO_PI*ns[:,None]/L
        d = c - grid[None,sl]
        sm = np.abs(d)<1e-12
        E = np.where(sm, np.sqrt(L)+0j,
                     (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
        arch += (E*(w[sl]*Om[sl]) @ E.conj().T).real / TWO_PI
        del E

    # B
    sh = np.sinh(L/4)**2
    ni = ns[:,None]; nj = ns[None,:]
    num = L**2 - 16*np.pi**2*ni*nj
    den = (L**2+16*np.pi**2*ni**2)*(L**2+16*np.pi**2*nj**2)
    B = 32*L*sh*num/den

    # P
    lkpp, wts_p = build_sieve(lsq)
    K = len(lkpp)
    S = np.zeros(D); S_diag = np.zeros(D)
    for k0 in range(0, K, K_BATCH):
        k1 = min(k0+K_BATCH, K)
        ph = TWO_PI*np.outer(ns, lkpp[k0:k1])/L
        S += np.sin(ph) @ wts_p[k0:k1]
        S_diag += np.cos(ph) @ (2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    Si = S[:,None]; Sj = S[None,:]
    dni = ns[:,None] - ns[None,:]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dni != 0, (Sj-Si)/(math.pi*dni), 0.0)
    np.fill_diagonal(P, S_diag)

    tau = arch+B-P
    return (tau+tau.T)/2

if __name__ == "__main__":
    LSQ_A = 14   # lsq_A = 14, L_A ≈ 2.639
    LSQ_B = 49   # lsq_B = 49, L_B ≈ 3.892

    print("Moiré overlay diff: tau_N(lsq=14) vs tau_N(lsq=49)")
    print("=======================================================")
    print(f"L_14 = {2*math.log(math.sqrt(LSQ_A)):.4f},  L_49 = {2*math.log(math.sqrt(LSQ_B)):.4f}")
    print()
    print(f"{'N':>4}  {'D':>4}  {'min_ev(14)':>13}  {'min_ev(49)':>13}  "
          f"{'||Delta||_F':>13}  {'min_ev(Delta)':>14}  {'max_ev(Delta)':>14}  time")
    print("-"*110)

    for N in [4, 6, 8, 10, 12, 16]:
        t0 = time.time()
        # Both operators in same index space (but different L, so indices mean different things)
        # We compare cell-by-cell using the matrix of the same DIMENSION
        tau14 = build_tau(N, LSQ_A, G=6000)
        tau49 = build_tau(N, LSQ_B, G=6000)
        D = 2*N+1

        # Eigenvalues of each
        evs14 = scipy.linalg.eigvalsh(tau14)
        evs49 = scipy.linalg.eigvalsh(tau49)

        # Seam matrix (difference)
        Delta = tau14 - tau49
        frob = np.sqrt(np.sum(Delta**2))
        evD = scipy.linalg.eigvalsh(Delta)

        dt = time.time() - t0
        print(f"{N:>4}  {D:>4}  {evs14[0]:>13.4e}  {evs49[0]:>13.4e}  "
              f"{frob:>13.4e}  {evD[0]:>14.4e}  {evD[-1]:>14.4e}  ({dt:.1f}s)")

    print()
    print("Seam analysis at N=16 (D=33):")
    N = 16
    tau14 = build_tau(N, LSQ_A, G=8000)
    tau49 = build_tau(N, LSQ_B, G=8000)
    Delta = tau14 - tau49
    evD = scipy.linalg.eigvalsh(Delta)
    evsA, vecsA = scipy.linalg.eigh(tau14)
    evsB, vecsB = scipy.linalg.eigh(tau49)

    print(f"  Delta eigenvalues (all {2*N+1}): {evD}")
    print()

    # Subspace overlap: cos²(angle) between near-zero eigenvectors
    # Near-zero subspace of tau14: eigenvectors with |ev| < threshold
    thresh = 1e-3
    idxA = np.where(np.abs(evsA) < thresh)[0]
    idxB = np.where(np.abs(evsB) < thresh)[0]
    print(f"  Near-zero (|ev|<1e-3) dimensions: tau14={len(idxA)}, tau49={len(idxB)}")

    if len(idxA) > 0 and len(idxB) > 0:
        VA = vecsA[:, idxA]
        VB = vecsB[:, idxB]
        overlap = VA.T @ VB  # shape: (|A|, |B|)
        # Singular values of overlap matrix = principal angles
        sv = scipy.linalg.svdvals(overlap)
        print(f"  Subspace principal angles (cos²): {sv**2}")
        print(f"  Max overlap: {sv[0]**2:.6f}  (1.0 = identical subspace)")

    print()
    # Spectral correlation: how correlated are tau14 and tau49 eigenvalues?
    D = 2*N+1
    evsA_full = scipy.linalg.eigvalsh(tau14)
    evsB_full = scipy.linalg.eigvalsh(tau49)
    corr = np.corrcoef(evsA_full, evsB_full)[0,1]
    print(f"  Eigenvalue correlation (tau14 vs tau49): r = {corr:.6f}")

    # Off-diagonal Delta structure
    diag_D = np.abs(np.diag(Delta))
    offdiag_D = np.abs(Delta - np.diag(np.diag(Delta)))
    print(f"  |Delta| diagonal mean: {diag_D.mean():.4e},  off-diagonal mean: {offdiag_D.mean():.4e}")
    print(f"  Delta is dominated by: {'diagonal' if diag_D.mean() > offdiag_D.mean() else 'off-diagonal'} entries")

    # Save results
    result = {
        "LSQ_A": LSQ_A, "LSQ_B": LSQ_B,
        "N_seam": N,
        "eigenvalue_corr_r": float(corr),
        "delta_diag_mean": float(diag_D.mean()),
        "delta_offdiag_mean": float(offdiag_D.mean()),
        "delta_dominant": "diagonal" if diag_D.mean() > offdiag_D.mean() else "off-diagonal",
        "delta_eigenvalues": evD.tolist(),
    }
    out = "/Users/scones/projects/riemann/results/moire_diff.json"
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nSaved → {out}")
