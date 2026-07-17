#!/usr/bin/env python3
"""
near_miss_tracker.py — PRISM item #4: Systematic closest-near-miss tracker.

Adversarial search: for non-Riemann lsq values (γ̃ NOT a Riemann zero),
what is the MINIMUM min_ev achievable, and at what N does it occur?

This directly tests: can a non-Riemann lsq achieve min_ev SMALLER than the
Riemann Dip1 for the SAME j, and if so, at what N relative to the Riemann
resonant N?

Method: for j=20 (small, fast), sweep:
  - Riemann lsq: round(exp(2π×20/γ₁)) = 154,552 → Riemann Dip1 at N≈2200
  - Non-Riemann lsq: grid of γ̃ values near γ₁ but NOT equal

For each γ̃, compute the MIN of min_ev over N ∈ {500, 1000, 1500, 2000, 2500, 3000}
and record: (γ̃, min_lsq, N_at_min, min_ev_achieved, vs_Riemann_Dip1)

Key hypothesis: for non-Riemann γ̃, ANY near-miss (min_ev ≈ Riemann Dip1) occurs
at a DIFFERENT N than the Riemann resonant N=2200.

Uses the fast streaming formula (batched arch + P, same as Limit Theorem scans).
"""

import math
import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import time

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
G      = 10000   # arch quadrature (reduced for speed at small N)
G_BATCH= 1000

RIEMANN_DIP1_j20 = 1.136e-7  # confirmed Dip1 at N=2200 for j=20
GAMMA1 = 14.13472514

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

def build_tau(N, lsq, gamma_tilde, G_arch=4000):
    L  = 2 * math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D  = len(ns)

    # Arch (reduced G for speed in sweep context)
    R    = TWO_PI*N/L + 40.0
    grid = np.linspace(-R, R, G_arch)
    h    = grid[1]-grid[0]
    w    = np.full(G_arch, h); w[0]=w[-1]=0.5*h
    Om   = omega_vec(grid)
    c    = TWO_PI*ns[:,None]/L
    d    = c - grid[None,:]
    small= np.abs(d)<1e-12
    E    = np.where(small, np.sqrt(L)+0j,
                    (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch = ((E*(w*Om)) @ E.conj().T).real/(TWO_PI)
    del E

    # B
    sh  = np.sinh(L/4)**2
    ni  = ns[:,None]; nj = ns[None,:]
    num = L**2 - 16*np.pi**2*ni*nj
    den = (L**2+16*np.pi**2*ni**2)*(L**2+16*np.pi**2*nj**2)
    B   = 32*L*sh*num/den

    # P (batched sieve, vectorized assembly)
    lkpp, wts = build_sieve(lsq)
    K      = len(lkpp)
    S      = np.zeros(D)
    S_diag = np.zeros(D)
    for k0 in range(0, K, G_BATCH):
        k1   = min(k0+G_BATCH, K)
        lk_b = lkpp[k0:k1]; wt_b = wts[k0:k1]
        ph   = TWO_PI*np.outer(ns, lk_b)/L
        S      += np.sin(ph) @ wt_b
        S_diag += np.cos(ph) @ (2*wt_b*(1-lk_b/L))
        del ph
    # Vectorized P assembly: P[i,j] = (S[i]-S[j])/(pi*(ns[i]-ns[j])) for i≠j
    Si  = S[:,None]; Sj  = S[None,:]
    dni = ns[:,None] - ns[None,:]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dni != 0, (Sj-Si)/(math.pi*dni), 0.0)
    np.fill_diagonal(P, S_diag)

    tau = arch+B-P
    return (tau+tau.T)/2

def run_sweep(j, gamma_tilde, N_list):
    lsq = round(math.exp(TWO_PI*j/gamma_tilde))
    results = []
    for N in N_list:
        t0 = time.time()
        tau = build_tau(N, lsq, gamma_tilde)
        evs = scipy.linalg.eigvalsh(tau, subset_by_index=[0,1])
        me  = min(abs(evs[0]), abs(evs[1]))
        results.append((N, me, time.time()-t0))
    return lsq, results

if __name__ == "__main__":
    j = 20
    N_LIST = [500, 1000, 1500, 2000]

    # Riemann baseline
    gamma1 = GAMMA1
    lsq_R  = round(math.exp(TWO_PI*j/gamma1))  # 154,552

    print(f"Near-miss tracker  j={j}  Riemann lsq={lsq_R}  Dip1={RIEMANN_DIP1_j20:.3e} at N=2200")
    print(f"γ̃ sweep: testing non-Riemann values near γ₁={gamma1:.5f}")
    print(f"N sweep: {N_LIST}")
    print()

    # Test γ̃ values: Riemann ± offsets
    gamma_tildes = [
        gamma1,          # exact Riemann (baseline)
        gamma1 - 0.3,    # below γ₁
        gamma1 - 0.1,
        gamma1 + 0.1,    # above γ₁
        gamma1 + 0.3,
        13.8,            # backward direction (tested in j=40 campaign)
        15.0,            # significantly off
        21.022,          # exactly γ₂ (different zero, different resonant N)
    ]

    print(f"{'γ̃':>8} {'lsq':>10} {'N_min':>7} {'min_ev':>12} {'vs_Dip1':>10} {'sub_Riemann?':>13}")
    print("-"*75)

    all_results = []
    for gt in gamma_tildes:
        lsq, sweep = run_sweep(j, gt, N_LIST)
        best_N, best_me, _ = min(sweep, key=lambda x: x[1])
        ratio = best_me / RIEMANN_DIP1_j20
        sub   = "YES <<<" if best_me < RIEMANN_DIP1_j20 else ""
        is_riemann = abs(gt - gamma1) < 0.001
        label = " ← RIEMANN" if is_riemann else ""
        print(f"{gt:>8.4f} {lsq:>10} {best_N:>7} {best_me:>12.3e} {ratio:>10.3f}×  {sub}{label}")
        all_results.append({"gamma": gt, "lsq": lsq, "best_N": best_N,
                            "best_me": best_me, "ratio": ratio})

    print()
    print("Near-miss candidates (min_ev < 10× Riemann Dip1):")
    near = [r for r in all_results if r["ratio"] < 10 and abs(r["gamma"]-gamma1)>0.001]
    if near:
        for r in near:
            print(f"  γ̃={r['gamma']:.4f}: min_ev={r['best_me']:.3e} at N={r['best_N']} ({r['ratio']:.3f}× Dip1)")
    else:
        print("  None — all non-Riemann cases ≥ 10× Riemann Dip1 at every N tested.")
    print()
    print("WHERE discrimination: Riemann Dip1 at N=2200; non-Riemann min at different N?")
    for r in all_results:
        if r["ratio"] < 100:
            print(f"  γ̃={r['gamma']:.4f}: deepest dip at N={r['best_N']} (Riemann at N=2200)")
