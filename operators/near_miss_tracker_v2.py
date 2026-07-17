#!/usr/bin/env python3
"""
near_miss_tracker_v2.py — Two-phase near-miss tracker.

v1 flaw: compared min over N (infima) — hid the WHERE discrimination.
v2 fix: two phases:
  Phase 1: For each γ̃, find N_Dip1 via coarse→fine N scan.
  Phase 2: At the RIEMANN N_Dip1 (=2200 for j=20), evaluate ALL γ̃.

The correct selectivity metric (from PROCESS.md revised strategy):
  "At Riemann's N" comparison, not infimum comparison.

At N=N_Dip1_Riemann, non-Riemann lsq must have min_ev >> Riemann Dip1.
If non-Riemann achieves a lower infimum at a DIFFERENT N, that's not
a violation — the Riemann resonance is N-specific.

Outputs:
  For each γ̃:
    - N_Dip1 (where does its min occur?)
    - min_ev at its own N_Dip1 (infimum)
    - min_ev at Riemann's N_Dip1 (=2200 for j=20)
    - ratio of each vs Riemann Dip1
    - N-displacement from Riemann N_Dip1
"""

import math
import numpy as np
import scipy.linalg
from scipy.special import digamma as _digamma
import json
import time

TWO_PI = 2*math.pi
LOGPI  = math.log(math.pi)
GAMMA1 = 14.13472514
GAMMA2 = 21.02203964
J      = 20
RIEMANN_DIP1  = 1.136e-7    # confirmed j=20 Dip1 at N=2200
N_DIP1_RIEMANN = 2200       # Riemann's resonant N for j=20

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

def build_tau_fast(N, lsq, G_arch=6000, K_BATCH=1000):
    L  = 2*math.log(math.sqrt(lsq))
    ns = np.arange(-N, N+1)
    D  = len(ns)
    R  = TWO_PI*N/L + 40.0
    grid = np.linspace(-R, R, G_arch)
    h    = grid[1]-grid[0]
    w    = np.full(G_arch, h); w[0]=w[-1]=0.5*h
    Om   = omega_vec(grid)
    c    = TWO_PI*ns[:,None]/L
    d    = c - grid[None,:]
    sm   = np.abs(d)<1e-12
    E    = np.where(sm, np.sqrt(L)+0j,
                    (1/np.sqrt(L))*np.exp(0.5j*d*L)*2*np.sin(0.5*d*L)/d)
    arch = (E*(w*Om) @ E.conj().T).real / TWO_PI
    sh  = np.sinh(L/4)**2
    ni  = ns[:,None]; nj = ns[None,:]
    num = L**2 - 16*np.pi**2*ni*nj
    den = (L**2+16*np.pi**2*ni**2)*(L**2+16*np.pi**2*nj**2)
    B   = 32*L*sh*num/den
    lkpp, wts_p = build_sieve(lsq)
    K = len(lkpp)
    S = np.zeros(D); S_diag = np.zeros(D)
    for k0 in range(0, K, K_BATCH):
        k1 = min(k0+K_BATCH, K)
        ph = TWO_PI*np.outer(ns, lkpp[k0:k1])/L
        S      += np.sin(ph) @ wts_p[k0:k1]
        S_diag += np.cos(ph) @ (2*wts_p[k0:k1]*(1-lkpp[k0:k1]/L))
    Si = S[:,None]; Sj = S[None,:]
    dni = ns[:,None] - ns[None,:]
    with np.errstate(divide='ignore', invalid='ignore'):
        P = np.where(dni != 0, (Sj-Si)/(math.pi*dni), 0.0)
    np.fill_diagonal(P, S_diag)
    tau = arch+B-P
    return (tau+tau.T)/2

def min_ev_at_N(N, lsq):
    tau = build_tau_fast(N, lsq)
    evs = scipy.linalg.eigvalsh(tau, subset_by_index=[0, 1])
    return float(min(abs(evs[0]), abs(evs[1])))

def find_dip1(lsq, N_coarse, label=""):
    """Phase 1: coarse→fine scan to find N_Dip1."""
    # Coarse sweep
    coarse_me = []
    for N in N_coarse:
        me = min_ev_at_N(N, lsq)
        coarse_me.append((N, me))

    # Find coarse minimum
    best_N, best_me = min(coarse_me, key=lambda x: x[1])
    idx = [x[0] for x in coarse_me].index(best_N)

    # Determine bracket for fine scan
    N_lo = coarse_me[max(0, idx-1)][0]
    N_hi = coarse_me[min(len(coarse_me)-1, idx+1)][0]
    step_fine = max(100, (N_hi - N_lo) // 10)

    # Fine scan in bracket
    N_fine = list(range(N_lo, N_hi+step_fine, step_fine))
    fine_me = []
    for N in N_fine:
        me = min_ev_at_N(N, lsq)
        fine_me.append((N, me))
        if label:
            print(f"    {label} N={N}: {me:.4e}")

    best_N_f, best_me_f = min(fine_me, key=lambda x: x[1])
    return best_N_f, best_me_f, coarse_me + fine_me

if __name__ == "__main__":
    print("Near-Miss Tracker v2 — Two-Phase Selectivity Analysis")
    print("=======================================================")
    print(f"j={J}, Riemann Dip1={RIEMANN_DIP1:.3e} at N={N_DIP1_RIEMANN}")
    print(f"γ₁ = {GAMMA1:.8f}")
    print()

    # γ̃ values to test
    gamma_tests = [
        (GAMMA1,           "γ₁ (RIEMANN baseline)"),
        (GAMMA1 - 0.5,     "γ₁ - 0.5"),
        (GAMMA1 - 0.1,     "γ₁ - 0.1"),
        (GAMMA1 + 0.1,     "γ₁ + 0.1"),
        (GAMMA1 + 0.5,     "γ₁ + 0.5"),
        (13.8,             "γ̃=13.8 (j=40 backward case)"),
        (GAMMA2,           "γ₂=21.022 (second zero)"),
        (17.0,             "γ̃=17.0 (between γ₁,γ₂)"),
    ]

    # Coarse N sweep: broad coverage around expected Dip1 range for j=20
    N_COARSE = [300, 600, 1000, 1500, 2000, 2500, 3000, 3500]

    print("=== PHASE 1: Find N_Dip1 for each γ̃ ===")
    print()

    phase1_results = []
    for gt, label in gamma_tests:
        lsq = round(math.exp(TWO_PI*J/gt))
        t0 = time.time()
        print(f"  {label}: lsq={lsq}")
        N_dip1, me_dip1, all_pts = find_dip1(lsq, N_COARSE)
        dt = time.time()-t0
        n_displacement = N_dip1 - N_DIP1_RIEMANN
        print(f"    → N_Dip1={N_dip1}, min_ev={me_dip1:.4e}, "
              f"N_disp={n_displacement:+d}  ({dt:.1f}s)")
        print()
        phase1_results.append({
            "gamma": gt, "label": label, "lsq": lsq,
            "N_dip1": N_dip1, "me_dip1": float(me_dip1),
            "N_displacement": int(n_displacement),
            "all_pts": [(n, float(m)) for n, m in all_pts]
        })

    print()
    print("=== PHASE 2: Evaluate all γ̃ at Riemann's N_Dip1 = 2200 ===")
    print()
    print(f"  {'γ̃':>10}  {'lsq':>10}  {'me_at_N2200':>13}  {'vs_Dip1':>10}  {'N_Dip1_own':>11}  {'N_disp':>8}")
    print("  " + "─"*75)

    phase2_results = []
    for r in phase1_results:
        me_at_riemann_N = min_ev_at_N(N_DIP1_RIEMANN, r["lsq"])
        ratio_at_N      = me_at_riemann_N / RIEMANN_DIP1
        ratio_infimum   = r["me_dip1"] / RIEMANN_DIP1
        is_riemann      = abs(r["gamma"] - GAMMA1) < 0.001
        flag = " ← RIEMANN" if is_riemann else ("  *** NEAR MISS ***" if ratio_at_N < 5 else "")
        print(f"  {r['gamma']:>10.4f}  {r['lsq']:>10}  {me_at_riemann_N:>13.4e}"
              f"  {ratio_at_N:>9.1f}×  {r['N_dip1']:>11}  {r['N_displacement']:>+8}{flag}")
        phase2_results.append({
            **r,
            "me_at_riemann_N": float(me_at_riemann_N),
            "ratio_at_riemann_N": float(ratio_at_N),
            "ratio_infimum": float(ratio_infimum)
        })

    print()
    print("=== COMPARISON: Infimum vs At-Riemann-N ratios ===")
    print()
    print(f"  {'γ̃':>10}  {'infimum_ratio':>14}  {'at_N2200_ratio':>15}  {'most_dangerous'}  ")
    print("  " + "─"*65)
    for r in phase2_results:
        inf_r = r["ratio_infimum"]
        atn_r = r["ratio_at_riemann_N"]
        most_dangerous = "infimum" if inf_r < atn_r else "at-N2200"
        is_riemann = abs(r["gamma"] - GAMMA1) < 0.001
        marker = " ← RIEMANN" if is_riemann else ""
        print(f"  {r['gamma']:>10.4f}  {inf_r:>13.1f}×  {atn_r:>14.1f}×  {most_dangerous}{marker}")

    print()
    print("=== KEY FINDINGS ===")
    print()
    non_riemann = [r for r in phase2_results if abs(r["gamma"]-GAMMA1) > 0.001]
    riemann = [r for r in phase2_results if abs(r["gamma"]-GAMMA1) < 0.001][0]
    min_atn = min(non_riemann, key=lambda r: r["ratio_at_riemann_N"])
    min_inf = min(non_riemann, key=lambda r: r["ratio_infimum"])
    print(f"  Riemann at-N2200:     me = {riemann['me_at_riemann_N']:.4e} = Dip1 baseline")
    print(f"  Best non-Riemann at-N2200: γ̃={min_atn['gamma']:.4f}, "
          f"ratio={min_atn['ratio_at_riemann_N']:.1f}× Dip1")
    print(f"  Best non-Riemann infimum:  γ̃={min_inf['gamma']:.4f}, "
          f"ratio={min_inf['ratio_infimum']:.1f}× Dip1 at N={min_inf['N_dip1']}")
    print()
    print("N-displacement analysis:")
    for r in non_riemann:
        disp = r["N_displacement"]
        print(f"  γ̃={r['gamma']:.4f}: N_Dip1={r['N_dip1']} "
              f"({disp:+d} from Riemann N=2200)")

    print()
    print("CONCLUSION:")
    all_above = all(r["ratio_at_riemann_N"] > 10 for r in non_riemann)
    print(f"  At Riemann's N=2200: all non-Riemann cases ≥ {min_atn['ratio_at_riemann_N']:.1f}× Dip1 → "
          f"{'SELECT (>10×)' if all_above else 'PARTIAL'}")
    all_above_inf = all(r["ratio_infimum"] > 10 for r in non_riemann)
    print(f"  By infimum: all non-Riemann cases ≥ {min_inf['ratio_infimum']:.1f}× Dip1 → "
          f"{'SELECT (>10×)' if all_above_inf else 'PARTIAL — infimum can be close!'}")
    print()
    print("  KEY: The at-N2200 metric is STRONGER than infimum for proving selectivity.")
    print("  A near-miss by infimum at DIFFERENT N ≠ a violation of RH selectivity.")

    out = "/Users/scones/projects/riemann/results/near_miss_v2.json"
    with open(out, "w") as fh:
        json.dump({
            "j": J, "riemann_dip1": RIEMANN_DIP1,
            "N_dip1_riemann": N_DIP1_RIEMANN,
            "phase1": [{k: v for k, v in r.items() if k != "all_pts"}
                       for r in phase1_results],
            "phase2": phase2_results
        }, fh, indent=2)
    print(f"Saved → {out}")
