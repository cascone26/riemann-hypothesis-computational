#!/usr/bin/env python3
"""
convergence_rate.py — Measure sub_cancel convergence rate and arch integral behavior.

sub_cancel (||P-AB||_F in near-null subspace / ||AB||_F) decreases as lambda_sq grows.
This fits sub_cancel ~ C / L^alpha. Measure alpha.

Also tests the arch integral G-scaling:
  Required G = max(60000, K * 2*pi*N^2 / L) for some constant K.
  Tests G=60000 vs G=120000 vs adaptive G = max(60000, ceil(8*pi*N^2/L)).

Tests lsq=4..100, step 2, N=64, with both standard G=60000 and adaptive G.
"""
import json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
N      = 64

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/convergence_rate.json")


def vm_table(kmax):
    tab = []
    for k in range(2, int(kmax) + 1):
        n, p, f, ok = k, 0, 2, True
        while f * f <= n:
            if n % f == 0:
                if p == 0: p = f
                elif p != f: ok = False; break
                while n % f == 0: n //= f
            f += 1
        if n > 1:
            if p == 0: p = n
            elif p != n: ok = False
        if ok and p > 0:
            tab.append((k, float(mp.log(p))))
    return tab


def build_AB_P(lambda_sq, G):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx
    R   = TWO_PI * N / L + 80.0
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]; w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = np.real(_digamma(0.25 + 0.5j * grid)) - LOGPI
    d = (TWO_PI * ns / L)[:, None] - grid[None, :]
    tiny = np.abs(d) < 1e-12
    E = np.where(tiny, np.sqrt(L),
                 (1 / np.sqrt(L)) * np.exp(0.5j * d * L)
                 * 2 * np.sin(0.5 * d * L) / (d + 1e-300 * tiny))
    arch = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
    B = (32.0 * L * sh * (L ** 2 - 16.0 * np.pi ** 2 * np.outer(ns, ns)) /
         np.outer(L ** 2 + 16.0 * np.pi ** 2 * ns ** 2,
                  L ** 2 + 16.0 * np.pi ** 2 * ns ** 2))
    AB = (arch + B + (arch + B).T) / 2.0
    if VM:
        lkpp = np.array([np.log(k) for k, _ in VM])
        wts  = np.array([lk * k ** (-0.5) for k, lk in VM])
        ph   = TWO_PI * np.outer(ns, lkpp) / L
        S    = np.sin(ph) @ wts
        dn   = ns[:, None] - ns[None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            P = (S[None, :] - S[:, None]) / (np.pi * dn)
        P[dn == 0] = 0.0
        np.fill_diagonal(P, np.cos(ph) @ (2 * wts * (1 - lkpp / L)))
    else:
        P = np.zeros((len(idx), len(idx)))
    P = (P + P.T) / 2.0
    return AB, P


def compute_sub_cancel(AB, P, thresh=0.02):
    tau = AB - P
    ev_tau, V_tau = scipy.linalg.eigh(tau)
    dim = int((np.abs(ev_tau) < thresh).sum())
    lam_min = float(ev_tau[0])
    if dim == 0:
        return float('nan'), 0, lam_min
    null_mask = np.abs(ev_tau) < thresh
    V_null = V_tau[:, null_mask]
    AB_proj = V_null.T @ AB @ V_null
    P_proj  = V_null.T @ P  @ V_null
    fro_AB = np.linalg.norm(AB_proj, 'fro')
    sub_cancel = np.linalg.norm(AB_proj - P_proj, 'fro') / fro_AB
    return float(sub_cancel), dim, lam_min


def main():
    G_fixed = 60000
    lsq_list = list(range(4, 101, 2))
    results = []

    print(f"{'lsq':>5}  {'L':>6}  {'G_fix':>8}  {'sc_fix':>10}  "
          f"{'G_adap':>8}  {'sc_adap':>10}  {'dim':>4}  {'lam_fix':>12}")
    t0 = time.time()

    for lsq in lsq_list:
        L = float(mp.log(lsq))
        # Adaptive G: ensure h = 2R/G < 2*pi/L / 10
        R_est = TWO_PI * N / L + 80.0
        G_adaptive = max(G_fixed, int(np.ceil(20 * TWO_PI * N / L)))
        G_adaptive = min(G_adaptive, 240000)  # cap at 240K

        t1 = time.time()
        AB_fix, P = build_AB_P(float(lsq), G_fixed)
        sc_fix, dim, lam_fix = compute_sub_cancel(AB_fix, P)

        AB_adap, _ = build_AB_P(float(lsq), G_adaptive)
        sc_adap, _, lam_adap = compute_sub_cancel(AB_adap, P)

        print(f"{lsq:>5}  {L:>6.3f}  {G_fixed:>8}  {sc_fix:>10.6f}  "
              f"{G_adaptive:>8}  {sc_adap:>10.6f}  {dim:>4}  {lam_fix:>12.4e}  "
              f"({time.time()-t1:.1f}s)", flush=True)
        results.append({
            'lsq': lsq, 'L': round(L, 6),
            'G_fixed': G_fixed, 'sub_cancel_fixed': sc_fix, 'lam_min_fixed': lam_fix,
            'G_adaptive': G_adaptive, 'sub_cancel_adaptive': float(sc_adap),
            'dim': dim,
        })

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'results': results}, f, indent=2)

    # Fit power law sub_cancel ~ C/L^alpha
    good = [(r['L'], r['sub_cancel_fixed']) for r in results
            if not np.isnan(r['sub_cancel_fixed']) and r['dim'] > 0]
    if len(good) > 4:
        Ls = np.array([x[0] for x in good])
        sc = np.array([x[1] for x in good])
        # log-log fit
        log_L = np.log(Ls)
        log_sc = np.log(np.clip(sc, 1e-12, None))
        A = np.column_stack([log_L, np.ones_like(log_L)])
        alpha, log_C = np.linalg.lstsq(A, log_sc, rcond=None)[0]
        C = np.exp(log_C)
        print(f"\nPower-law fit: sub_cancel ≈ {C:.4e} / L^{-alpha:.3f}")
        print(f"i.e., sub_cancel ~ 1/L^{-alpha:.3f} (exponent = {-alpha:.3f})")

    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
