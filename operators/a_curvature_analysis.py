#!/usr/bin/env python3
"""
a_curvature_analysis.py — Characterize the sharpness of a=1/4 as a maximum of lambda_min.

The control experiment showed: ONLY a=1/4 in Omega(r) = Re psi(a + ir/2) - log pi
gives tau_N >= 0. Moving a to any other value makes lambda_min strongly negative.

This script quantifies the CURVATURE:
  d^2 lambda_min / d a^2  at a = 1/4

A sharper maximum (larger |curvature|) means the criticality at a=1/4 is harder to
approach from "nearby" values. If the curvature has an analytic form connected to
known quantities (gamma, log pi, etc.), that might give a proof.

Method: compute lambda_min(tau_N(a)) for a in {0.25 - 2h, 0.25-h, 0.25, 0.25+h, 0.25+2h}
and use central differences to estimate the 1st and 2nd derivatives.
"""

import json, os
import numpy as np
from scipy.special import digamma as _digamma
import mpmath as mp

mp.mp.dps = 30

LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/a_curvature_analysis.json")


def vm_table(kmax):
    tab = []
    for k in range(2, kmax + 1):
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

VM = vm_table(LAMBDA_SQ)


def build_tau_at_a(N, a):
    """Build tau_N with Omega(r) = Re psi(a + ir/2) - log pi."""
    idx = list(range(-N, N+1))
    D = len(idx)
    R = TWO_PI * N / L + 80.0
    G = max(60_000, 3000 * N)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5*h

    # Modified Omega with parameter a
    Om = np.real(_digamma(a + 0.5j * grid)) - LOGPI

    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        s = TWO_PI * n / L
        d = s - grid
        sm = np.abs(d) < 1e-12
        E[k, sm] = np.sqrt(L)
        dd = d[~sm]
        E[k, ~sm] = (1/np.sqrt(L)) * np.exp(0.5j*dd*L) * 2*np.sin(0.5*dd*L) / dd

    arch = ((E * (w * Om)) @ E.conj().T).real / (2*np.pi)

    sh = np.sinh(L/4)**2
    ns = np.array(idx, dtype=float)
    B = 32*L*sh*(L**2 - 16*np.pi**2*np.outer(ns,ns)) / np.outer(L**2+16*np.pi**2*ns**2, L**2+16*np.pi**2*ns**2)

    log_k = np.array([np.log(k) for k, _ in VM])
    wts = np.array([lk * k**(-0.5) for k, lk in VM])
    P = np.zeros((D, D))
    for ii, n in enumerate(idx):
        for jj in range(ii, D):
            m = idx[jj]
            q = np.where(n == m,
                         2*(1-np.abs(log_k)/L)*np.cos(TWO_PI*n*log_k/L),
                         (np.sin(TWO_PI*m*log_k/L)-np.sin(TWO_PI*n*log_k/L))/(np.pi*(n-m)+1e-300))
            P[ii,jj] = P[jj,ii] = float(wts @ q)

    tau = arch + B - P
    return (tau + tau.T)/2


def lambda_min_at_a(N, a):
    """Compute minimum eigenvalue of tau_N(a)."""
    tau = build_tau_at_a(N, a)
    ev = np.linalg.eigvalsh(tau)
    return float(ev[0])


def analyze_curvature(N):
    print(f"\n=== N={N} ===", flush=True)

    # Fine grid around a=1/4 to map the landscape
    a_fine = np.array([0.10, 0.15, 0.20, 0.22, 0.24, 0.25, 0.26, 0.28, 0.30, 0.35, 0.40, 0.50])
    lm_fine = []
    for a in a_fine:
        lm = lambda_min_at_a(N, a)
        lm_fine.append(float(lm))
        print(f"  a={a:.3f}: lambda_min={lm:.6e}", flush=True)

    # Central difference at a=1/4 with two step sizes
    a0 = 0.25
    results_h = []
    for h in [0.005, 0.01, 0.02]:
        lm_m2 = lambda_min_at_a(N, a0 - 2*h)
        lm_m1 = lambda_min_at_a(N, a0 - h)
        lm_0  = lambda_min_at_a(N, a0)
        lm_p1 = lambda_min_at_a(N, a0 + h)
        lm_p2 = lambda_min_at_a(N, a0 + 2*h)

        # 4th order central difference for 2nd derivative
        d2 = (-lm_m2 + 16*lm_m1 - 30*lm_0 + 16*lm_p1 - lm_p2) / (12 * h**2)
        # 1st derivative
        d1 = (-lm_m2 + 8*lm_m1 - 8*lm_p1 + lm_p2) / (12 * h)

        print(f"  h={h}: d1={d1:.4e}  d2={d2:.4e}  lm_0={lm_0:.4e}", flush=True)
        results_h.append({
            "h": h, "d1": float(d1), "d2": float(d2),
            "lm_at_a0": float(lm_0),
            "values": [float(x) for x in [lm_m2, lm_m1, lm_0, lm_p1, lm_p2]]
        })

    return {
        "N": N,
        "a_fine_scan": {"a_values": a_fine.tolist(), "lambda_min": lm_fine},
        "derivative_estimates": results_h,
    }


def main():
    print("A-CURVATURE ANALYSIS", flush=True)
    print(f"lambda^2={LAMBDA_SQ}, L={L:.4f}", flush=True)

    N_values = [16, 32, 48, 64]
    results = []

    for N in N_values:
        r = analyze_curvature(N)
        results.append(r)
        with open(OUTFILE, "w") as f:
            json.dump({"lambda_sq": LAMBDA_SQ, "L": L, "results": results}, f, indent=2)

    # Summary
    print("\n=== CURVATURE SUMMARY ===", flush=True)
    for r in results:
        best_d2 = r["derivative_estimates"][0]["d2"]  # h=0.005
        print(f"  N={r['N']:>4}: d^2 lambda_min/da^2|_{{a=1/4}} ≈ {best_d2:.6e}", flush=True)

    print(f"\nResults: {OUTFILE}", flush=True)


if __name__ == "__main__":
    main()
