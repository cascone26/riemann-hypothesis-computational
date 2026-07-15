#!/usr/bin/env python3
"""
ev_trajectory.py — Track individual eigenvalue trajectories through k=4..5.

Computes the 3 smallest |eigenvalues| of tau as lambda_sq increases from 4.0 to 5.0
in very fine steps (0.002). Determines:
1. Whether crossing at ~k+0.75 is exactly 3/4 or threshold-dependent
2. The actual shape of the eigenvalue curve (where does it hit zero vs 0.02)
3. Whether the curve is monotone or has a minimum
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
N      = 64

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/ev_trajectory.json")


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


def build_ev5(lambda_sq):
    """Return 5 smallest eigenvalues of tau at lambda_sq."""
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
    tau = arch + B - P
    tau = (tau + tau.T) / 2.0
    ev = scipy.linalg.eigvalsh(tau)
    # Return 5 smallest (may be negative at machine prec)
    order = np.argsort(np.abs(ev))
    return [float(ev[order[i]]) for i in range(5)]


def main():
    print("Eigenvalue trajectories through k=4..5 and k=5..6", flush=True)

    results = {}
    for k in [4, 5]:
        lsq_arr = np.arange(float(k), float(k) + 1.01, 0.002)
        pts = []
        print(f"\nk={k}: scanning {len(lsq_arr)} points...", flush=True)
        t0 = time.time()
        for lsq in lsq_arr:
            ev5 = build_ev5(lsq)
            pts.append({
                'lsq': round(float(lsq), 4),
                'ev': [round(x, 8) for x in ev5],
            })
            # Print progress at key values
            if abs(lsq - (k + 0.25)) < 0.002 or abs(lsq - (k + 0.50)) < 0.002 \
               or abs(lsq - (k + 0.75)) < 0.002 or abs(lsq - (k + 1.0)) < 0.002 \
               or abs(lsq - k) < 0.002:
                ev5s = ' '.join(f'{x:.4e}' for x in ev5)
                print(f"  lsq={lsq:.3f}: [{ev5s}]", flush=True)
        dt = time.time() - t0

        # Find where each eigenvalue trajectory crosses 0.02 (threshold) vs 0 (true)
        print(f"\n  Threshold crossings for k={k}:")
        for thresh in [0.05, 0.02, 0.005, 0.001]:
            # Count how many |ev[0]| < thresh at each point
            transitions = []
            prev_count = sum(1 for x in pts[0]['ev'] if abs(x) < thresh)
            for pt in pts[1:]:
                count = sum(1 for x in pt['ev'] if abs(x) < thresh)
                if count != prev_count:
                    transitions.append((pt['lsq'], prev_count, count))
                    prev_count = count
            offsets = [round(lsq - k, 4) for lsq, _, _ in transitions]
            print(f"  thresh={thresh:.3f}: transitions at offsets={offsets}")

        results[k] = pts
        print(f"  Done in {dt:.1f}s")

    with open(OUTFILE, 'w') as f:
        json.dump({'N': N, 'G': G, 'trajectories': results}, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
