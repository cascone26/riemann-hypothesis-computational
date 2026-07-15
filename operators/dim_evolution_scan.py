#!/usr/bin/env python3
"""
dim_evolution_scan.py — Fine-grained scan of near-null dim vs lambda^2.

Step size 0.05 from lambda^2 = 4 to 22.
Questions:
  1. Does dim jump by exactly 2 at each integer lambda^2?
  2. Do jumps happen AT the integer or slightly before/after?
  3. Is there any signature at Riemann zeros (gamma_1=14.135, gamma_2=21.022)?
  4. What is the exact eigenvalue of the "25th mode" as lambda^2 varies near 14?

Also tracks:
  - The 5 smallest-magnitude eigenvalues at each lambda^2
  - dim at three thresholds: 1e-4, 1e-3, 1e-2
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
N      = 96
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/dim_evolution_scan.json")

ZEROS = [14.134725141734693, 21.022039638771555, 25.010857580145688]


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


def build_tau(lambda_sq, N):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
    D, ns = len(idx), idx
    R = TWO_PI * N / L + 80.0
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
        P = np.zeros((D, D))
    tau = arch + B - P
    return (tau + tau.T) / 2.0


def analyze(lsq):
    t0  = time.time()
    tau = build_tau(lsq, N)
    ev  = scipy.linalg.eigvalsh(tau)
    dt  = time.time() - t0

    # Dims at three thresholds
    dims = {'1e-04': int((np.abs(ev) < 1e-4).sum()),
            '0.001': int((np.abs(ev) < 1e-3).sum()),
            '0.01':  int((np.abs(ev) < 1e-2).sum())}

    # Five smallest-magnitude eigenvalues
    order = np.argsort(np.abs(ev))
    top5  = [float(ev[i]) for i in order[:5]]

    # First eigenvalue above each threshold
    nearest_zero = min(ZEROS, key=lambda z: abs(z - lsq))
    return {
        'lambda_sq': lsq,
        'L': float(mp.log(lsq)),
        'n_pp': len(vm_table(lsq)),
        'dims': dims,
        'lam_min': float(ev[0]),
        'top5_ev': top5,
        'nearest_zero': nearest_zero,
        'dist_to_zero': lsq - nearest_zero,
        't': round(dt, 2),
    }


def main():
    # Fine scan: step 0.05 from 4 to 22
    lsqs = list(np.arange(4.0, 22.05, 0.05))
    # Add exact zeros and nearby
    for z in ZEROS[:2]:
        for d in [-0.001, 0.0, 0.001]:
            lsqs.append(z + d)
    # Add exact integers
    for n in range(4, 23):
        lsqs.append(float(n))
    lsqs = sorted(set(round(x, 4) for x in lsqs))

    print(f"dim_evolution_scan | N={N} | G={G} | {len(lsqs)} points")
    print(f"{'lam_sq':>9}  {'n_pp':>4}  {'d1e-4':>5}  {'d1e-3':>5}  {'d1e-2':>5}  "
          f"{'lam_min':>12}  {'ev[1]':>12}  {'t':>5}")

    results = []
    for lsq in lsqs:
        try:
            r = analyze(lsq)
            d4 = r['dims']['1e-04']
            d3 = r['dims']['0.001']
            d2 = r['dims']['0.01']
            top = r['top5_ev']
            print(f"{lsq:>9.4f}  {r['n_pp']:>4}  {d4:>5}  {d3:>5}  {d2:>5}  "
                  f"{top[0]:>12.4e}  {top[1]:>12.4e}  {r['t']:>5.2f}s")
            sys.stdout.flush()
            results.append(r)
        except Exception as e:
            print(f"lsq={lsq:.4f} ERR: {e}")

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
