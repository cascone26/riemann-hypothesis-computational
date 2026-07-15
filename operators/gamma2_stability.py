#!/usr/bin/env python3
"""
gamma2_stability.py — Is the dim_wide=39 at lambda^2=gamma_2 stable as N increases?

At gamma_2=21.022: dim_wide jumped from 38 to 39 right at the zero.
This checks: is this a genuine spectral feature or N=64 artifact?

Also scans very finely around gamma_1 (14.135) and gamma_2 (21.022) at larger N.
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
OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/gamma2_stability.json")

GAMMA1 = 14.134725141734693
GAMMA2 = 21.022039638771555


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


def analyze(lsq, N):
    t0  = time.time()
    tau = build_tau(lsq, N)
    ev  = scipy.linalg.eigvalsh(tau)
    dt  = time.time() - t0
    dims = {'1e-04': int((np.abs(ev) < 1e-4).sum()),
            '0.001': int((np.abs(ev) < 1e-3).sum()),
            '0.01':  int((np.abs(ev) < 1e-2).sum())}
    order = np.argsort(np.abs(ev))
    top8  = [float(ev[i]) for i in order[:8]]
    return {
        'lambda_sq': lsq, 'N': N,
        'dims': dims,
        'lam_min': float(ev[0]),
        'top8_ev': top8,
        't': round(dt, 2),
    }


def main():
    results = []

    # 1. N-stability at gamma_1 and gamma_2
    print("=== N-stability at gamma_1 and gamma_2 ===")
    for lsq, label in [(GAMMA1, 'gamma_1'), (GAMMA2, 'gamma_2')]:
        print(f"\nlambda_sq = {lsq:.6f} ({label})")
        print(f"{'N':>5}  {'d1e-4':>6}  {'d1e-3':>6}  {'d1e-2':>6}  {'lam_min':>12}  "
              f"{'ev[0]':>12}  {'ev[1]':>12}  {'t':>6}")
        for N in [32, 48, 64, 96, 128, 192]:
            r = analyze(lsq, N)
            top = r['top8_ev']
            print(f"{N:>5}  {r['dims']['1e-04']:>6}  {r['dims']['0.001']:>6}  "
                  f"{r['dims']['0.01']:>6}  {r['lam_min']:>12.4e}  "
                  f"{top[0]:>12.4e}  {top[1]:>12.4e}  {r['t']:>6.2f}s")
            sys.stdout.flush()
            r['label'] = label
            results.append(r)

    # 2. Fine lambda_sq scan around gamma_2 at large N
    print(f"\n=== Fine scan around gamma_2={GAMMA2:.4f} at N=128 ===")
    N2 = 128
    lsqs2 = sorted([
        GAMMA2 - 0.1, GAMMA2 - 0.05, GAMMA2 - 0.02, GAMMA2 - 0.01,
        GAMMA2 - 0.005, GAMMA2 - 0.002, GAMMA2 - 0.001,
        GAMMA2, GAMMA2 + 0.001, GAMMA2 + 0.002, GAMMA2 + 0.005,
        GAMMA2 + 0.01, GAMMA2 + 0.02, GAMMA2 + 0.05, GAMMA2 + 0.1,
        21.0, 21.5, 22.0,
    ])
    print(f"{'lam_sq':>10}  {'dist':>9}  {'d1e-3':>6}  {'d1e-2':>6}  "
          f"{'lam_min':>12}  {'ev[0]':>12}  {'t':>6}")
    for lsq in lsqs2:
        r = analyze(lsq, N2)
        top = r['top8_ev']
        dist = lsq - GAMMA2
        print(f"{lsq:>10.6f}  {dist:>+9.4f}  {r['dims']['0.001']:>6}  "
              f"{r['dims']['0.01']:>6}  {r['lam_min']:>12.4e}  "
              f"{top[0]:>12.4e}  {r['t']:>6.2f}s")
        sys.stdout.flush()
        r['label'] = f'gamma2_fine_N{N2}'
        results.append(r)

    with open(OUTFILE, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
