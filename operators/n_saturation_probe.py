#!/usr/bin/env python3
"""
n_saturation_probe.py — Map the N-saturation transition at γ₁ resonance.

At lsq=200 (γ₁ main resonance):
  N=400 (2*lsq): ratio=2.891 SPIKE
  N=600 (3*lsq): ratio=0.579 sub
  N=1000 (5*lsq, Olympus): ratio=0.033 DEEP MIN

This documents the exact N at which the resonance minimum emerges.
Also probe lsq=250 (γ₁ n=1) and lsq=160 (γ₁ n=-1 spike, should stay spike at all N).

Prediction:
  - Resonance minima NEED N ≥ 5*lsq to be fully visible
  - Precursor spikes are ROBUST to N (stay spikes at any reasonable N)
  - N-saturation threshold: N_sat ≈ 4-5*lsq for resonance points
"""
import time, json
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 50000
THRESH = 0.02
OUTFILE = '/Users/scones/projects/riemann/results/n_saturation_probe.out'

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

def measure(lsq, N):
    L  = float(mp.log(lsq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lsq)
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
    AB = arch + B
    tau = (AB - P); tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    lam0 = float(ev[0])
    null_mask = np.abs(ev) < THRESH
    dim = int(null_mask.sum())
    if dim > 0:
        Vn = V[:, null_mask]
        ABp = Vn.T @ AB @ Vn
        Pp  = Vn.T @ P  @ Vn
        sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
    else:
        sc = float('nan')
    return sc, lam0, dim

def main():
    # Probe lsq=200 at multiple N values: map the N-saturation transition
    # Also lsq=160 (spike, should be robust) and lsq=80 (control)
    probes = []

    # lsq=200: N = 1*lsq to 6*lsq (max 1200 on Mac)
    for mult in [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5]:
        N = int(mult * 200)
        if N > 1200: N = 1200
        probes.append((200, N, f'γ₁n=0 N={mult}*lsq'))

    # lsq=160: N = 2*lsq to 5*lsq (spike should be stable)
    for mult in [2, 3, 4, 5]:
        N = int(mult * 160)
        if N > 1200: N = 1200
        probes.append((160, N, f'γ₁n=-1 N={mult}*lsq'))

    # lsq=250: N = 2*lsq to 5*lsq
    for mult in [2, 3, 4]:
        N = int(mult * 250)
        if N > 1200: N = 1200
        probes.append((250, N, f'γ₁n=1 N={mult}*lsq'))

    print(f"N-saturation probe: tracking sc vs N for lsq=200, 160, 250")
    print(f"{'lsq':>5} {'N':>5} {'label':25} {'sc':>12} {'pred':>12} {'ratio':>7} {'dt':>5}")
    print("-" * 80)
    with open(OUTFILE, 'w') as f:
        f.write(f"N-saturation probe {time.strftime('%Y-%m-%d %H:%M:%S')}\n")

    results = []
    t0 = time.time()
    for lsq, N, label in probes:
        L = float(mp.log(lsq))
        t1 = time.time()
        sc, lam0, dim = measure(lsq, N)
        dt = time.time() - t1
        sc_pred = 3.073e-3 / L ** 1.9151
        ratio = sc / sc_pred if sc_pred > 0 and not np.isnan(sc) else float('nan')
        line = f"{lsq:5d} {N:5d} {label:25}: sc={sc:.4e}  pred={sc_pred:.4e}  ratio={ratio:7.3f}  ({dt:.0f}s)"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")
        results.append({'lsq': lsq, 'N': N, 'label': label, 'sc': sc, 'ratio': round(ratio,4) if not np.isnan(ratio) else None})

    with open(OUTFILE.replace('.out', '.json'), 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nTotal: {time.time()-t0:.1f}s  →  {OUTFILE}")

if __name__ == "__main__":
    main()
