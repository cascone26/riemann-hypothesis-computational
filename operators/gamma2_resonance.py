#!/usr/bin/env python3
"""
gamma2_resonance.py — Test the gamma-resonance theory at lsq ≈ γ₂² ≈ 442.

From §52: sc has deep minima at lsq ≈ γ_k² (when λ=√lsq ≈ γ_k).
Confirmed: γ₁ resonance at lsq=200 (ratio 0.033, 30× below power law).

Prediction: γ₂ resonance at lsq=442 (γ₂=21.022, γ₂²=441.9).
Expected pattern: spike at lsq=400-430, deep minimum at lsq≈442.

Scan: lsq=380,400,410,420,425,430,435,440,442,445,450,460,480,500
N=min(5*lsq, 2000) — capped for memory.
"""
import os, json, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
THRESH = 0.02
OUTFILE = os.path.expanduser("~/riemann/results/gamma2_resonance.json")


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
    exp_dim = 2 * (lsq - 2) + 1
    if dim > 0:
        Vn = V[:, null_mask]
        ABp = Vn.T @ AB @ Vn
        Pp  = Vn.T @ P  @ Vn
        sc = float(np.linalg.norm(ABp - Pp, 'fro') / (np.linalg.norm(ABp, 'fro') + 1e-100))
    else:
        sc = float('nan')
    return sc, lam0, dim, exp_dim


def main():
    # Dense scan around γ₂²=441.9
    # Also include wider context to see spike before and recovery after
    lsq_list = [380, 400, 410, 420, 425, 430, 435, 440, 441, 442, 443,
                445, 450, 460, 480, 500]

    def N_for(lsq):
        return min(5 * lsq, 2000)  # cap at N=2000 (matrix 4001x4001, ~128MB)

    results = []
    t0 = time.time()
    print("Gamma-2 resonance test: scanning lsq=380..500 for deep minimum at lsq≈442")
    print(f"γ₁=14.135 → lsq_resonance=200 (confirmed)")
    print(f"γ₂=21.022 → lsq_resonance=442 (PREDICTION)")
    print(f"{'lsq':>5} {'L':>6} {'N':>5} {'sc':>12} {'pred':>12} {'ratio':>7} {'lam0':>12} {'dim':>10} {'dt':>6}")
    for lsq in lsq_list:
        N = N_for(lsq)
        L = float(mp.log(lsq))
        t1 = time.time()
        sc, lam0, dim, exp_dim = measure(lsq, N)
        dt = time.time() - t1
        sc_pred = 3.073e-3 / L ** 1.9151
        ratio = sc / sc_pred if sc_pred > 0 and not np.isnan(sc) else float('nan')
        print(f"{lsq:5d} {L:6.4f} {N:5d}: sc={sc:.4e}  pred={sc_pred:.4e}  ratio={ratio:.3f}  lam0={lam0:.2e}  dim={dim}/{exp_dim}  ({dt:.1f}s)", flush=True)
        results.append({'lsq': lsq, 'L': round(L,5), 'N': N, 'sc': sc,
                        'sc_pred': round(sc_pred, 8), 'ratio': round(ratio, 5) if not np.isnan(ratio) else None,
                        'lam0': lam0, 'dim': dim, 'exp_dim': exp_dim, 'dt': round(dt,1)})

    out = {'lsq_list': lsq_list, 'results': results, 'total_time': round(time.time()-t0,1),
           'gamma1': 14.134725, 'gamma2': 21.022040,
           'resonance_lsq1': 199.8, 'resonance_lsq2': 441.9}
    with open(OUTFILE,'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")

if __name__ == "__main__":
    main()
