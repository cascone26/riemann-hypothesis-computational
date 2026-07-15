#!/usr/bin/env python3
"""
gamma3_fast.py — Fast γ₃ scan using sparse eigensolver for lam_min only.

At lsq=580-680 with N=600 (lower N but faster):
- N=600 gives N/lsq ≈ 1.0 (undersaturated but fast)
- We only need lam0 trend, not sc
- Use scipy.sparse.linalg.eigsh for just the smallest eigenvalue

Also try: N=400 at lsq=580-680 (even faster)
The γ₃ minimum at lsq=625 should show as lam0 approaching 0 from below.
"""
import time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 30000
OUTFILE = '/Users/scones/projects/riemann/results/gamma3_fast.out'

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

def lam_min(lsq, N):
    """Compute the smallest eigenvalue of τ using full diag on smaller matrix."""
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
        lkpp = np.array([float(mp.log(k)) for k, _ in VM])
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

    tau = (arch + B - P); tau = (tau + tau.T) / 2.0

    # Use sparse eigsh for just the smallest eigenvalue (much faster)
    tau_sp = scipy.sparse.csr_matrix(tau)
    try:
        ev, _ = scipy.sparse.linalg.eigsh(tau_sp, k=1, which='SA',
                                           tol=1e-8, maxiter=5000)
        return float(ev[0])
    except Exception:
        # Fall back to full eigh
        ev = scipy.linalg.eigvalsh(tau)
        return float(ev[0])

def main():
    # Scan lsq=560-700 with N=400 for speed
    lsq_vals = list(range(560, 705, 5))
    N = 400

    print(f"γ₃ fast scan: lsq=560-700, N={N} (sparse eigsh)")
    print(f"γ₃²=625.5 → minimum at lsq≈625")
    print(f"{'lsq':>5} {'N/lsq':>7} {'lam0':>14} {'dt':>5}")
    print("-" * 40)
    with open(OUTFILE, 'w') as f:
        f.write(f"gamma3_fast N={N} {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"γ₃²=625.5, sparse eigsh, N={N}\n\n")

    for lsq in lsq_vals:
        t0 = time.time()
        lam0 = lam_min(lsq, N)
        dt = time.time() - t0
        ratio = N / lsq
        label = " ← γ₃ n=0!" if abs(lsq - 625) <= 5 else ""
        line = f"{lsq:5d} {ratio:7.3f} {lam0:14.4e}  ({dt:.0f}s){label}"
        print(line, flush=True)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")

if __name__ == "__main__":
    main()
