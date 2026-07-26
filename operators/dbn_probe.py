#!/usr/bin/env python3
"""dbn_probe.py -- de Bruijn-Newman probe on the finite Weil form tau.

STRUCTURAL QUESTION (Sam, room discussion 2026-07-24):
  RH <=> Lambda = 0 (de Bruijn-Newman), with Lambda >= 0 PROVEN (Rodgers-Tao 2020).
  Lambda = 0 means tau sits EXACTLY on the PSD boundary: strict positive margin for
  the "stable" flow direction t>0, marginal (near-null kernel) at t=0.
  CLAIM TO TEST: our near-null tower IS the DBN-marginal kernel at t=0.
  FALSIFIABLE PREDICTION: there is a sign of t for which lambda_min(t) lifts
  monotonically off ~0 and the near-null DIMENSION collapses to 0 for |t|>0.
  If neither sign does this (margin non-monotone, or kernel persists), the near-null
  is a robust interior degeneracy, NOT DBN-marginal -> cleanly SEPARATES them.

DEFORMATION (spectral r-domain, consistent across arch+prime):
  e_hat_n(r) -> exp(-t r^2 / 2) e_hat_n(r)  =>  Gram gets exp(-t r^2).
  arch(t)  = (1/2pi) int Omega(r) exp(-t r^2) e_n(r) e_m*(r) dr
  prime(t) = sum_k Lam(k) k^-1/2 * q_t(n,m,log k),
             q_t(n,m,y) = Re[ (1/2pi) int exp(-t r^2) e_n(r) e_m*(r) e^{-i r y} dr ]
  W02 boundary/pole term held fixed (r-independent closed form) -- flagged, v1.

SANITY GATE (must pass or the object is wrong): grid-built prime at t=0 must
reproduce the closed-form prime_mat (connes_decisive_v3) to ~1e-6. Only then are
the t!=0 numbers trustworthy. Reported explicitly.
"""
import numpy as np, mpmath as mp, json
from scipy.special import digamma as _digamma

mp.mp.dps = 30
LAMBDA_SQ = 14
L = float(2 * mp.log(mp.sqrt(LAMBDA_SQ)))
TWO_PI = 2.0 * np.pi
LOGPI = float(np.log(np.pi))


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


def omega_vec(r):
    return np.real(_digamma(0.25 + 0.5j * r)) - LOGPI


def e_hat_vec(r, n):
    s = TWO_PI * n / L
    d = s - r
    out = np.empty_like(d, dtype=complex)
    small = np.abs(d) < 1e-12
    out[small] = np.sqrt(L) + 0j
    dd = d[~small]
    out[~small] = (1.0 / np.sqrt(L)) * np.exp(0.5j * dd * L) * 2.0 * np.sin(0.5 * dd * L) / dd
    return out


def q_kernel_closed(n, m, y):
    if n == m:
        return 2.0 * (1.0 - abs(y) / L) * np.cos(TWO_PI * n * y / L)
    return (np.sin(TWO_PI * m * y / L) - np.sin(TWO_PI * n * y / L)) / (np.pi * (n - m))


def w02_mat(idx):
    D = len(idx)
    W = np.zeros((D, D))
    sh = np.sinh(L / 4) ** 2
    for ii, n in enumerate(idx):
        for jj, m in enumerate(idx):
            num = (L ** 2 - 16 * np.pi ** 2 * m * n)
            den = (L ** 2 + 16 * np.pi ** 2 * m ** 2) * (L ** 2 + 16 * np.pi ** 2 * n ** 2)
            W[ii, jj] = 32 * L * sh * num / den
    return W


def build_tau_dbn(N, R, G, t):
    """tau(t) with DBN spectral deformation exp(-t r^2) on the Gram. Returns tau, prime_check."""
    idx = list(range(-N, N + 1))
    D = len(idx)
    grid = np.linspace(-R, R, G)
    h = grid[1] - grid[0]
    w = np.full(G, h); w[0] = w[-1] = 0.5 * h
    Om = omega_vec(grid)
    heat = np.exp(-t * grid ** 2)               # DBN deformation factor on the Gram
    E = np.empty((D, G), dtype=complex)
    for k, n in enumerate(idx):
        E[k] = e_hat_vec(grid, n)
    # arch(t)
    A = E * (w * Om * heat)
    arch = (A @ E.conj().T).real / (2.0 * np.pi)
    # prime(t): grid-built q_t for every k, deformed by heat, and a t=0 closed-form check
    Pt = np.zeros((D, D))
    prime_err = 0.0
    Wk = E * (w * heat)                          # (D,G) weighted rows
    for k, Lam in VM:
        y = np.log(k)
        phase = np.exp(-1j * grid * y)           # (G,)
        # q_t(n,m,y) = Re[(1/2pi) sum_grid w*heat * e_n * conj(e_m) * e^{-i r y}]
        qk = (Wk * phase) @ E.conj().T
        qk = qk.real / (2.0 * np.pi)
        Pt += Lam * k ** (-0.5) * qk
        if t == 0.0:                             # sanity: compare to closed form
            for ii, n in enumerate(idx):
                for jj, m in enumerate(idx):
                    prime_err = max(prime_err, abs(qk[ii, jj] - q_kernel_closed(n, m, y)))
    Pt = (Pt + Pt.T) / 2
    tau = arch + w02_mat(idx) - Pt
    tau = (tau + tau.T) / 2
    return tau, prime_err


def near_null_dim(ev, tol_rel=1e-8):
    mx = max(abs(ev[0]), abs(ev[-1]))
    return int(np.sum(np.abs(ev) < tol_rel * mx))


if __name__ == "__main__":
    N, R, G = 24, 60.0, 24000
    print(f"DBN probe: N={N} (dim {2*N+1}), grid R={R} G={G}, LAMBDA_SQ={LAMBDA_SQ}, L={L:.5f}")
    print("Deformation exp(-t r^2) on Gram (both signs). W02 boundary held fixed (v1).\n")

    # t=0 sanity gate first
    tau0, perr = build_tau_dbn(N, R, G, 0.0)
    ev0 = np.linalg.eigvalsh(tau0)
    print(f"[SANITY GATE] grid-prime vs closed-form max|diff| at t=0: {perr:.3e}  "
          f"({'PASS' if perr < 1e-6 else 'FAIL -- object wrong, do not trust t!=0'})")
    print(f"[t=0] lambda_min={ev0[0]:.6e}  near-null dim={near_null_dim(ev0)}  "
          f"smallest 6 |eig|={[f'{x:.2e}' for x in sorted(np.abs(ev0))[:6]]}\n")

    ts = [-1e-1, -1e-2, -1e-3, -1e-4, 1e-4, 1e-3, 1e-2, 1e-1]
    print(f"{'t':>10} {'lambda_min':>14} {'nn_dim':>7} {'smallest 4 |eig|'}")
    print(f"{0.0:>10.0e} {ev0[0]:>14.6e} {near_null_dim(ev0):>7d}  {[f'{x:.2e}' for x in sorted(np.abs(ev0))[:4]]}")
    rows = [(0.0, float(ev0[0]), near_null_dim(ev0))]
    for t in ts:
        tau, _ = build_tau_dbn(N, R, G, t)
        ev = np.linalg.eigvalsh(tau)
        nn = near_null_dim(ev)
        print(f"{t:>10.0e} {ev[0]:>14.6e} {nn:>7d}  {[f'{x:.2e}' for x in sorted(np.abs(ev))[:4]]}")
        rows.append((t, float(ev[0]), nn))

    rows.sort()
    print("\nREAD: if one sign of t lifts lambda_min monotonically off ~0 AND collapses nn_dim -> 0,")
    print("the near-null is the DBN-marginal kernel (RH<=>Lambda=0 knife-edge). If neither: SEPARATED.")
    json.dump({"N": N, "R": R, "G": G, "sanity_prime_err": perr,
               "rows": rows}, open("results/dbn_probe_out.json", "w"), indent=2)
    print("-> results/dbn_probe_out.json")
