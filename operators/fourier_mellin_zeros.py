#!/usr/bin/env python3
"""
fourier_mellin_zeros.py — Extract Riemann zero approximations from tau ground state
via Fourier-Mellin transform (Groskin 2605.20224 / Suzuki 2606.09096).

Method:
  Build tau at (lambda_sq, N), extract a good "ground state" eigenvector c,
  then compute F(z) = sqrt(L) * integral_0^L xi_N(t) e^{-izt} dt
  and find zeros/minima of |F(z)| on real z-axis.

Three eigenvector selection strategies:
  1. min_ev:       eigenvector with smallest eigenvalue (most negative)
  2. min_kinetic:  linear combination in near-null subspace with min Sigma_n n^2 |c_n|^2
  3. dc_dominant:  near-null eigenvector with maximum |c_0|^2 power

Known Riemann zeros: gamma_1=14.135, gamma_2=21.022, gamma_3=25.011, ...
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp
from scipy.signal import argrelmin

mp.mp.dps = 50
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000

KNOWN_ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302,
]

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/fourier_mellin_zeros.json")


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
    return (tau + tau.T) / 2.0, idx.astype(int)


def fourier_mellin(c, idx, L, z_grid):
    """
    F(z) = sqrt(L) * Integral_0^L xi_N(t) e^{-izt} dt
         = sqrt(L) * sum_n c_n * e^{i*alpha_n/2} * sin(alpha_n/2) / (alpha_n/2)
    where alpha_n = 2*pi*n - z*L.
    Returns complex array shape (len(z_grid),).
    """
    # alpha_n(z) = 2*pi*n - z*L, shape (2N+1, Nz)
    alpha = TWO_PI * idx[:, None] - L * z_grid[None, :]
    half  = 0.5 * alpha
    with np.errstate(divide='ignore', invalid='ignore'):
        kern = np.where(np.abs(half) < 1e-10, 1.0, np.sin(half) / half)
    F = np.sqrt(L) * np.sum(c[:, None] * np.exp(1j * half) * kern, axis=0)
    return F


def select_eigenvectors(ev, evec, idx, thresh=1e-2):
    """
    Return three candidate eigenvectors from the near-null subspace:
      - 'min_ev':       index with most negative eigenvalue
      - 'min_kinetic':  linear combination minimizing Sigma n^2 |c_n|^2
      - 'dc_dominant':  near-null eigenvec with max |c_0|^2 (DC power)
    """
    # near-null mask (|ev| < thresh)
    mask = np.abs(ev) < thresh
    nn_ev   = ev[mask]
    nn_evec = evec[:, mask]   # shape (2N+1, K)
    K = int(mask.sum())

    # 1. min_ev: eigenvector with most negative eigenvalue (classical ground state)
    c_min = evec[:, 0]  # already sorted by eigh

    results = {'min_ev': (c_min, float(ev[0]), 'min_ev')}

    if K == 0:
        return results, 0

    # 2. min_kinetic: min Sigma n^2 |c|^2 in the near-null subspace
    #    = min eigenvalue of the projected kinetic matrix K_proj = V^T diag(n^2) V
    n2 = (idx ** 2).astype(float)
    K_proj = (nn_evec * n2[:, None]).T @ nn_evec  # (K, K)
    K_proj = (K_proj + K_proj.T) / 2
    kev, kevec = np.linalg.eigh(K_proj)
    coeff_kin = kevec[:, 0]  # min kinetic energy direction in near-null subspace
    c_kin = nn_evec @ coeff_kin
    c_kin /= np.linalg.norm(c_kin)
    results['min_kinetic'] = (c_kin, float(kev[0]), 'min_kinetic')

    # 3. dc_dominant: near-null eigvec with max |c_{n=0}|^2
    dc_idx = np.argmax(idx == 0)
    dc_power = nn_evec[dc_idx, :] ** 2
    best_dc = np.argmax(dc_power)
    c_dc = nn_evec[:, best_dc]
    results['dc_dominant'] = (c_dc, float(nn_ev[best_dc]), 'dc_dominant')

    return results, K


def analyze_fm(c, label, L, idx, z_grid, F_abs_grid):
    """Analyze Fourier-Mellin on given z_grid and return zero candidates."""
    local_min_idx = argrelmin(F_abs_grid, order=10)[0]
    global_max = F_abs_grid.max()
    candidates = []
    for mi in local_min_idx:
        zc = z_grid[mi]
        fval = F_abs_grid[mi]
        if fval < 0.05 * global_max:
            nearest = min(KNOWN_ZEROS, key=lambda z: abs(z - zc))
            dist = zc - nearest
            candidates.append({
                'z': float(zc), 'F_abs': float(fval),
                'F_abs_rel': float(fval / global_max),
                'nearest_zero': float(nearest), 'dist': float(dist),
            })
    return candidates


def run_case(lambda_sq, N, z_min=5.0, z_max=55.0, Nz=100000):
    L = float(np.log(lambda_sq))
    print(f"\n{'='*60}")
    print(f"tau build: lambda_sq={lambda_sq}, N={N}, L={L:.4f}")
    t0 = time.time()
    tau, idx = build_tau(lambda_sq, N)
    t1 = time.time()
    ev, evec = scipy.linalg.eigh(tau)
    t2 = time.time()
    print(f"  build={t1-t0:.1f}s  eigh={t2-t1:.1f}s")
    print(f"  ev[0]={ev[0]:.4e}  ev[1]={ev[1]:.4e}  ev[-1]={ev[-1]:.4e}")

    # Count near-null modes
    for thresh in [1e-4, 1e-3, 1e-2]:
        print(f"  dim(|ev|<{thresh:.0e})={int((np.abs(ev)<thresh).sum())}", end='')
    print()

    # z grid for Fourier-Mellin
    z_grid = np.linspace(z_min, z_max, Nz)

    # Get three eigenvector candidates
    ev_cands, K = select_eigenvectors(ev, evec, idx, thresh=1e-2)
    print(f"  near-null K={K}")

    # Evaluate |F(z)| for each candidate
    case_results = {
        'lambda_sq': lambda_sq, 'N': N, 'L': L,
        'lam_min': float(ev[0]), 'K': K,
        'strategies': {}
    }

    for strat_name, (c_vec, ev_val, _) in ev_cands.items():
        F = fourier_mellin(c_vec, idx, L, z_grid)
        F_abs = np.abs(F)
        candidates = analyze_fm(c_vec, strat_name, L, idx, z_grid, F_abs)

        print(f"\n  Strategy '{strat_name}' (ev={ev_val:.4e}):")
        if not candidates:
            print("    No zero candidates found.")
        else:
            print(f"    {'z_cand':>10}  {'|F|/max':>10}  {'nearest_zero':>12}  {'dist':>10}")
            for cand in candidates[:15]:
                print(f"    {cand['z']:>10.5f}  {cand['F_abs_rel']:>10.4e}"
                      f"  {cand['nearest_zero']:>12.5f}  {cand['dist']:>+10.5f}")

        # Also evaluate at known zeros
        zero_check = []
        for gk in KNOWN_ZEROS[:8]:
            iz = np.argmin(np.abs(z_grid - gk))
            fval = float(F_abs[iz])
            lo, hi = max(0, iz-200), min(len(F_abs)-1, iz+200)
            is_min = bool(fval == F_abs[lo:hi+1].min())
            zero_check.append({'gamma': gk, 'F_abs': fval,
                                'F_abs_rel': fval / float(F_abs.max()),
                                'is_local_min': is_min})

        print(f"    |F| at known zeros:")
        for zc in zero_check:
            mark = '*' if zc['is_local_min'] else ' '
            print(f"   {mark}  gamma={zc['gamma']:.4f}: |F|/max={zc['F_abs_rel']:.4e}")

        case_results['strategies'][strat_name] = {
            'ev_val': float(ev_val),
            'candidates': candidates,
            'zero_check': zero_check,
        }

    return case_results


def main():
    all_results = []

    # Case 1: lambda_sq=50, N=192 (feasible on Mac)
    r = run_case(50.0, 192, z_min=5.0, z_max=55.0)
    all_results.append(r)
    sys.stdout.flush()

    # Case 2: lambda_sq=100, N=256
    r = run_case(100.0, 256, z_min=5.0, z_max=55.0)
    all_results.append(r)
    sys.stdout.flush()

    # Case 3: lambda_sq=14, N=192 (near gamma_1, tight scan)
    r = run_case(14.134725, 192, z_min=5.0, z_max=30.0)
    all_results.append(r)
    sys.stdout.flush()

    with open(OUTFILE, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
