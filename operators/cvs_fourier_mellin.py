#!/usr/bin/env python3
"""
cvs_fourier_mellin.py — Use connes-cvs ground state + our Fourier-Mellin formula.

Compares:
  (A) extract_zeros from the package (high-precision findroot)
  (B) our |F(z)| local-minimum method applied to the CvS eigenvector
  (C) our |F(z)| method applied to our tau's min_kinetic eigenvector

Goal: see if our B formula gives same zero locations as CvS, and check
which formula gives the deepest/most precise zero detection.
"""
import sys, json, os, time
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp
from scipy.signal import argrelmin

mp.mp.dps = 100
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000

try:
    from connes_cvs import build_galerkin_matrix, compute_ground_state, extract_zeros
    HAVE_CVS = True
except ImportError:
    HAVE_CVS = False
    print("WARNING: connes-cvs not found, skipping CvS eigenvector")

KNOWN_ZEROS = [
    14.134725141734693, 21.022039638771555, 25.010857580145688,
    30.424876125859513, 32.935061587739189, 37.586178158825671,
    40.918719012147495, 43.327073280914999, 48.005150881167159,
    49.773832477672302,
]

OUTFILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "../results/cvs_fourier_mellin.json")


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


def build_our_tau(lambda_sq, N):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N, N + 1, dtype=float)
    ns  = idx
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
        P = np.zeros((len(idx), len(idx)))
    tau = arch + B - P
    return (tau + tau.T) / 2.0, idx.astype(int)


def fourier_mellin(c, idx, L, z_grid):
    alpha = TWO_PI * idx[:, None] - L * z_grid[None, :]
    half  = 0.5 * alpha
    with np.errstate(divide='ignore', invalid='ignore'):
        kern = np.where(np.abs(half) < 1e-10, 1.0, np.sin(half) / half)
    return np.sqrt(L) * np.sum(c[:, None] * np.exp(1j * half) * kern, axis=0)


def min_kinetic_eigvec(ev, evec, idx, thresh=1e-2):
    mask = np.abs(ev) < thresh
    nn_evec = evec[:, mask]
    K = int(mask.sum())
    if K == 0:
        return evec[:, 0], 0
    n2 = (idx ** 2).astype(float)
    K_proj = (nn_evec * n2[:, None]).T @ nn_evec
    K_proj = (K_proj + K_proj.T) / 2
    kev, kevec = np.linalg.eigh(K_proj)
    coeff = kevec[:, 0]
    c = nn_evec @ coeff
    return c / np.linalg.norm(c), K


def analyze_zeros(c_vec, L, idx, z_grid, label):
    F = fourier_mellin(c_vec, idx, L, z_grid)
    F_abs = np.abs(F)
    Fmax = F_abs.max()
    lmin_idx = argrelmin(F_abs, order=10)[0]

    print(f"\n  [{label}] |F(gamma_k)|/max at known zeros:")
    checks = []
    for gk in KNOWN_ZEROS[:8]:
        iz = np.argmin(np.abs(z_grid - gk))
        fval = float(F_abs[iz])
        lo, hi = max(0, iz-200), min(len(F_abs)-1, iz+200)
        is_min = bool(fval == F_abs[lo:hi+1].min())
        mark = '*' if is_min else ' '
        print(f"   {mark}  gamma={gk:.4f}: |F|/max={fval/Fmax:.4e}")
        checks.append({'gamma': gk, 'F_abs': fval, 'F_abs_rel': fval/Fmax, 'is_local_min': is_min})

    # Candidate zero locations
    candidates = [float(z_grid[i]) for i in lmin_idx if F_abs[i] < 0.01 * Fmax]
    return {'label': label, 'zero_checks': checks,
            'n_candidates': len(candidates),
            'candidates_near_zeros': [c for c in candidates
                                      if min(abs(c - g) for g in KNOWN_ZEROS[:8]) < 0.5]}


def run(c, N_our=128, dps_cvs=80, T_cvs=400):
    L = float(np.log(c))
    print(f"\n{'='*60}")
    print(f"c={c}, N_our={N_our}, dps={dps_cvs}, T={T_cvs}")
    z_grid = np.linspace(5.0, 55.0, 100000)

    results = {'c': c, 'L': L, 'N_our': N_our}

    # --- (A) connes-cvs eigenvector ---
    if HAVE_CVS:
        mp.mp.dps = dps_cvs + 20
        t0 = time.time()
        Q = build_galerkin_matrix(c=c, N=N_our, T=T_cvs, dps=dps_cvs)
        lam_cvs, v_cvs = compute_ground_state(Q)
        t1 = time.time()
        print(f"  CvS: lam_min={mp.nstr(lam_cvs, 6)}  build+eig={t1-t0:.1f}s")

        # Extract the eigenvector as numpy array (it's in the EVEN sector, size N+1)
        # Full basis size 2N+1, so the eigenvector from compute_ground_state is (2N+1)
        c_cvs = np.array([float(v_cvs[i]) for i in range(len(v_cvs))])
        idx_cvs = np.arange(-N_our, N_our + 1, dtype=float)

        z_cvs = analyze_zeros(c_cvs, L, idx_cvs, z_grid, f"cvs_eigvec(c={c})")
        z_cvs['lam_min'] = str(lam_cvs)
        results['cvs'] = z_cvs

        # Also run package extract_zeros
        print(f"  Running extract_zeros...")
        t2 = time.time()
        try:
            zeros = extract_zeros(v_cvs, L=L, n_zeros=6, dps=dps_cvs)
            t3 = time.time()
            print(f"  extract_zeros: {t3-t2:.1f}s")
            pkg_zeros = []
            for z in zeros:
                err = z.get('error')
                gd = z.get('gamma_detected')
                print(f"    k={z['k']}: detected={mp.nstr(gd, 8) if gd else None}  err={mp.nstr(err, 6) if err else None}")
                pkg_zeros.append({'k': z['k'], 'gamma_detected': str(gd), 'error': str(err)})
            results['cvs_extract_zeros'] = pkg_zeros
        except Exception as e:
            print(f"  extract_zeros failed: {e}")
            results['cvs_extract_zeros'] = str(e)

    # --- (B) our tau, min_kinetic ---
    t4 = time.time()
    tau, idx_our = build_our_tau(c, N_our)
    ev, evec = scipy.linalg.eigh(tau)
    t5 = time.time()
    print(f"\n  Our tau: lam_min={ev[0]:.4e}  build+eig={t5-t4:.1f}s")
    K = int((np.abs(ev) < 1e-2).sum())
    print(f"  near-null K={K}")

    c_kin, K_kin = min_kinetic_eigvec(ev, evec, idx_our, thresh=1e-2)
    z_our = analyze_zeros(c_kin, L, idx_our, z_grid, f"our_min_kinetic(c={c})")
    z_our['lam_min'] = float(ev[0])
    z_our['K'] = K
    results['our_min_kinetic'] = z_our

    return results


def main():
    all_results = []

    # Compare at multiple cutoffs
    for c, N_our, dps_cvs, T_cvs in [
        (13, 100, 80, 400),
        (23, 100, 80, 400),
        (50, 100, 80, 400),
    ]:
        try:
            r = run(c, N_our=N_our, dps_cvs=dps_cvs, T_cvs=T_cvs)
            all_results.append(r)
        except Exception as e:
            print(f"ERROR c={c}: {e}")
            import traceback; traceback.print_exc()

    with open(OUTFILE, 'w') as f:
        json.dump(all_results, f, indent=2)
    print(f"\nSaved: {OUTFILE}")


if __name__ == "__main__":
    main()
