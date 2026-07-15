#!/usr/bin/env python3
"""
moved_zero_test.py — MOVE a Riemann zero from σ=0.5 to σ≠0.5, test sub_cancel.

The ghost_zero_test added an EXTRA zero (even on critical line, sc jumped 14x).
That was measuring "sensitivity to any A modification," not σ-specificity.

This test does the CORRECT thing: MOVE an existing zero from σ=0.5+iγ to σ₀+iγ.
  ΔA = [A contribution of σ₀+iγ] - [A contribution of 0.5+iγ]

If σ₀=0.5: ΔA=0, sc unchanged (control)
If σ₀≠0.5: ΔA≠0, question is whether sc grows MORE than for on-critical move

Also compare: moving zero vs adding zero.

First 3 Riemann zeros: γ₁=14.135, γ₂=21.022, γ₃=25.011

We move γ₁ from σ=0.5 to σ₀ ∈ {0.5, 0.51, 0.55, 0.6, 0.7}
and similarly for γ₂, γ₃.
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

OUTFILE = os.path.expanduser("~/riemann/results/moved_zero_test.json")

RIEMANN_ZEROS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
                 37.586178, 40.918719, 43.327073, 48.005151, 49.773832]


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


def build_components(lambda_sq, N_use):
    L  = float(mp.log(lambda_sq))
    sh = float(np.sinh(L / 4) ** 2)
    VM = vm_table(lambda_sq)
    idx = np.arange(-N_use, N_use + 1, dtype=float)
    ns  = idx
    R   = TWO_PI * N_use / L + 80.0
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
    return arch, B, P, E, grid, w, L, ns


def zero_delta_A(E, grid, w, sigma_new, gamma, sigma_old=0.5):
    """Contribution ΔA from MOVING zero γ from σ_old to σ_new.

    ΔΩ(xi) = Ω_contribution(sigma_new, gamma) - Ω_contribution(sigma_old, gamma)
    where Ω_contribution(σ, γ) = Re[1/(1/4+iξ/2 - (σ+iγ)) + 1/(σ+iγ)]
                                 + Re[1/(1/4+iξ/2 - (σ-iγ)) + 1/(σ-iγ)]
    """
    s0 = 0.25 + 0.5j * grid
    def contrib(sigma):
        rho = sigma + 1j * gamma
        rho_bar = sigma - 1j * gamma
        return (np.real(1.0 / (s0 - rho) + 1.0 / rho) +
                np.real(1.0 / (s0 - rho_bar) + 1.0 / rho_bar))
    delta_Om = contrib(sigma_new) - contrib(sigma_old)
    delta_A = ((E * (w * delta_Om)) @ E.conj().T).real / (2.0 * np.pi)
    return delta_A


def measure_sc(arch, B, P, E, grid, w, delta_A=None):
    AB = arch + B
    if delta_A is not None:
        AB = AB + delta_A
    tau = AB - P
    tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    null_mask = np.abs(ev) < THRESH
    dim = null_mask.sum()
    if dim == 0:
        return float('nan'), 0, float(ev[0])
    V_null = V[:, null_mask]
    AB_proj = V_null.T @ AB @ V_null
    P_proj  = V_null.T @ P  @ V_null
    sc = float(np.linalg.norm(AB_proj - P_proj, 'fro') /
               (np.linalg.norm(AB_proj, 'fro') + 1e-100))
    return sc, int(dim), float(ev[0])


def main():
    lsq = 50
    N_use = 250
    print(f"Moved zero test: lsq={lsq}, N={N_use}")
    t0 = time.time()
    arch, B, P, E, grid, w, L, ns = build_components(lsq, N_use)

    # Standard baseline
    sc_std, dim_std, lam_std = measure_sc(arch, B, P, E, grid, w)
    print(f"Standard: sc={sc_std:.4e}  dim={dim_std}  lam={lam_std:.3e}")

    results = {'standard': {'sc': sc_std, 'dim': dim_std, 'lam': lam_std}}

    # Test: move zero γ₁=14.135 from 0.5 to various σ
    sigmas = [0.50, 0.501, 0.505, 0.51, 0.55, 0.6, 0.7, 0.8]
    results['move_gamma1'] = {}
    print(f"\n=== Move γ₁=14.135 from σ=0.5 to σ_new ===")
    for sig in sigmas:
        dA = zero_delta_A(E, grid, w, sig, 14.134725)
        sc, dim, lam = measure_sc(arch, B, P, E, grid, w, dA)
        ratio = sc / sc_std
        print(f"  σ={sig:.4f}: sc={sc:.4e}  ratio={ratio:.2f}x  dim={dim}  lam={lam:.3e}")
        results['move_gamma1'][str(sig)] = {'sc': sc, 'dim': dim, 'lam': lam, 'ratio': ratio}

    # Test: move γ₂=21.022
    results['move_gamma2'] = {}
    print(f"\n=== Move γ₂=21.022 from σ=0.5 to σ_new ===")
    for sig in sigmas:
        dA = zero_delta_A(E, grid, w, sig, 21.022040)
        sc, dim, lam = measure_sc(arch, B, P, E, grid, w, dA)
        ratio = sc / sc_std
        print(f"  σ={sig:.4f}: sc={sc:.4e}  ratio={ratio:.2f}x  dim={dim}  lam={lam:.3e}")
        results['move_gamma2'][str(sig)] = {'sc': sc, 'dim': dim, 'lam': lam, 'ratio': ratio}

    # Test: move both γ₁ and γ₂ simultaneously
    results['move_gamma1_2'] = {}
    print(f"\n=== Move BOTH γ₁ AND γ₂ from σ=0.5 to σ_new ===")
    for sig in [0.5, 0.51, 0.55, 0.6, 0.7]:
        dA1 = zero_delta_A(E, grid, w, sig, 14.134725)
        dA2 = zero_delta_A(E, grid, w, sig, 21.022040)
        sc, dim, lam = measure_sc(arch, B, P, E, grid, w, dA1 + dA2)
        ratio = sc / sc_std
        print(f"  σ={sig:.4f}: sc={sc:.4e}  ratio={ratio:.2f}x  dim={dim}  lam={lam:.3e}")
        results['move_gamma1_2'][str(sig)] = {'sc': sc, 'dim': dim, 'lam': lam, 'ratio': ratio}

    # Test: move γ₁ at lsq=5 (fine resolution near γ₁)
    print(f"\n=== Scan Δσ near 0.5 (γ₁=14.135) ===")
    fine_sigmas = [0.5 + d for d in [0, 0.001, 0.002, 0.005, 0.010, 0.020, 0.050, 0.100]]
    results['sigma_scan_gamma1'] = {}
    for sig in fine_sigmas:
        dA = zero_delta_A(E, grid, w, sig, 14.134725)
        sc, dim, lam = measure_sc(arch, B, P, E, grid, w, dA)
        ratio = sc / sc_std
        print(f"  Δσ={sig-0.5:.4f}: sc={sc:.4e}  ratio={ratio:.2f}x  lam={lam:.3e}")
        results['sigma_scan_gamma1'][str(sig)] = {'sc': sc, 'dim': dim, 'lam': lam, 'ratio': ratio}

    out = {'lsq': lsq, 'N': N_use, 'L': round(L, 4), 'results': results,
           'total_time': round(time.time() - t0, 1)}
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
