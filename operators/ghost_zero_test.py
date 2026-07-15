#!/usr/bin/env python3
"""
ghost_zero_test.py — Add a synthetic off-critical zero to A+B and check sub_cancel.

THE QUESTION: if τ = A+B-P and we ADD a fake zero at ρ₀=σ₀+iγ₀ (with σ₀≠0.5)
to the A matrix (via its contribution to the Weil digamma function), does sub_cancel
increase (CANNOT be canceled by P)? If YES, this confirms the backward direction:
sub_cancel → 0 ⟺ all zeros on critical line.

The modification: add contribution of zero at ρ₀ to the Ω integrand:
  ΔΩ(ξ) = Re[1/(1/4 + iξ/2 - ρ₀) + 1/ρ₀] + Re[1/(1/4 + iξ/2 - ρ₀^* + 1/ρ₀^*)]

This represents a "deformed" ζ function with an EXTRA zero at ρ₀ ≠ 1/2 + iγ.
The P matrix remains the REAL prime weights (unchanged).

Tests:
1. Standard τ (all real zeros on σ=1/2): sc_standard
2. Modified τ with extra zero at ρ₀ = 0.6 + i×14.13 (near first real zero, off-critical)
3. Modified τ with extra zero at ρ₀ = 0.51 + i×14.13 (tiny deviation)
4. Modified τ with extra zero at ρ₀ = 0.7 + i×100 (far off-critical, large imaginary part)
5. Control: extra zero at ρ₀ = 0.5 + i×14.13 (on critical line — should NOT increase sc)
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

OUTFILE = os.path.expanduser("~/riemann/results/ghost_zero_test.json")


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
    """Build A, B, P components separately for ghost-zero modification."""
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


def ghost_delta_A(E, grid, w, L, rho_0):
    """Compute the A-matrix modification from an extra zero at rho_0.

    The Weil digamma sum includes Re[1/(s - rho) + 1/rho] for each zero rho.
    Adding a zero at rho_0 adds:
      ΔΩ(xi) = Re[1/(1/4 + i*xi/2 - rho_0)] + Re[1/rho_0]
             + Re[1/(1/4 + i*xi/2 - conj(rho_0))] + Re[1/conj(rho_0)]

    Then ΔA_{mn} = ∫ ΔΩ(xi) * E_m(xi) * E_n^*(xi) dxi / (2*pi)

    Note: rho_0 and its complex conjugate both contribute (from xi<->xi symmetry
    of the Riemann xi function, zeros come in conjugate pairs).
    """
    s0 = 0.25 + 0.5j * grid  # 1/4 + i*xi/2 on grid
    # Contribution from rho_0
    delta_Om = np.real(1.0 / (s0 - rho_0) + 1.0 / rho_0)
    # Contribution from conj(rho_0)
    rho_bar = np.conj(rho_0)
    delta_Om += np.real(1.0 / (s0 - rho_bar) + 1.0 / rho_bar)
    delta_A = ((E * (w * delta_Om)) @ E.conj().T).real / (2.0 * np.pi)
    return delta_A


def measure_sc_modified(arch, B, P, E, grid, w, rho_0=None, scale=1.0):
    """Measure sub_cancel for τ = (A + ΔA) + B - P."""
    AB = arch + B
    if rho_0 is not None:
        dA = ghost_delta_A(E, grid, w, None, rho_0) * scale
        AB = AB + dA
    tau = AB - P
    tau = (tau + tau.T) / 2.0
    ev, V = scipy.linalg.eigh(tau)
    null_mask = np.abs(ev) < THRESH
    dim = null_mask.sum()
    if dim == 0:
        return float('nan'), 0, float(ev[0]), float(np.nan)
    V_null = V[:, null_mask]
    AB_proj = V_null.T @ AB @ V_null
    P_proj  = V_null.T @ P  @ V_null
    sc = float(np.linalg.norm(AB_proj - P_proj, 'fro') / (np.linalg.norm(AB_proj, 'fro') + 1e-100))
    dA_proj = float(np.linalg.norm(V_null.T @ (AB - arch - B) @ V_null, 'fro')) if rho_0 is not None else 0.0
    return sc, int(dim), float(ev[0]), dA_proj


def main():
    lsq = 50
    N_use = 250  # 5 * lsq

    print(f"Ghost zero test: lsq={lsq}, N={N_use}")
    print(f"Building components...", flush=True)
    t0 = time.time()
    arch, B, P, E, grid, w, L, ns = build_components(lsq, N_use)

    # Test cases: (label, rho_0)
    cases = [
        ("standard (no ghost)", None),
        ("ghost σ=0.5+0i γ=14.13 (on-critical)", 0.5 + 14.13j),   # on critical → should NOT change much
        ("ghost σ=0.501 γ=14.13 (barely off)", 0.501 + 14.13j),
        ("ghost σ=0.51 γ=14.13 (slightly off)", 0.51 + 14.13j),
        ("ghost σ=0.55 γ=14.13 (moderately off)", 0.55 + 14.13j),
        ("ghost σ=0.6 γ=14.13 (off-critical, near first zero)", 0.6 + 14.13j),
        ("ghost σ=0.7 γ=14.13 (far off-critical)", 0.7 + 14.13j),
        ("ghost σ=0.6 γ=100 (off-critical, large γ)", 0.6 + 100j),
        ("ghost σ=0.6 γ=50 (off-critical, med γ)", 0.6 + 50j),
    ]

    results = []
    for label, rho_0 in cases:
        t1 = time.time()
        sc, dim, lam_min, dA_proj = measure_sc_modified(arch, B, P, E, grid, w, rho_0)
        dt = time.time() - t1
        print(f"  [{label}]: sc={sc:.4e}  dim={dim}  lam={lam_min:.3e}  ({dt:.1f}s)", flush=True)
        results.append({
            'label': label,
            'rho_0': str(rho_0) if rho_0 else None,
            'sc': sc,
            'dim': dim,
            'lam_min': lam_min,
            'dA_proj': dA_proj,
            'dt': round(dt, 2),
        })

    # Test scaling with ghost amplitude
    print(f"\n=== Ghost amplitude scan (σ=0.6, γ=14.13) ===")
    rho_ghost = 0.6 + 14.13j
    for scale in [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]:
        t1 = time.time()
        sc, dim, lam_min, _ = measure_sc_modified(arch, B, P, E, grid, w, rho_ghost, scale)
        dt = time.time() - t1
        print(f"  scale={scale:5.2f}: sc={sc:.4e}  lam={lam_min:.3e}  ({dt:.1f}s)", flush=True)

    out = {
        'lsq': lsq, 'N': N_use, 'L': round(L, 4),
        'results': results,
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
