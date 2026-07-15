#!/usr/bin/env python3
"""
analytic_slope_correct.py — Correct the H-F discrepancy.

The analytic_slope_verify.py computed ⟨ψ₀|dA/dσ|ψ₀⟩ for the GROUND STATE
and got ~-7e-3, while numeric slope is ~-1.076.

INSIGHT: The lam_min response is NOT the ground state H-F slope.
It is: lam_g - lam0 ≈ min_k ⟨ψ_k|ΔA|ψ_k⟩ where min is over all null modes k.

Or equivalently: min eigenvalue of (V_null^T × ΔA × V_null) × ds.

This script:
1. Builds tau, finds null modes V_null
2. Computes ΔA_fd (finite difference, actual moved-zero delta)
3. Projects: ΔA_null = V_null^T × ΔA × V_null
4. Computes min_null_eig = min eigenvalue of ΔA_null
5. Predicted lam_shift = min_null_eig × ds (should match lam_g - lam0)
6. Also computes: diagonal elements ⟨ψ_k|ΔA|ψ_k⟩ and their minimum (first-order H-F)
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
GAMMA1 = 14.134725
DELTA_SIGMA = 0.001
OUTFILE = os.path.expanduser("~/riemann/results/analytic_slope_correct.json")


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


def build_and_analyze(lsq, N):
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

    # Null subspace
    null_mask = np.abs(ev) < THRESH
    dim = int(null_mask.sum())
    Vn = V[:, null_mask]  # (D, dim)

    # Numerical ΔA from moving γ₁ by Δσ
    ds = DELTA_SIGMA
    s0_g = 0.25 + 0.5j * grid
    def contrib(sigma):
        rho = sigma + 1j * GAMMA1
        rho_bar = sigma - 1j * GAMMA1
        return (np.real(1.0 / (s0_g - rho) + 1.0 / rho) +
                np.real(1.0 / (s0_g - rho_bar) + 1.0 / rho_bar))
    dOm_fd = contrib(0.5 + ds) - contrib(0.5)
    dA_fd = ((E * (w * dOm_fd)) @ E.conj().T).real / (2.0 * np.pi)

    # Actual numerical lam_min response
    tau2 = (AB + dA_fd - P); tau2 = (tau2 + tau2.T) / 2.0
    lam_g = float(scipy.linalg.eigvalsh(tau2, subset_by_index=[0, 0])[0])
    actual_lam_shift = lam_g - lam0
    numeric_slope = actual_lam_shift / ds

    # Method 1: min diagonal H-F (first-order for each null mode individually)
    hf_diag = np.array([float(Vn[:, k] @ dA_fd @ Vn[:, k]) for k in range(dim)])
    min_hf_diag = float(hf_diag.min())
    hf_diag_slope = min_hf_diag / ds

    # Method 2: min eigenvalue of ΔA_null (includes off-diagonal coupling within null space)
    dA_null = Vn.T @ dA_fd @ Vn  # (dim, dim) matrix
    eig_null = scipy.linalg.eigvalsh(dA_null)
    min_null_eig = float(eig_null[0])
    null_eig_slope = min_null_eig / ds

    # Method 3: ground state H-F (what the original script computed)
    psi0 = V[:, 0]
    gs_hf = float(psi0 @ dA_fd @ psi0)
    gs_slope = gs_hf / ds

    # Method 4: trace of ΔA_null / dim (average null-space H-F)
    trace_null = float(np.trace(dA_null)) / max(1, dim)
    trace_slope = trace_null / ds

    print(f"lsq={lsq:4d} L={L:.3f} N={N} dim={dim}")
    print(f"  numeric_slope   = {numeric_slope:+.5e}")
    print(f"  gs_hf_slope     = {gs_slope:+.5e}  ratio={gs_slope/numeric_slope:.4f}")
    print(f"  min_diag_slope  = {hf_diag_slope:+.5e}  ratio={hf_diag_slope/numeric_slope:.4f}")
    print(f"  null_eig_slope  = {null_eig_slope:+.5e}  ratio={null_eig_slope/numeric_slope:.4f}")
    print(f"  trace_slope     = {trace_slope:+.5e}  ratio={trace_slope/numeric_slope:.4f}")
    print(f"  actual shift    = {actual_lam_shift:+.4e}  (lam0={lam0:.3e} lam_g={lam_g:.3e})")

    # Also report distribution of null H-F diagonal values
    hf_sorted = np.sort(hf_diag)
    print(f"  H-F diag null spread: min={hf_sorted[0]:.3e}  p25={np.percentile(hf_diag,25):.3e}  "
          f"p50={np.percentile(hf_diag,50):.3e}  max={hf_sorted[-1]:.3e}")

    return {
        'lsq': lsq, 'N': N, 'L': round(L, 5), 'dim': dim,
        'lam0': lam0, 'lam_g': lam_g,
        'numeric_slope': float(numeric_slope),
        'gs_hf_slope': float(gs_slope),
        'min_diag_hf_slope': float(hf_diag_slope),
        'null_eig_slope': float(null_eig_slope),
        'trace_slope': float(trace_slope),
        'null_eig_all': eig_null.tolist()[:10],  # first 10 eigenvalues of ΔA_null
        'hf_diag_min': float(hf_sorted[0]),
        'hf_diag_max': float(hf_sorted[-1]),
        'hf_diag_mean': float(np.mean(hf_diag)),
    }


def main():
    lsq_list = [20, 50, 100]
    results = []
    t0 = time.time()
    for lsq in lsq_list:
        N = 5 * lsq
        t1 = time.time()
        r = build_and_analyze(lsq, N)
        r['dt'] = round(time.time() - t1, 1)
        results.append(r)
        print(f"  time: {r['dt']}s\n", flush=True)

    out = {'gamma1': GAMMA1, 'delta_sigma': DELTA_SIGMA, 'results': results,
           'total_time': round(time.time() - t0, 1)}
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"Saved: {OUTFILE}")


if __name__ == "__main__":
    main()
