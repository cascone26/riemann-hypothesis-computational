#!/usr/bin/env python3
"""
sensitivity_scaling.py — How does dlam_min/dσ scale with lsq?

The moved_zero_test showed: moving γ₁=14.135 by Δσ=0.001 gives lam_min=-1.08e-3 at lsq=50.
This is first-order perturbation. The KEY question: does |dlam_min/dσ| grow with lsq?

If |dlam_min/dσ| ~ L^β for some β>0, then in the limit lsq→∞:
ANY nonzero Δσ produces lam_min → -∞.
This would PROVE that τ ≥ 0 requires all zeros at σ=1/2 in the limit.

Test at lsq = 5, 8, 12, 20, 30, 50 with γ₁=14.135:
  - Δσ values: 0, 0.001, 0.002, 0.005 (to get linear slope)
  - Compute dlam_min/dσ by linear fit to small Δσ points
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
OUTFILE = os.path.expanduser("~/riemann/results/sensitivity_scaling.json")


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
    s0 = 0.25 + 0.5j * grid
    def contrib(sigma):
        rho = sigma + 1j * gamma
        rho_bar = sigma - 1j * gamma
        return (np.real(1.0 / (s0 - rho) + 1.0 / rho) +
                np.real(1.0 / (s0 - rho_bar) + 1.0 / rho_bar))
    delta_Om = contrib(sigma_new) - contrib(sigma_old)
    return ((E * (w * delta_Om)) @ E.conj().T).real / (2.0 * np.pi)


def measure_lam(arch, B, P, E, grid, w, delta_A=None):
    AB = arch + B + (delta_A if delta_A is not None else 0)
    tau = (AB - P)
    tau = (tau + tau.T) / 2.0
    ev = scipy.linalg.eigvalsh(tau)
    return float(ev[0])


def main():
    lsq_list = [5, 8, 12, 20, 30, 50, 80]
    # Use N=5*lsq (optimal from sc_vs_N study), min N=20
    delta_sigmas = [0.001, 0.002, 0.005, 0.010]  # linear regime

    results = []
    t0 = time.time()

    for lsq in lsq_list:
        N = max(20, 5 * lsq)
        L = float(mp.log(lsq))
        print(f"\nlsq={lsq}, N={N}, L={L:.4f}", flush=True)

        t1 = time.time()
        arch, B, P, E, grid, w, L, ns = build_components(lsq, N)

        # Baseline
        lam0 = measure_lam(arch, B, P, E, grid, w)
        print(f"  lam0 = {lam0:.4e}", flush=True)

        # Small Δσ values
        lam_vals = [lam0]
        for ds in delta_sigmas:
            dA = zero_delta_A(E, grid, w, 0.5 + ds, GAMMA1)
            lam = measure_lam(arch, B, P, E, grid, w, dA)
            lam_vals.append(lam)
            print(f"  Δσ={ds}: lam={lam:.4e}  Δlam/Δσ={(lam-lam0)/ds:.4e}", flush=True)

        # Linear fit to get dlam/dσ
        # Use points Δσ=0.001, 0.002, 0.005 for linearity (avoid 0.010 which might be non-linear)
        ds_arr = np.array([0] + delta_sigmas)
        lam_arr = np.array(lam_vals)
        # Linear fit: lam(σ₀+Δσ) = lam0 + slope*Δσ
        slope = np.polyfit(ds_arr[:4], lam_arr[:4], 1)[0]  # first 4 points (0, 0.001, 0.002, 0.005)
        print(f"  slope dlam/dσ = {slope:.4e}  (log|slope|/logL = {np.log10(abs(slope))/np.log10(L):.3f})", flush=True)

        results.append({
            'lsq': lsq, 'L': round(L, 6), 'N': N,
            'lam0': lam0,
            'lam_vals': lam_vals,
            'delta_sigmas': [0] + delta_sigmas,
            'dlam_dsigma': float(slope),
            'dt': round(time.time() - t1, 1),
        })

    # Fit slope vs L
    Ls = np.array([r['L'] for r in results])
    slopes = np.abs([r['dlam_dsigma'] for r in results])
    try:
        from scipy.optimize import curve_fit
        popt, pcov = curve_fit(lambda x, C, b: C * x**b, Ls, slopes, p0=[1.0, 2.0])
        perr = np.sqrt(np.diag(pcov))
        print(f"\nSlope scaling: |dlam/dσ| ~ {popt[0]:.4f} × L^{popt[1]:.4f} ± {perr[1]:.4f}")
    except Exception as e:
        print(f"Slope fit failed: {e}")
        popt = [None, None]

    out = {
        'lsq_list': lsq_list,
        'gamma1': GAMMA1,
        'results': results,
        'slope_fit': {'C': float(popt[0]) if popt[0] else None,
                      'beta': float(popt[1]) if popt[1] else None},
        'total_time': round(time.time() - t0, 1),
    }
    with open(OUTFILE, 'w') as f:
        json.dump(out, f, indent=2)
    print(f"\nSaved: {OUTFILE} in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
