#!/usr/bin/env python3
"""
witness_exact.py — Exact witness test using the Hadamard-product formula for δΩ.

DERIVATION:
The A matrix uses Ω(ξ) = Re[ψ(1/4 + iξ/2)] - log π.

Each Riemann zero ρ = 1/2 + iγ (and its pair 1 - ρ̄ = 1/2 - iγ) contributes
to Ω(ξ) via the Hadamard product:

  Ω_ρ(ξ) = Re[1/(1/4 + iξ/2 - ρ)] + Re[1/(1/4 + iξ/2 - (1-ρ̄))]

For ρ = 1/2 + iγ:
  = Re[1/(-1/4 + i(ξ/2 - γ))] + Re[1/(-1/4 + i(ξ/2 + γ))]
  = -1/4 / [(1/4)² + (ξ/2 - γ)²] + -1/4 / [(1/4)² + (ξ/2 + γ)²]

For off-line ρ = 1/2 + Δσ + iγ (and pair 1/2 - Δσ - iγ):
  = (-1/4 - Δσ) / [(1/4 + Δσ)² + (ξ/2 - γ)²]
  + (-1/4 + Δσ) / [(1/4 - Δσ)² + (ξ/2 + γ)²]

The change δΩ = Ω_off - Ω_on is the EXACT perturbation.

KEY RESULT:
Near ξ = 2γ: δΩ ≈ +16Δσ > 0  (A becomes LESS negative → τ becomes LARGER)
Near ξ = -2γ: δΩ ≈ -16Δσ < 0 (A becomes MORE negative → τ becomes SMALLER)

Net effect depends on which modes (positive vs negative n) dominate in the near-null subspace.

IF the most destabilizing mode has n < 0 (receives the -2γ negative contribution),
then λ_min(τ̃) < λ_min(τ) — the witness IS accessible.

IF modes are symmetric, effects cancel. If n > 0 dominates, τ̃ would be LARGER.
"""
import time, json
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 60000
OUTFILE = '/Users/scones/projects/riemann/results/witness_exact.out'

GAMMA = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918720, 43.327073, 48.005151, 49.773832]

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

def build_base_and_modes(lsq, N):
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

    A = ((E * (w * Om)) @ E.conj().T).real / (2.0 * np.pi)
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

    tau = (A + B - P); tau = (tau + tau.T) / 2.0
    return tau, A, B, P, grid, w, E, L, ns

def delta_omega_exact(grid, gamma, dsigma):
    """
    Exact change δΩ(ξ) when zero γ moves from σ=1/2 to σ=1/2+Δσ.

    δΩ = [(-1/4 - Δσ) / ((1/4+Δσ)² + (ξ/2-γ)²)]
       + [(-1/4 + Δσ) / ((1/4-Δσ)² + (ξ/2+γ)²)]
       - [(-1/4) / ((1/4)² + (ξ/2-γ)²)]
       - [(-1/4) / ((1/4)² + (ξ/2+γ)²)]
    """
    xi2 = grid / 2.0
    a1 = 0.25 + dsigma; b1 = xi2 - gamma
    a2 = 0.25 - dsigma; b2 = xi2 + gamma
    a0 = 0.25

    off1 = (-a1) / (a1**2 + b1**2)
    off2 = (-a2) / (a2**2 + b2**2)
    on1  = (-a0) / (a0**2 + b1**2)
    on2  = (-a0) / (a0**2 + b2**2)

    return (off1 + off2) - (on1 + on2)

def main():
    print("Exact witness test: perturb first Riemann zero off critical line")
    print(f"Using EXACT Hadamard-product formula for δΩ")
    print()

    with open(OUTFILE, 'w') as f:
        f.write(f"witness_exact {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write("EXACT Hadamard δΩ formula\n\n")

    # Test at multiple (lsq, N) values
    configs = [
        (50, 250, "small"),
        (100, 500, "medium"),
        (200, 800, "γ₁-resonance"),
    ]

    delta_sigmas = [0.0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.4]

    for lsq, N, label in configs:
        print(f"\n--- lsq={lsq}, N={N} ({label}) ---")
        t0 = time.time()
        tau_base, A, B, P, grid, w, E, L, ns = build_base_and_modes(lsq, N)
        ev_base = scipy.linalg.eigvalsh(tau_base)
        lam_min_base = float(ev_base[0])
        print(f"  τ base: lam_min={lam_min_base:.4e}")

        with open(OUTFILE, 'a') as f:
            f.write(f"\nlsq={lsq} N={N} ({label})\n")
            f.write(f"  base lam_min={lam_min_base:.4e}\n")
            f.write(f"  {'Δσ':>8} {'lam_min':>12} {'Δlam':>12} {'sign':>10}\n")
            f.write("  " + "-" * 48 + "\n")

        for dsigma in delta_sigmas:
            if abs(dsigma) < 1e-10:
                lam_t = lam_min_base; dlam = 0.0
            else:
                # Exact perturbation matrix from moving γ₁
                dOm = delta_omega_exact(grid, GAMMA[0], dsigma)
                dA = ((E * (w * dOm)) @ E.conj().T).real / (2.0 * np.pi)
                tau_t = tau_base + dA
                tau_t = (tau_t + tau_t.T) / 2.0
                ev_t = scipy.linalg.eigvalsh(tau_t)
                lam_t = float(ev_t[0]); dlam = lam_t - lam_min_base

            sign = "NEGATIVE!" if lam_t < -1e-10 else "pos"
            print(f"  Δσ={dsigma:.3f}: lam_min={lam_t:.4e}  Δlam={dlam:+.4e}  {sign}")
            with open(OUTFILE, 'a') as f:
                f.write(f"  {dsigma:>8.4f} {lam_t:>12.4e} {dlam:>+12.4e} {sign:>10}\n")

        print(f"  [{time.time()-t0:.0f}s]")

    # Second test: perturb EACH zero individually to see which is most destabilizing
    print(f"\n--- Zero-by-zero destabilization test (lsq=100, N=500, Δσ=0.05) ---")
    lsq, N = 100, 500; dsigma = 0.05
    t0 = time.time()
    tau_base2, A2, B2, P2, grid2, w2, E2, L2, ns2 = build_base_and_modes(lsq, N)
    ev2 = scipy.linalg.eigvalsh(tau_base2)
    lam_base2 = float(ev2[0])
    print(f"  base lam_min={lam_base2:.4e}")

    with open(OUTFILE, 'a') as f:
        f.write(f"\nZero-by-zero test (lsq={lsq}, N={N}, Δσ={dsigma})\n")
        f.write(f"  base lam_min={lam_base2:.4e}\n")

    for k, gamma in enumerate(GAMMA):
        dOm = delta_omega_exact(grid2, gamma, dsigma)
        dA = ((E2 * (w2 * dOm)) @ E2.conj().T).real / (2.0 * np.pi)
        tau_t = tau_base2 + dA; tau_t = (tau_t + tau_t.T) / 2.0
        ev_t = scipy.linalg.eigvalsh(tau_t)
        lam_t = float(ev_t[0]); dlam = lam_t - lam_base2
        sign = "NEG!" if lam_t < -1e-10 else "pos"
        line = f"  γ_{k+1}={gamma:.3f}: lam_min={lam_t:.4e}  Δlam={dlam:+.4e}  {sign}"
        print(line)
        with open(OUTFILE, 'a') as f:
            f.write(line + "\n")

    print(f"  [{time.time()-t0:.0f}s]")
    print(f"\nOutput: {OUTFILE}")

if __name__ == "__main__":
    main()
