#!/usr/bin/env python3
"""
tau_spectrum.py — Extract the low-eigenvalue spectrum of τ and compare to Riemann zeros.

HYPOTHESIS: The small positive eigenvalues of τ encode the Riemann zero imaginary parts.

If τ's eigenvalue at index k (near zero) is proportional to 1/γ_k², this would be
a new spectral encoding of the zeros.

Comparison approach:
1. Compute full eigenvalue spectrum of τ at several (lsq, N) values
2. Extract the ~dim(null_space) smallest eigenvalues
3. Map them to predicted zero locations via n_k = γ_k * L / (2π)
4. Check if eigenvalue ordering matches zero ordering

Also check: do eigenvalue GAPS match zero-spacing statistics (GUE behavior)?

This would be a completely new result: the zeros encoded in the spectrum of τ.
"""
import time, json
import numpy as np
from scipy.special import digamma as _digamma
import scipy.linalg
import mpmath as mp

mp.mp.dps = 30
TWO_PI = 2.0 * np.pi
LOGPI  = float(np.log(np.pi))
G      = 50000
OUTFILE = '/Users/scones/projects/riemann/results/tau_spectrum.out'

# Known Riemann zero imaginary parts
GAMMA = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062,
         37.586178, 40.918720, 43.327073, 48.005151, 49.773832,
         52.970321, 56.446248, 59.347044, 60.831779, 65.112544,
         67.079811, 69.546402, 72.067158, 75.704691, 77.144840]

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

def build_tau(lsq, N):
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
    return tau, L, ns

def main():
    print("τ spectrum analysis: small eigenvalues vs Riemann zeros")
    print()
    with open(OUTFILE, 'w') as f:
        f.write(f"tau_spectrum {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    configs = [(50, 250), (100, 500)]

    for lsq, N in configs:
        print(f"--- lsq={lsq}, N={N} ---")
        t0 = time.time()
        tau, L, ns = build_tau(lsq, N)
        ev = scipy.linalg.eigvalsh(tau)

        # Predicted mode indices for each zero: n_γ = γ * L / (2π)
        n_gammas = [g * L / TWO_PI for g in GAMMA]

        # Extract 40 smallest eigenvalues (both positive and those near zero)
        ev_sorted = np.sort(ev)
        ev_pos = ev_sorted[ev_sorted > 0][:40]
        ev_neg = ev_sorted[ev_sorted < 0]
        lam_min = float(ev_sorted[0])
        dim_near_zero = int(np.sum(np.abs(ev_sorted) < 0.02))

        print(f"  N_modes={len(ev)}, lam_min={lam_min:.4e}, dim(near-zero)={dim_near_zero}")
        print(f"  First 10 positive eigenvalues:")
        for i, e in enumerate(ev_pos[:10]):
            print(f"    λ[{i}] = {e:.6e}")
        print()
        print(f"  Zero-mode correspondence: n_γ = γ × L/(2π)")
        for k, (g, n_g) in enumerate(zip(GAMMA[:10], n_gammas[:10])):
            print(f"  γ_{k+1}={g:.3f}: n_γ={n_g:.2f}  (should correspond to mode n≈{n_g:.0f})")

        # Check if eigenvalue distribution has GUE statistics
        # GUE prediction: level spacing distribution follows Wigner surmise p(s) = (π/2)s·exp(-πs²/4)
        ev_in_range = ev_pos[:min(30, len(ev_pos))]
        if len(ev_in_range) > 2:
            spacings = np.diff(ev_in_range)
            mean_sp = np.mean(spacings)
            norm_sp = spacings / mean_sp  # normalized spacings
            # Kolmogorov-Smirnov test vs Wigner surmise (rough)
            # Expected: mean_spacing ≈ mean(γ_k - γ_{k-1}) × some factor
            print(f"\n  Eigenvalue spacing statistics (first 30):")
            print(f"  mean spacing = {mean_sp:.4e}")
            print(f"  std/mean = {np.std(spacings)/mean_sp:.3f}  (GUE: ~0.52, Poisson: ~1.0)")

        with open(OUTFILE, 'a') as f:
            f.write(f"lsq={lsq} N={N} L={L:.4f}\n")
            f.write(f"  lam_min={lam_min:.4e}  dim(near-zero)={dim_near_zero}\n")
            f.write(f"  First 20 positive eigenvalues:\n")
            for i, e in enumerate(ev_pos[:20]):
                f.write(f"    λ[{i:2d}] = {e:.8e}\n")
            f.write(f"  n_γ (zero mode indices): {[f'{n:.2f}' for n in n_gammas[:10]]}\n")
            if len(ev_in_range) > 2:
                f.write(f"  mean_spacing={mean_sp:.4e}  std/mean={np.std(spacings)/mean_sp:.3f}\n")
            f.write(f"  [{time.time()-t0:.0f}s]\n\n")

        print(f"  [{time.time()-t0:.0f}s]")

    # Check eigenvalue ratio test: if λ_k ~ 1/γ_k², then λ_k × γ_k² = const
    print(f"\n--- Eigenvalue-zero ratio test at lsq=100 ---")
    lsq, N = 100, 500
    tau, L, ns = build_tau(lsq, N)
    ev = scipy.linalg.eigvalsh(tau)
    ev_pos = np.sort(ev[ev > 0])[:20]

    print("  λ_k × γ_k²  (constant = Riemann zeros encoded in spectrum):")
    with open(OUTFILE, 'a') as f:
        f.write("Eigenvalue-zero ratio test (lsq=100, N=500):\n")
        f.write("  k  γ_k     λ_k         λ_k × γ_k²\n")
    for k in range(min(10, len(ev_pos), len(GAMMA))):
        ratio = ev_pos[k] * GAMMA[k]**2
        line = f"  k={k+1}: γ_{k+1}={GAMMA[k]:.3f}  λ_{k+1}={ev_pos[k]:.4e}  λ×γ²={ratio:.6f}"
        print(line)
        with open(OUTFILE, 'a') as f:
            f.write(f"  {k+1}  {GAMMA[k]:.3f}  {ev_pos[k]:.6e}  {ratio:.6f}\n")

    print(f"\nOutput: {OUTFILE}")

if __name__ == "__main__":
    main()
