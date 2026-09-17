#!/usr/bin/env python3
"""
Analyze the per-prime Hellmann-Feynman cascade from S2624 to understand
the structure of the margin bound.

Key insight from S2624: the running residual
  R_k = <ψ|W02-WR|ψ> + Sum_{q₁..q_k} h(q_j)
telescopes to match λ_min exactly, and the remaining term |(last prime's h value)|
predicts the next step.

This suggests the h(q) contributions form a SPECIFIC SEQUENCE where each term
is determined by the structure of the matrix and the eigenvector.

Can we bound the RATE OF CONVERGENCE and relate it to the empirical 0.44 exponent?

Per S2624, the key structural finding:
  - rank(Wp) <= 2M where M = #{prime powers <= c}
  - But ALL of that rank is USED to cancel the (W02-WR) negativity
  - The degree of cancellation precision determines the margin

Hypothesis: the margin is roughly
  margin ~ (W02-WR negative part magnitude) / (how well Wp can cancel it)

If (W02-WR) has negative eigenvalues ~ -O(1) (say, O(1) in absolute value),
and Wp must cancel them with M rank-2 matrices, then each per-prime
contribution is ~O(1/M) on average.

But the actual cancellation is exponentially precise due to the specific
structure (Fourier/Dirichlet form + prime powers). This could lead to
exponential decay in the margin.

Let's model:
  - |arch_negative_spectrum| ~ O(1)
  - M prime powers to cancel it
  - Perfect cancellation structure due to Fourier basis + prime powers
  - Margin = residual after telescoping = exp(-f(c,N))

The empirical law margin ~ K * c^(-0.44*N) in log scale:
  log(margin) ~ -0.44*N*log(c) = -0.44*N*log(c)
  This is POWER-LAW decay, not exponential decay.

Could it be margin ~ c^(-a*N) for some a ~ 0.44?
Then log(margin) = -a*N*log(c), which matches.

Structure analysis: what determines 'a'?
- a relates to how many "effective prime bits" are needed to cancel the arch negativity
- In Fourier space, 2*pi*n*log(q)/L determines the oscillation; there are M primes
- The dimension 2N+1 limits frequency resolution
- The "effective rank" of Wp is at most 2M, and it needs to cancel ~O(1) negativity over ~2N+1 dimensions

Possible scaling: a ~ (# primes)/(dimension) ~ (c/log(c))/(N)
At c=30: #primes ~= 10, N=8 -> a ~ 10/8 = 1.25, but empirically 0.44
This is too crude.

Better idea: the margin is determined by the SMALLEST singular value of the
"effective cancellation operator" that Wp can achieve. If this scales as
σ_min(Wp) / (spectrum amplitude) ~ 1/(c*N) or similar, we get power-law decay.
"""

import numpy as np
from scipy import linalg
import sys
import math

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms, prime_powers_up_to

def per_prime_cascade_analysis(c, N):
    """Detailed analysis of the per-prime Hellmann-Feynman cascade."""
    print(f"\n{'='*100}")
    print(f"Per-Prime Cascade Analysis: c={c}, N={N}")
    print(f"{'='*100}")

    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    # Solve for eigenvector of tau
    evals, evecs = np.linalg.eigh(tau)
    idx_min = np.argmin(np.abs(evals))
    lam_min = evals[idx_min]
    psi = evecs[:, idx_min]
    psi = psi / np.linalg.norm(psi)

    # Cascade computation
    h_arch = psi @ W02_WR @ psi
    running = h_arch

    print(f"\nTrue λ_min = {lam_min:.6e}")
    print(f"<ψ|W02-WR|ψ> = {h_arch:.6e}")

    # Per-prime contributions
    print(f"\n{'prime^k':>8} {'y=log(q)':>10} {'h(q)':>16} {'|h(q)|':>12} {'residual':>16} {'rate':>12}")
    print("-" * 100)

    prev_residual = h_arch
    cascade_residuals = []

    for (q, p), Wp_mat in zip(pdata, Wp_terms):
        h = -(psi @ Wp_mat @ psi)
        running += h
        residual = running - lam_min
        cascade_residuals.append(abs(residual))

        if abs(h) > 1e-100:
            rate = math.log(abs(residual) / abs(prev_residual)) if abs(prev_residual) > 1e-300 else 0
        else:
            rate = 0

        print(f"{q:>8} {math.log(q):>10.3f} {h:>16.6e} {abs(h):>12.6e} {residual:>16.6e} {rate:>12.4f}")
        prev_residual = residual

    # Analyze the decay rate
    cascade_residuals = np.array(cascade_residuals)
    cascade_residuals = cascade_residuals[cascade_residuals > 1e-300]  # Filter out numerical noise

    if len(cascade_residuals) > 1:
        log_residuals = np.log(cascade_residuals)
        # Fit to see if decay is exponential or power-law
        indices = np.arange(len(log_residuals))

        # Linear regression: log(residual) vs index
        slope_index = np.polyfit(indices, log_residuals, 1)[0]

        # Linear regression: log(residual) vs log(prime index)
        prime_indices = np.arange(len(pdata))
        slope_log_index = np.polyfit(np.log(prime_indices + 1), log_residuals, 1)[0]

        print(f"\n--- Decay Rate Analysis ---")
        print(f"Slope of log(residual) vs prime index: {slope_index:.4f} (exponential ~ -const*k)")
        print(f"Slope of log(residual) vs log(prime index): {slope_log_index:.4f} (power-law ~ -const*log(k))")

    # Compare with the empirical 0.44*N*log(c) law
    empirical_scale = -0.44 * N * math.log(c)
    print(f"\nEmpirical decay scale: -0.44*N*log(c) = {empirical_scale:.4f}")
    print(f"True log(λ_min): {math.log(abs(lam_min)) if lam_min > 1e-300 else 'near zero'}")

    # Analyze M (number of primes) vs N
    M = len(pdata)
    print(f"\nStructural parameters:")
    print(f"  M (# prime powers <= c): {M}")
    print(f"  N (frequency cutoff): {N}")
    print(f"  dim = 2N+1: {2*N+1}")
    print(f"  log(c): {math.log(c):.4f}")
    print(f"  M/N: {M/N:.4f}")
    print(f"  log(c)/N: {math.log(c)/N:.4f}")
    print(f"  N*log(c): {N * math.log(c):.4f}")

    return lam_min, pdata

def collect_many_points():
    """Collect margin data across many (c,N) points to fit the exponent."""
    print("\n\n" + "="*100)
    print("SCALING LAW FIT: margin ~ K * c^(-α*N)")
    print("="*100)

    data_points = []

    for c in [10, 13, 20, 30, 50]:
        for N in [4, 6, 8, 10]:
            try:
                W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
                Wp_total = sum(Wp_terms)
                tau = W02_WR - Wp_total
                evals = np.linalg.eigvalsh(tau)
                lam_min = evals[0]

                M = len(pdata)
                data_points.append({
                    'c': c, 'N': N, 'M': M, 'dim': 2*N+1, 'lambda_min': lam_min,
                    'log_lam': math.log(abs(lam_min)) if abs(lam_min) > 1e-300 else -300,
                })
            except Exception as e:
                print(f"  Skipped c={c}, N={N}: {e}")

    # Fit the exponent
    print(f"\n{'c':>4} {'N':>4} {'M':>4} {'dim':>4} {'λ_min':>16} {'log(λ)':>12}")
    for pt in data_points:
        print(f"{pt['c']:>4} {pt['N']:>4} {pt['M']:>4} {pt['dim']:>4} {pt['lambda_min']:>16.6e} {pt['log_lam']:>12.4f}")

    # Regression: fit log(lambda_min) = a0 + a1*N*log(c) to determine alpha
    if len(data_points) > 2:
        N_vals = np.array([p['N'] for p in data_points])
        c_vals = np.array([p['c'] for p in data_points])
        log_lam_vals = np.array([p['log_lam'] for p in data_points])

        X = np.column_stack([np.ones(len(data_points)), N_vals * np.log(c_vals)])
        coeffs = np.linalg.lstsq(X, log_lam_vals, rcond=None)[0]
        a0, a1 = coeffs

        print(f"\n--- Regression Results ---")
        print(f"log(λ_min) ≈ {a0:.4f} + {a1:.4f} * N*log(c)")
        print(f"Estimated exponent α: {-a1:.4f} (empirical 0.44)")
        print(f"This means: λ_min ~ exp({a0:.4f}) * c^({-a1:.4f}*N)")
        print(f"Or: λ_min ~ {np.exp(a0):.6e} * c^({-a1:.4f}*N)")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--c', type=int, default=30)
    parser.add_argument('--N', type=int, default=8)
    parser.add_argument('--sweep', action='store_true', help='Run sweep across multiple (c,N)')
    args = parser.parse_args()

    if args.sweep:
        collect_many_points()
    else:
        per_prime_cascade_analysis(args.c, args.N)
