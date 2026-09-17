#!/usr/bin/env python3
"""
EIGENVECTOR STRUCTURE ANALYSIS for Analytical Exponent Derivation
==================================================================

Key insight: The per-prime Hellmann-Feynman term h(q) = -<ψ|Wp_q|ψ>
depends critically on how the eigenvector ψ is distributed across
frequency modes n ∈ [-N, N].

If ψ is concentrated in a few low frequencies (like a prolate eigenfunction),
then h(q) will be small for most primes, making the cascade decay slow.

If ψ is spread across many frequencies with certain decay pattern,
the cascade might follow a power-law: log|h(k)| ~ -β*log(k).

This analysis measures:
1. The eigenvector concentration (how many "effective" frequencies it uses)
2. The correlation between ψ components and the prime/Fourier structure
3. How this relates to the observed power-law decay rate β
4. Whether the decay rate β determines the final exponent 0.44
"""

import numpy as np
import sys
import math
from math import log, exp, pi

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms

def analyze_eigenvector_concentration(c, N):
    """
    Analyze the Fourier concentration of the smallest-margin eigenvector.

    Returns:
    - ψ_components: values of ψ at each frequency n
    - effective_width: Rényi entropy estimate of how many frequencies matter
    - concentration_ratio: measure of spread (0=delta-like, 1=uniform)
    - frequency_decay: fit to |ψ_n| ~ n^(-γ) or exponential
    """
    print(f"\n{'='*100}")
    print(f"EIGENVECTOR CONCENTRATION: c={c}, N={N}")
    print(f"{'='*100}")

    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    # Solve for min eigenvector
    evals, evecs = np.linalg.eigh(tau)
    idx_min = np.argmin(np.abs(evals))
    lam_min = evals[idx_min]
    psi = evecs[:, idx_min]
    psi = psi / np.linalg.norm(psi)

    # Extract components by frequency
    psi_by_freq = {}
    for i in range(len(psi)):
        n = i - N
        psi_by_freq[n] = psi[i]

    # Analyze concentration
    abs_psi = np.abs(list(psi_by_freq.values()))
    psi_sq = abs_psi ** 2

    # Rényi entropy (order 2): H_2 = -log(sum_n psi_n^4)
    # Effective support: D_eff = exp(H_2) estimates how many frequencies "effectively" contribute
    sum_psi_4 = np.sum(psi_sq ** 2)
    renyi_entropy = -log(sum_psi_4 + 1e-300)
    effective_width = exp(renyi_entropy)

    # Concentration: ratio of power in central modes vs total
    central_width = N // 2 + 1
    power_central = np.sum(psi_sq[N - central_width:N + central_width + 1])
    concentration_ratio = power_central

    print(f"Eigenvector properties:")
    print(f"  Eigenvalue: λ_min = {lam_min:.6e}")
    print(f"  Norm: ||ψ|| = {np.linalg.norm(psi):.6f}")
    print(f"  Rényi entropy (order 2): H_2 = {renyi_entropy:.4f}")
    print(f"  Effective width (D_eff = exp(H_2)): {effective_width:.2f} frequencies")
    print(f"  Concentration (power in central {central_width} freqs): {concentration_ratio:.4f}")
    print(f"  Max component: |ψ_{max(psi_by_freq, key=lambda k: abs(psi_by_freq[k]))}| = {max(abs_psi):.6e}")

    # Print top components
    print(f"\n  Top 10 Fourier components (by power):")
    freqs_by_power = sorted(psi_by_freq.items(), key=lambda x: abs(x[1]), reverse=True)[:10]
    for n, val in freqs_by_power:
        print(f"    ψ_{n:+3d} = {val:+.6e}  (power {val**2:.6e})")

    # Try to fit decay of |ψ_n| with n
    # Model 1: |ψ_n| ~ A * n^(-α)  (power-law)
    # Model 2: |ψ_n| ~ A * exp(-β*n)  (exponential)

    nonzero_freqs = [(n, val) for n, val in psi_by_freq.items() if abs(val) > 1e-15]
    if len(nonzero_freqs) > 3:
        n_vals = np.array([n for n, _ in nonzero_freqs])
        psi_vals = np.array([val for _, val in nonzero_freqs])
        abs_psi_vals = np.abs(psi_vals)

        # Power-law fit: log|ψ| ~ α*log|n| + const
        # Only use |n| > 0 to avoid log(0)
        nonzero_mask = n_vals != 0
        if np.sum(nonzero_mask) > 2:
            n_nz = np.abs(n_vals[nonzero_mask])
            psi_nz = abs_psi_vals[nonzero_mask]
            log_n = np.log(n_nz)
            log_psi = np.log(psi_nz + 1e-300)

            # Linear regression: log(psi) vs log(n)
            coeffs_power = np.polyfit(log_n, log_psi, 1)
            alpha_fit = coeffs_power[0]
            const_fit = coeffs_power[1]

            print(f"\n  Power-law fit: log|ψ_n| ≈ {alpha_fit:.3f}*log|n| + {const_fit:.3f}")
            print(f"    i.e., |ψ_n| ~ {exp(const_fit):.6e} * n^({alpha_fit:.3f})")

            # Predict: what is the effective bandwidth given this decay?
            # Integral ∫_1^∞ n^(-2α) dn converges if α > 1/2
            if alpha_fit > 0.5:
                print(f"    → Eigenvector decays faster than n^(-1/2), so effective bandwidth is finite")
                # Effective bandwidth: ∫_1^N n^(-2α) dn ~ N^(1-2α) / (1-2α)
                eff_bw_integral = N**(1 - 2*alpha_fit) / abs(1 - 2*alpha_fit)
                print(f"    → Effective bandwidth ~ {eff_bw_integral:.2f}")
            else:
                print(f"    → Eigenvector decays slowly; effective bandwidth is large")

    return psi_by_freq, effective_width, concentration_ratio, lam_min

def analyze_per_prime_contributions(c, N):
    """
    Compute the actual per-prime Hellmann-Feynman contributions h(q)
    and analyze their decay pattern.

    Hypothesis: h(q) depends on the Fourier structure of ψ and the
    oscillation pattern of the q-th prime's contribution in Fourier space.
    """
    print(f"\n{'='*100}")
    print(f"PER-PRIME CONTRIBUTIONS: c={c}, N={N}")
    print(f"{'='*100}")

    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    # Solve for min eigenvector
    evals, evecs = np.linalg.eigh(tau)
    idx_min = np.argmin(np.abs(evals))
    lam_min = evals[idx_min]
    psi = evecs[:, idx_min]
    psi = psi / np.linalg.norm(psi)

    print(f"\nλ_min = {lam_min:.6e}")

    # Compute per-prime contributions
    print(f"\n{'q':>8} {'log(q)':>10} {'p':>5} {'log(p)':>8} {'h(q)':>16} {'log|h(q)|':>12} {'Rate (β)':>10}")
    print("-" * 100)

    h_values = []
    log_q_values = []
    for (q, p), Wp_mat in zip(pdata, Wp_terms):
        h = -(psi @ Wp_mat @ psi)
        h_values.append(abs(h))
        log_q_values.append(log(q))

        if abs(h) > 1e-300:
            log_h = log(abs(h))
        else:
            log_h = -700

        print(f"{q:>8} {log(q):>10.3f} {p:>5} {log(p):>8.3f} {h:>16.6e} {log_h:>12.2f}")

    # Fit the decay: log|h(q)| vs log(q)
    h_values = np.array(h_values)
    log_q_values = np.array(log_q_values)

    nonzero_mask = h_values > 1e-300
    if np.sum(nonzero_mask) > 2:
        h_nz = h_values[nonzero_mask]
        log_q_nz = log_q_values[nonzero_mask]
        log_h_nz = np.log(h_nz)

        # Linear regression: log|h(q)| vs log(q)
        # Model: log|h| = a + b*log(q)  =>  h ~ q^b
        coeffs = np.polyfit(log_q_nz, log_h_nz, 1)
        b = coeffs[0]  # power-law exponent
        a = coeffs[1]

        print(f"\n--- Decay Analysis ---")
        print(f"Per-prime h(q) follows: log|h(q)| ≈ {a:.3f} + {b:.3f}*log(q)")
        print(f"  i.e., h(q) ~ {exp(a):.6e} * q^({b:.3f})")
        print(f"  Power-law exponent for per-prime decay: {b:.3f}")

        # Now, if the cascade sums h(q) over all q, what's the total decay?
        # Σ_q h(q) ~ Σ_q q^b
        # For b < -1: Σ_q q^b ~ converges, sum ~ ∫ q^b dq from 1 to c ~ c^(b+1)/(b+1)
        # For -1 < b < 0: Σ_q q^b ~ ∫ q^b dq ~ c^(b+1)/(b+1), also power-law in c
        # For b > 0: Σ_q q^b diverges (each prime contributes positively)

        print(f"\nTotal cascade contribution: Σ_{{q≤c}} h(q)")
        if b < -1:
            total_scale_exponent = b + 1
            print(f"  Since b = {b:.3f} < -1, the sum converges")
            print(f"  Total ~ ∫_1^c q^{b:.3f} dq ~ c^({total_scale_exponent:.3f})")
        elif -1 < b < 0:
            total_scale_exponent = b + 1
            print(f"  Since -1 < b = {b:.3f} < 0, the sum is dominated by large q")
            print(f"  Total ~ ∫_1^c q^{b:.3f} dq ~ c^({total_scale_exponent:.3f})")
        else:
            print(f"  Since b = {b:.3f} >= 0, individual h(q) do not decay with q")
            print(f"  This is unexpected! h(q) should decay for the cascade to converge.")

        # Relate to the final margin
        # margin = background + Σ h(q)
        # ~ background + c^(b+1)
        background = psi @ W02_WR @ psi
        print(f"\nMargin reconstruction:")
        print(f"  Background <ψ|W02-WR|ψ> = {background:.6e}")
        print(f"  Cascade sum ~ c^({total_scale_exponent:.3f}) = {c**(total_scale_exponent):.6e}")
        print(f"  Final margin ~ {background:.6e} + {c**(total_scale_exponent):.6e} = ...")
        print(f"  True λ_min = {lam_min:.6e}")

    return h_values, log_q_values

def relate_to_exponent(c, N):
    """
    Try to relate the per-prime decay rate to the empirical -0.44*N exponent.

    Key idea: If h(q) ~ q^(-1.5) (for example), then:
      Σ_q h(q) ~ Σ_q q^(-1.5) ~ c^(-0.5)

    If this pattern holds across different (c,N), and if the decay exponent
    depends on N (say, b = -1.5 - 0.5*N), then:
      Σ_q h(q) ~ c^(b+1) = c^(-0.5-0.5*N)

    This would give log(margin) ~ (-0.5-0.5*N)*log(c), a close match!
    """
    print(f"\n{'='*100}")
    print(f"EXPONENT RECONSTRUCTION: c={c}, N={N}")
    print(f"{'='*100}")

    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    evals, evecs = np.linalg.eigh(tau)
    idx_min = np.argmin(np.abs(evals))
    lam_min = evals[idx_min]
    psi = evecs[:, idx_min]
    psi = psi / np.linalg.norm(psi)

    # Per-prime contributions
    h_values = []
    indices = []
    for idx, ((q, p), Wp_mat) in enumerate(zip(pdata, Wp_terms)):
        h = abs(-(psi @ Wp_mat @ psi))
        if h > 1e-300:
            h_values.append(h)
            indices.append(idx + 1)  # 1-indexed

    # The cascade sum: running sum of h values
    cascade_running = np.array(h_values)
    cascade_running = np.cumsum(cascade_running)

    # Fit power-law: cascade_sum(k) ~ A * k^(-α)
    # where k is the number of primes included so far
    if len(cascade_running) > 2:
        log_indices = np.log(np.array(indices))
        log_cascade = np.log(cascade_running)

        # Linear regression
        coeffs = np.polyfit(log_indices, log_cascade, 1)
        slope = coeffs[0]  # This is -α in the power-law decay
        intercept = coeffs[1]

        print(f"\nCascade sum decay:")
        print(f"  log(sum_k h(k)) ≈ {intercept:.3f} + {slope:.3f}*log(k)")
        print(f"  i.e., sum_k h(k) ~ {exp(intercept):.6e} * k^({slope:.3f})")
        print(f"  Power-law decay exponent: {slope:.3f}")

        # This is the "cascade decay rate β" mentioned in S2626
        cascade_decay_beta = -slope  # Negative slope means decay

        print(f"\n  Cascade decay rate β = {cascade_decay_beta:.3f}")
        print(f"  Empirical β (c={c},N={N}) ≈ 11.6")

        # Relate to margin:
        # If cascade sum ~ k^(slope), and k ~ M ~ c/ln(c), then:
        # cascade sum ~ (c/ln(c))^(slope) ~ c^(slope)
        # So log(margin) ~ slope * log(c)

        predicted_exponent = slope
        empirical_exponent = -0.44 * N

        print(f"\nExponent comparison:")
        print(f"  Predicted (from cascade decay): {predicted_exponent:.4f}")
        print(f"  Empirical (-0.44*N at N={N}): {empirical_exponent:.4f}")
        print(f"  Match: {abs(predicted_exponent - empirical_exponent) < 0.1}")

def main():
    test_points = [(30, 4), (30, 6), (30, 8), (30, 10), (13, 8), (100, 8)]

    for c, N in test_points:
        psi_by_freq, eff_width, conc, lam = analyze_eigenvector_concentration(c, N)
        h_vals, log_q = analyze_per_prime_contributions(c, N)
        relate_to_exponent(c, N)

if __name__ == '__main__':
    main()
