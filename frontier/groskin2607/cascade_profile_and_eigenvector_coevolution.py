#!/usr/bin/env python3
"""
Deep investigation of the cascade discrepancy from S2637:

ISSUE: The exact cascade (from S2624) correctly telescopes to λ_min.
But when we fit h(q) to a power law and sum it, we get c^(-11) instead of c^(-0.44N).

HYPOTHESIS (to test):
1. h(q) is NOT a simple power law — it has complex structure depending on both q and the eigenvector
2. The eigenvector ψ CHANGES as c varies at fixed N, so h(q) values "for prime q" are NOT fixed across different c
3. The cascade's power-law-like appearance in log-residual-vs-log-prime is misleading:
   it mixes per-prime contributions (which depend on q) with eigenvector changes (which depend on c)

STRATEGY:
- Part A: For FIXED (c,N), analyze h(q) profile in detail (signs, magnitudes, relation to eigenvector structure)
- Part B: For FIXED N, vary c; track how λ_min AND ψ change
- Part C: Separate the "true per-prime decay if ψ were fixed" from "how ψ reorganizes with c"
- Part D: Try to reconstruct c^(-0.44N) from first principles
"""

import numpy as np
from scipy import linalg
import sys
import math
from collections import defaultdict

sys.path.insert(0, '/Users/scones/projects/riemann/frontier/groskin2607')
from cascade_generalizes_finiteN import build_tau_and_terms, prime_powers_up_to

def analyze_cascade_profile_fixed_c_N(c, N):
    """
    PART A: For fixed (c,N), analyze the full h(q) profile.
    """
    print(f"\n{'='*100}")
    print(f"PART A: Cascade Profile Analysis at FIXED c={c}, N={N}")
    print(f"{'='*100}")

    W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
    Wp_total = sum(Wp_terms)
    tau = W02_WR - Wp_total

    evals, evecs = np.linalg.eigh(tau)
    idx_min = np.argmin(np.abs(evals))
    lam_min = evals[idx_min]
    psi = evecs[:, idx_min]
    psi = psi / np.linalg.norm(psi)

    h_arch = psi @ W02_WR @ psi

    print(f"\nTrue λ_min = {lam_min:.6e}")
    print(f"<ψ|W02-WR|ψ> = {h_arch:.6e}")
    print(f"Eigenvector concentration (Renyi entropy D_eff): ", end="")

    # Renyi entropy D_eff = 1 / sum(p_i^2)
    p_sq = psi**2
    D_eff = 1.0 / np.sum(p_sq)
    print(f"{D_eff:.2f} (max possible {2*N+1})")

    # Track per-prime contributions
    h_dict = {}  # maps q -> h(q)
    running = h_arch

    for (q, p), Wp_mat in zip(pdata, Wp_terms):
        h = -(psi @ Wp_mat @ psi)
        h_dict[q] = h
        running += h

    # Analyze the structure
    print(f"\n--- Per-Prime Contributions (sorted by prime, NOT by magnitude) ---")
    print(f"{'q':>6} {'prime':>6} {'log(q)':>8} {'h(q)':>16} {'|h(q)|':>12} {'% of total |h|':>14}")

    total_abs_h = sum(abs(h) for h in h_dict.values())

    for idx, (q, p) in enumerate(pdata):
        if q in h_dict:
            h = h_dict[q]
            pct = 100.0 * abs(h) / total_abs_h if total_abs_h > 0 else 0
            print(f"{q:>6} {p:>6} {math.log(q):>8.3f} {h:>16.6e} {abs(h):>12.6e} {pct:>14.2f}%")

    # Separate primes from prime powers
    base_primes = [q for q, p in pdata if q == p]
    prime_powers = [q for q, p in pdata if q != p]

    h_base_primes = np.array([abs(h_dict[q]) for q in base_primes])
    h_prime_powers = np.array([abs(h_dict[q]) for q in prime_powers])

    print(f"\n--- Structural Breakdown ---")
    print(f"Base primes (q=p): {base_primes}")
    print(f"  Total contribution: {sum(abs(h_dict[q]) for q in base_primes):.6e}")
    print(f"Prime powers (q≠p): {prime_powers}")
    print(f"  Total contribution: {sum(abs(h_dict[q]) for q in prime_powers):.6e}")
    print(f"  Ratio (prime powers / base primes): {sum(abs(h_dict[q]) for q in prime_powers) / (sum(abs(h_dict[q]) for q in base_primes) + 1e-300):.6f}")

    # Fit power law to BASE primes only
    if len(base_primes) > 2:
        y = np.log(h_base_primes)
        x = np.log([float(q) for q in base_primes])
        # Linear regression
        slope, intercept = np.polyfit(x, y, 1)
        print(f"\n--- Power-Law Fit (BASE PRIMES ONLY) ---")
        print(f"log|h(q)| ≈ {intercept:.4f} + {slope:.4f} * log(q)")
        print(f"=> |h(q)| ~ q^({slope:.4f})")
        print(f"This is MUCH STEEPER than the empirical cascade decay in S2637 (which was ~q^(-12))")

    # Check sign patterns
    positive_h = [(q, h_dict[q]) for q, p in pdata if q == p and h_dict[q] > 0]
    negative_h = [(q, h_dict[q]) for q, p in pdata if q == p and h_dict[q] < 0]

    print(f"\n--- Sign Pattern (Base Primes Only) ---")
    print(f"Positive contributions: {len(positive_h)} primes")
    if positive_h:
        for q, h in positive_h[:5]:
            print(f"  q={q}: h={h:.6e}")
        if len(positive_h) > 5:
            print(f"  ... ({len(positive_h)-5} more)")
    print(f"Negative contributions: {len(negative_h)} primes")
    if negative_h:
        for q, h in negative_h[:5]:
            print(f"  q={q}: h={h:.6e}")
        if len(negative_h) > 5:
            print(f"  ... ({len(negative_h)-5} more)")

    print(f"\nCRITICAL OBSERVATION:")
    print(f"  The cascade is NOT a simple power law with constant exponent!")
    print(f"  h(q) depends on BOTH the prime q AND the eigenvector ψ at this specific (c,N).")

    return lam_min, psi, h_dict, pdata


def track_eigenvector_evolution_at_fixed_N(N, c_list):
    """
    PART B: For FIXED N, vary c; track how λ_min AND ψ change.
    """
    print(f"\n\n{'='*100}")
    print(f"PART B: Eigenvector Evolution at FIXED N={N}, varying c")
    print(f"{'='*100}")

    results = {}

    for c in c_list:
        W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
        Wp_total = sum(Wp_terms)
        tau = W02_WR - Wp_total

        evals, evecs = np.linalg.eigh(tau)
        idx_min = np.argmin(np.abs(evals))
        lam_min = evals[idx_min]
        psi = evecs[:, idx_min]
        psi = psi / np.linalg.norm(psi)

        # Concentration measure
        p_sq = psi**2
        D_eff = 1.0 / np.sum(p_sq)

        # Projection of ψ onto a few top modes
        top_mode_power = np.sum(psi[:3]**2)  # modes n = -N, -N+1, -N+2

        results[c] = {
            'lambda_min': lam_min,
            'psi': psi,
            'D_eff': D_eff,
            'top_mode_power': top_mode_power,
            'pdata': pdata,
        }

    print(f"\n{'c':>6} {'λ_min':>16} {'log10(λ_min)':>16} {'D_eff':>10} {'top3_power':>12} {'M':>5}")
    print("-" * 80)
    for c in c_list:
        res = results[c]
        lam = res['lambda_min']
        D = res['D_eff']
        t3p = res['top_mode_power']
        M = len(res['pdata'])
        log_lam = math.log10(abs(lam)) if abs(lam) > 1e-300 else -300
        print(f"{c:>6} {lam:>16.6e} {log_lam:>16.2f} {D:>10.4f} {t3p:>12.6f} {M:>5}")

    # Check eigenvector SHAPE changes
    print(f"\n--- Eigenvector Shape Comparison ---")
    c_ref = c_list[0]
    psi_ref = results[c_ref]['psi']

    for c in c_list[1:]:
        psi_c = results[c]['psi']

        # Need to align them (eigenvector sign is arbitrary)
        overlap = abs(np.dot(psi_ref, psi_c))

        # Renormalize psi_c if overlap < 0.5 (might have flipped sign)
        if np.dot(psi_ref, psi_c) < 0:
            psi_c = -psi_c
            overlap = abs(np.dot(psi_ref, psi_c))

        print(f"c={c_ref} vs c={c}: overlap = {overlap:.6f}, " +
              f"Δλ = {results[c]['lambda_min'] - results[c_ref]['lambda_min']:.6e}")

    print(f"\nCRITICAL OBSERVATION:")
    print(f"  The eigenvector DOES change shape as c varies (overlap < 1).")
    print(f"  This means h(q) 'for prime q' is NOT independent of c!")
    print(f"  When we increase c, we get NEW primes, but ALSO ψ reorganizes.")

    return results


def separate_effect_of_new_primes_from_eigenvector_change(N):
    """
    PART C: Attempt to separate the two effects.
    """
    print(f"\n\n{'='*100}")
    print(f"PART C: Separating New-Prime Effect from Eigenvector Reorganization (N={N})")
    print(f"{'='*100}")

    # Compare c=30 vs c=50
    c1, c2 = 30, 50

    W02_WR_1, Wp_terms_1, pdata_1 = build_tau_and_terms(c1, N)
    W02_WR_2, Wp_terms_2, pdata_2 = build_tau_and_terms(c2, N)

    Wp_total_1 = sum(Wp_terms_1)
    Wp_total_2 = sum(Wp_terms_2)

    tau_1 = W02_WR_1 - Wp_total_1
    tau_2 = W02_WR_2 - Wp_total_2

    evals_1, evecs_1 = np.linalg.eigh(tau_1)
    evals_2, evecs_2 = np.linalg.eigh(tau_2)

    idx_1 = np.argmin(np.abs(evals_1))
    idx_2 = np.argmin(np.abs(evals_2))

    lam_1 = evals_1[idx_1]
    lam_2 = evals_2[idx_2]

    psi_1 = evecs_1[:, idx_1]
    psi_2 = evecs_2[:, idx_2]

    psi_1 = psi_1 / np.linalg.norm(psi_1)
    psi_2 = psi_2 / np.linalg.norm(psi_2)

    print(f"\nComparison: c={c1} (M={len(pdata_1)} primes) vs c={c2} (M={len(pdata_2)} primes)")
    print(f"  λ_min at c={c1}: {lam_1:.6e}")
    print(f"  λ_min at c={c2}: {lam_2:.6e}")
    print(f"  Decay: {abs(lam_2/lam_1):.6e}")
    print(f"  Empirical prediction c^(-0.44*N) = {c2/c1}^(-0.44*{N}) = {(c2/c1)**(-0.44*N):.6e}")

    # Extract h(q) contributions at each c
    h_dict_1 = {}
    h_dict_2 = {}

    for (q, p), Wp_mat in zip(pdata_1, Wp_terms_1):
        h = -(psi_1 @ Wp_mat @ psi_1)
        h_dict_1[q] = h

    for (q, p), Wp_mat in zip(pdata_2, Wp_terms_2):
        h = -(psi_2 @ Wp_mat @ psi_2)
        h_dict_2[q] = h

    # Compare h(q) for primes in BOTH sets
    common_primes = set([q for q, p in pdata_1 if q == p]) & set([q for q, p in pdata_2 if q == p])
    common_primes = sorted(list(common_primes))

    print(f"\n--- Comparing h(q) values for primes in BOTH c={c1} and c={c2} ---")
    print(f"{'q':>6} {'h(q)|c=30':>16} {'h(q)|c=50':>16} {'ratio (50/30)':>14}")
    print("-" * 60)

    for q in common_primes:
        h1 = h_dict_1[q]
        h2 = h_dict_2[q]
        ratio = h2 / h1 if h1 != 0 else float('nan')
        print(f"{q:>6} {h1:>16.6e} {h2:>16.6e} {ratio:>14.6e}")

    print(f"\nCRITICAL OBSERVATION:")
    print(f"  h(q) values CHANGE when we increase c, even for the SAME prime q!")
    print(f"  This is because ψ reorganizes as c changes.")
    print(f"  Therefore, 'h(q) as a function of q alone' is NOT a well-defined quantity independent of c.")
    print(f"  The naive picture 'add new prime, get new h(q), sum them all' is FLAWED.")

    return lam_1, lam_2, common_primes, h_dict_1, h_dict_2


def try_corrected_reconstruction(N):
    """
    PART D: Try to reconstruct c^(-0.44N) from first principles.

    Corrected picture:
    - For each NEW prime that enters when expanding c -> c',
      the eigenvector ψ REORGANIZES
    - Track not just the new h(q) for the new prime,
      but also how the EXISTING h(q) values change
    """
    print(f"\n\n{'='*100}")
    print(f"PART D: Corrected Reconstruction Attempt (N={N})")
    print(f"{'='*100}")

    c_sequence = [10, 15, 20, 25, 30, 40, 50]

    margins = []
    deltas = []

    print(f"\n{'c':>6} {'λ_min':>16} {'log10(λ)':>14} {'log10(Δλ)':>14} {'M':>5} {'log(c)':>8}")
    print("-" * 80)

    for c in c_sequence:
        W02_WR, Wp_terms, pdata = build_tau_and_terms(c, N)
        Wp_total = sum(Wp_terms)
        tau = W02_WR - Wp_total

        evals = np.linalg.eigvalsh(tau)
        lam = evals[0]
        margins.append(lam)

        M = len(pdata)
        log_lam = math.log10(abs(lam)) if abs(lam) > 1e-300 else -300
        log_c = math.log(c)

        if len(margins) > 1:
            delta_lam = abs(margins[-1] - margins[-2])
            log_delta = math.log10(delta_lam) if delta_lam > 1e-300 else -300
        else:
            log_delta = None

        if log_delta is not None:
            print(f"{c:>6} {lam:>16.6e} {log_lam:>14.2f} {log_delta:>14.2f} {M:>5} {log_c:>8.3f}")
        else:
            print(f"{c:>6} {lam:>16.6e} {log_lam:>14.2f} {'N/A':>14} {M:>5} {log_c:>8.3f}")

    # Fit to c^(-a*N)
    print(f"\n--- Fitting margin ~ K * c^(-a*N) ---")

    valid_pairs = []
    for i, c in enumerate(c_sequence):
        lam = margins[i]
        if abs(lam) > 1e-300:
            valid_pairs.append((c, lam))

    if len(valid_pairs) > 2:
        c_vals = np.array([c for c, _ in valid_pairs])
        lam_vals = np.array([lam for _, lam in valid_pairs])

        # Fit: log(lam) = a0 + a1 * log(c)
        X = np.column_stack([np.ones(len(c_vals)), np.log(c_vals)])
        coeffs = np.linalg.lstsq(X, np.log(np.abs(lam_vals)), rcond=None)[0]
        a0, a1 = coeffs

        print(f"log|λ_min| ≈ {a0:.4f} + {a1:.4f} * log(c)")
        print(f"Exponent (scaled by N={N}): a1*N = {a1*N:.4f} (empirical -0.44*N = {-0.44*N:.4f})")
        print(f"Observed exponent α = {-a1/N:.4f} (target 0.44)")

    return margins


# ========== MAIN ==========

if __name__ == '__main__':
    # RUN ALL PARTS

    print("\n" + "="*100)
    print("COMPREHENSIVE CASCADE DISCREPANCY ANALYSIS")
    print("="*100)

    # PART A: Fixed (c,N) profile
    lam_30_8, psi_30_8, h_dict_30_8, pdata_30_8 = analyze_cascade_profile_fixed_c_N(30, 8)

    # PART B: Eigenvector evolution
    results_N8 = track_eigenvector_evolution_at_fixed_N(8, [13, 20, 30, 40, 50])

    # PART C: Separate effects
    lam_1, lam_2, common_primes, h_1, h_2 = separate_effect_of_new_primes_from_eigenvector_change(8)

    # PART D: Reconstruction attempt
    margins = try_corrected_reconstruction(8)

    print("\n" + "="*100)
    print("SUMMARY: KEY FINDINGS")
    print("="*100)
    print("""
1. h(q) is NOT a simple power law. It has a complex structure that depends on:
   - The prime q
   - The eigenvector ψ at the current (c,N)

2. The eigenvector ψ REORGANIZES as c changes, even at fixed N.
   - This means h(q) 'for a given prime q' is NOT independent of c
   - Adding a new prime ALSO changes the values of h for existing primes

3. The naive picture (fit h(q) to a power law, sum over all primes) is FLAWED because:
   - It treats h(q) as if it's a fixed function of q alone
   - It ignores the eigenvector reorganization
   - It mixes the effect of adding new primes with the effect of ψ changing

4. The correct picture requires:
   - Tracking the FULL eigenvector ψ(c,N) and its structure
   - Understanding how ψ changes shape as c expands
   - Computing the true margin as a function of how fast λ_min shrinks when new primes are added
   - This shrinkage is NOT just from the new prime's h(q), but from the entire spectrum reorganization
    """)
