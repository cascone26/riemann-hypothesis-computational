#!/usr/bin/env python3
"""
ANALYTICAL DERIVATION OF THE -0.44*N*log(c) EXPONENT
======================================================

Goal: Starting from the parity-sector cascade mechanism, derive why
log(λ_min) ~ -0.44*N*log(c) analytically, not just empirically.

Key inputs:
1. Parity decomposition: τ splits into independent even/odd sectors
2. Each sector has dimension ~ N (even) and ~ N (odd)
3. W02-WR has negative spectrum in each sector
4. Prime powers provide cancellation via Hellmann-Feynman cascade
5. Each prime power q contributes h(q) proportional to log(p)/sqrt(q)

The cascade is power-law: log|residual_k| ~ -β*log(k)
where β ≈ 11.6, which suggests structured near-nullness.

Hypothesis: The margin is determined by the "cancellation bandwidth" of the primes
relative to the "negativity span" across the frequency dimension.

Step 1: Estimate the negative spectrum of W02-WR in one parity sector
Step 2: Estimate the effective rank and spectral density of Wp_sector
Step 3: Model the cascade telescoping as power-law decay
Step 4: Compute the total accumulated decay and predict log(λ_min)
Step 5: Extract the -0.44*N exponent analytically
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special, linalg
import math
from math import log, exp, pi

# Empirical findings to anchor against:
EMPIRICAL_EXPONENT = 0.44  # λ_min ~ c^(-0.44*N)
EMPIRICAL_SLOPE_vs_logc = -3.428  # at N=8
EMPIRICAL_CASCADE_BETA = 11.6  # log|residual| ~ -11.6*log(prime_index) at c=30,N=8

def prime_counting_pi(c):
    """Approximate number of primes ≤ c via prime number theorem."""
    if c < 2:
        return 0
    return c / (log(c) - 1)

def prime_powers_count_M(c):
    """Approximate number of prime powers ≤ c.

    M = sum_{p prime, p<=c} floor(log_p(c)) ≈ π(c) * avg_power
    ≈ π(c) * log(log(c)) / log(2) for rough estimate
    Actually, M ~ π(c) for small c, more precisely:
    M = π(c) + π(c^(1/2)) + π(c^(1/3)) + ...
    For practical c (30-1000), this is dominated by π(c).
    """
    M = 0
    # Primes and their powers
    for q in range(2, int(c) + 1):
        # Check if q is a prime power
        for p in range(2, int(q**0.5) + 2):
            if q == p:
                M += 1
                break
            # Check powers
            pk = p
            while pk <= q:
                if pk == q:
                    M += 1
                    break
                pk *= p
    return M

def negative_eigenvalue_count_per_sector(c, N):
    """
    Estimate the number of negative eigenvalues in (W02-WR)_sector.

    From empirical observation (S2627):
    - (W02-WR)_even: ~2 negative eigenvalues
    - (W02-WR)_odd: ~4 negative eigenvalues

    Hypothesis: related to low-frequency modes coupling to the archimedean block.
    Very roughly, we expect ~ log(N) or ~ sqrt(N) negative eigenvalues per sector,
    but empirically it's small and roughly constant.

    For now, use empirical values but scale with N.
    """
    # Observed: at N=8, even has ~2 neg, odd has ~4 neg
    # Assume this scales roughly with ln(N) or stays constant
    n_neg_even = max(2, int(math.ceil(log(N + 1))))
    n_neg_odd = max(3, int(math.ceil(1.5 * log(N + 1))))
    return n_neg_even, n_neg_odd

def spectral_density_and_norm(c, N):
    """
    Estimate the spectral norm and sparsity of W02-WR and Wp in one sector.

    W02-WR: Toeplitz-like, bounded eigenvalues
    - Archimedean part: ~O(1) in scale
    - Regular part: similar scale
    - Spectral norm σ_max ~ O(1)

    Wp: weighted sum of rank-2 (or low-rank) matrices per prime
    - Each matrix has weight ~ log(p)/sqrt(q)
    - Summed over M primes
    - σ_max(Wp) ~ sum of weights ~ ∫_2^c log(p)/sqrt(p) dp ~ O(log(c))

    Intuitively: σ_max(Wp_sector) ~ O(log(c))
    """
    # Spectral norms (rough, from theory + numerics)
    sigma_max_arch = 2.0  # W02-WR is bounded, σ_max ~ O(1) from Arb computations

    # Wp: sum of log(p)/sqrt(q) over prime powers
    # ∫_2^c log(p)/sqrt(p) dp ~ O(log(c)) asymptotically
    sigma_max_wp = 0.5 * log(c) + 1.0  # empirical scaling

    return sigma_max_arch, sigma_max_wp

def cascade_decay_rate(c, N, M):
    """
    Model the per-prime Hellmann-Feynman cascade decay.

    Empirical finding: log|residual_k| ~ -β*log(k) where β ≈ 11.6 (at c=30,N=8)

    Hypothesis: β depends on the "effective dimension" and "cancellation efficiency".

    β ~ (number of negative evals to cancel) * (efficiency factor)
      ~ n_neg * efficiency

    If efficiency ~ log(N) (from Fourier concentration), then:
    β ~ n_neg * log(N)

    At c=30, N=8: n_neg ~ 2-4, log(N) ~ 2, so β ~ 4-8... but we observe 11.6.

    Alternative: β ~ M / sqrt(N) or similar.
    M/sqrt(N) ~ 10/sqrt(8) ~ 3.5... still too small.

    Let me just use the empirical β for now and model its scaling.
    """
    # Empirical: β ≈ 11.6 at c=30, N=8
    # Hypothesize: β ~ constant * M / N or β ~ constant * log(c) * N

    # Simple model: β proportional to M and log(N)
    # β ≈ 11.6 (observed), try to match
    # If β ~ a * (c/ln(c)) * ln(N), then at c=30,N=8:
    # β ~ a * (30/ln(30)) * ln(8) ~ a * 11 * 2 ~ a * 22, so a ~ 0.5

    beta_empirical = EMPIRICAL_CASCADE_BETA  # 11.6

    # Model: β(c,N) ~ β_0 * f(c) * g(N)
    # From the data: β is roughly independent of c? or grows with ln(c)?
    # Let's assume β ~ ln(c) * N for a test
    beta_model = log(c + 1) * N * 0.3  # empirical tuning

    return beta_empirical, beta_model

def total_cascade_decay(c, N, M, beta):
    """
    Model the total decay from the cascade.

    If log|residual_k| ~ -β*log(k) for k=1..M, then:
    log|residual_final| ~ -β * log(M!) ~ -β * (M*log(M) - M)

    At c=30, N=8: M≈10, β≈11.6
    Total decay ~ -11.6 * (10*ln(10) - 10) ~ -11.6 * (23 - 10) ~ -151

    So |residual| ~ exp(-151) which is tiny.

    But we observe λ_min ~ 10^(-26.73) = exp(-61.5*ln(10)) ≈ exp(-141) in natural log.

    Actually wait, let me recalculate. If log|residual_k| ~ -β*log(k), this is log_e.

    sum_{k=1}^M (-β*log(k)) = -β * log(M!) ≈ -β * M * log(M)

    At c=30, N=8: -11.6 * 10 * log(10) ≈ -11.6 * 10 * 2.3 ≈ -267

    But the true log10(λ_min) is only -26.73, which is ln(10) * 26.73 ≈ 61.5.

    There's a factor-of-4 discrepancy. This suggests either:
    1. The cascade decay rate β is not constant but changes with k
    2. The logarithm base is different (log10 vs ln)
    3. Only a subset of M primes contribute effectively
    """
    # In natural log
    M_fact = M * log(M + 0.1) - M + 0.5  # Stirling's approximation for log(M!)
    total_decay_ln = -beta * M_fact

    # Convert to log10
    total_decay_log10 = total_decay_ln / log(10)

    return total_decay_log10

def model_margin_scaling(c, N):
    """
    Comprehensive model of the margin.

    Components:
    1. Background negativity: (W02-WR)_sector has eigenvalue ~ -O(1)
    2. Number of primes: M ~ c/log(c)
    3. Parity sectors: even and odd decouple
    4. Cascade decay: each prime contributes via Hellmann-Feynman
    5. Final margin: balance between negativity and cancellation

    Predicted form: log10(λ_min) ~ -α*N*log(c) + β(N)
    where α ≈ 0.44 is to be derived.
    """
    print(f"\n{'='*100}")
    print(f"ANALYTICAL MARGIN DERIVATION: c={c}, N={N}")
    print(f"{'='*100}")

    M = prime_powers_count_M(c)
    M_pnt = prime_counting_pi(c)

    print(f"\nStructural parameters:")
    print(f"  c (prime cutoff): {c}")
    print(f"  N (frequency cutoff): {N}")
    print(f"  dim per sector: {N+1} (even) and {N} (odd)")
    print(f"  M (# prime powers ≤ c): {M}")
    print(f"  π(c) (# primes ≤ c): {M_pnt:.1f}")
    print(f"  log(c): {log(c):.4f}")
    print(f"  M/N: {M/N:.4f}")
    print(f"  log(c)*N: {log(c)*N:.4f}")

    # Spectral estimates
    sigma_arch, sigma_wp = spectral_density_and_norm(c, N)
    print(f"\nSpectral scales:")
    print(f"  σ_max(W02-WR) ≈ {sigma_arch:.2f}")
    print(f"  σ_max(Wp) ≈ {sigma_wp:.2f}")
    print(f"  σ_max(Wp) / σ_max(W02-WR) ≈ {sigma_wp/sigma_arch:.2f}")

    # Negative eigenvalue count
    n_neg_even, n_neg_odd = negative_eigenvalue_count_per_sector(c, N)
    print(f"\nNegative eigenvalue estimates:")
    print(f"  (W02-WR)_even: ~{n_neg_even} negative evals")
    print(f"  (W02-WR)_odd: ~{n_neg_odd} negative evals")
    print(f"  Total negativity to cancel: ~{n_neg_even + n_neg_odd} evals")

    # Cascade decay
    beta_emp, beta_model = cascade_decay_rate(c, N, M)
    print(f"\nCascade decay rate:")
    print(f"  Empirical β: {beta_emp:.2f}")
    print(f"  Model β(c,N): {beta_model:.2f}")

    # Total decay from cascade
    decay_log10_emp = total_cascade_decay(c, N, M, beta_emp)
    decay_log10_model = total_cascade_decay(c, N, M, beta_model)
    print(f"\nTotal cascade decay (log10):")
    print(f"  Via empirical β: {decay_log10_emp:.2f}")
    print(f"  Via model β: {decay_log10_model:.2f}")

    # Predicted margin
    print(f"\nPredicted margin (first-order):")
    print(f"  If M primes create ~decay_log10 total decay, then λ_min ~ 10^({decay_log10_emp:.1f})")

    # Empirical reference
    empirical_log10_margin = -EMPIRICAL_SLOPE_vs_logc * log(c) - 14.5  # From margin_law.py at N=8 for one point
    # Better: use the full scaling law
    empirical_log10_margin_general = -0.44 * N * log(c)

    print(f"\nEmpirical reference (λ_min ~ c^(-0.44*N)):")
    print(f"  log10(λ_min) ~ -0.44*{N}*log({c}) = {empirical_log10_margin_general:.2f}")
    print(f"  (At N=8: -0.44*8*log(c) = -3.52*log(c), observed -3.428*log(c))")

    # Analysis of the 0.44 exponent
    print(f"\n{'='*100}")
    print(f"EXPONENT ANALYSIS:")
    print(f"{'='*100}")

    print(f"\nHypothesis 1: Exponent from prime density")
    print(f"  If M ~ c/ln(c) and margin ~ 1/M, then:")
    print(f"    log(margin) ~ -log(M) ~ -log(c/ln(c)) ~ -log(c)")
    print(f"  This gives exponent α=1, not 0.44 ✗")

    print(f"\nHypothesis 2: Exponent from dimension scaling")
    print(f"  If margin ~ 1/(M/N)^p = (N/M)^p, then:")
    print(f"    log(margin) ~ p*log(N/M) ~ p*log(N/(c/ln(c))) ~ p*[log(N*ln(c)/c)]")
    print(f"  For p=0.44: log(margin) ~ 0.44*[log(N) + log(ln(c)) - log(c)]")
    print(f"  At c=30, N=8: -0.44*[log(8) - log(30)] = -0.44*[0.90 - 1.40] = +0.22")
    print(f"  This doesn't match either ✗")

    print(f"\nHypothesis 3: Exponent from eigenvalue count vs frequency resolution")
    print(f"  # frequencies = 2N+1")
    print(f"  # negative evals ~ O(1) or O(log N)")
    print(f"  If each negative eval needs M/log(N) primes to cancel (logarithmic effectiveness),")
    print(f"  then effective_M ~ M/log(N)")
    print(f"  margin ~ (1/effective_M)^{0.44*N} ~ [(ln(N)/M)^{0.44*N}]")
    print(f"  This is getting closer but still complex ...")

    print(f"\nHypothesis 4: Exponent from cascade decay structure (MOST LIKELY)")
    print(f"  The cascade decay is POWER-LAW: log|residual| ~ -β*log(k)")
    print(f"  where β ~ (some scale) * N")
    print(f"  ")
    print(f"  When summed over all M primes:")
    print(f"    log(final residual) ~ -β*log(M!) ~ -β*M*log(M)")
    print(f"  ")
    print(f"  If β scales as β ~ c^α/N (or similar), then:")
    print(f"    log(residual) ~ -(c^α/N) * (c/ln(c)) * ln(c/ln(c))")
    print(f"                  ~ -c^(1+α) * ln(c) / [N*ln(c)]")
    print(f"                  ~ -c^(1+α) / N")
    print(f"  ")
    print(f"  For this to match c^(-0.44*N), we'd need a different dependence...")
    print(f"  ")
    print(f"  REVISED: If each prime k contributes -β*log(k) and β itself depends on k:")
    print(f"    β(k) ~ f(k) where f is to be determined")
    print(f"  Then the structure becomes much richer.")

    print(f"\nHypothesis 5: Per-prime contribution scales with prime's 'reach' in Fourier space")
    print(f"  Prime p creates oscillations at frequency 2π*p/L")
    print(f"  In frequency domain n ∈ [-N,N], this reaches up to ~N frequencies")
    print(f"  So each prime 'covers' ~ min(p, N) frequencies")
    print(f"  ")
    print(f"  Effective rank of Wp ~ (# primes covering distinct frequencies)")
    print(f"  For large c and fixed N, this plateaus at ~N (all frequencies covered many times)")
    print(f"  ")
    print(f"  If effective_rank ~ min(M, N^2), then:")
    print(f"    margin ~ 1/(effective_rank)^0.44 ~ (min(M,N^2))^(-0.44)")
    print(f"  ")
    print(f"  When M ~ c/ln(c) >> N^2 (which is true for N~8, c~30):")
    print(f"    margin ~ M^(-0.44) ~ (c/ln(c))^(-0.44) ~ c^(-0.44)")
    print(f"  But this ignores the N-dependence!")

    print(f"\n" + "="*100)
    print(f"CRITICAL INSIGHT:")
    print(f"="*100)
    print(f"\nThe empirical law log10(λ) ~ -0.44*N*log(c) has BOTH c and N in the exponent.")
    print(f"This suggests the exponent is a PRODUCT of two independent scalings:")
    print(f"  • One dimension from N (frequency resolution, # frequencies to worry about)")
    print(f"  • One dimension from log(c) (# primes available to do the work)")
    print(f"")
    print(f"PRODUCT FORM: If margin ~ (f(c))^(g(N)), then:")
    print(f"  log(margin) ~ g(N)*log(f(c))")
    print(f"  = (something depending on N) * (something depending on log(c))")
    print(f"")
    print(f"With λ_min ~ c^(-0.44*N):")
    print(f"  log(λ_min) = -0.44*N*log(c)")
    print(f"  So g(N) = 0.44*N  and log(f(c)) = log(c), meaning f(c) = c")
    print(f"")
    print(f"This means: margin ~ c^(-0.44*N)")
    print(f"")
    print(f"The 0.44 coefficient arises because:")
    print(f"  • Negative eigenvalue count in (W02-WR)_sector scales with N")
    print(f"    (more frequencies → more potential for mode coupling → more negatives)")
    print(f"  • Prime effectiveness scales with c (more primes = more cancellation, but with diminishing returns)")
    print(f"  • The PRODUCT of these two effects gives margin ~ c^(-const*N)")
    print(f"")
    print(f"The 0.44 specifically may come from:")
    print(f"  • Empirical trade-off between negativity growth (N) and cancellation capacity (ln(c))")
    print(f"  • Or from fundamental Fourier analysis: if M ~ c/ln(c) and dim ~ N,")
    print(f"    then M/N ~ c/(N*ln(c)), and log(margin) ~ -(M/N)*log(something) ~ -c/(N*ln(c))*log(...)")
    print(f"  • Exact coefficient requires computation of the Fourier basis decay rates")


def test_against_empirical():
    """Test the model predictions against empirical data from margin_law.py"""
    print(f"\n\n{'='*100}")
    print(f"MODEL vs EMPIRICAL COMPARISON")
    print(f"{'='*100}\n")

    test_points = [
        (30, 4),
        (30, 6),
        (30, 8),
        (30, 10),
        (13, 8),
        (100, 8),
    ]

    for c, N in test_points:
        model_margin_scaling(c, N)


if __name__ == '__main__':
    # Test the analytical model
    test_against_empirical()
