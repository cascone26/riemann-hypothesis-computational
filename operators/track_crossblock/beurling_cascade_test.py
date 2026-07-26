#!/usr/bin/env python3
"""
BEURLING CASCADE TEST — L2 Probe §214

Hypothesis: the 58-order telescoping cascade in the ℤ case relies on Poisson
summation over the additive structure of integers. If true, the cascade should
HALT or stall for a Beurling generalized-number system where the primes are
perturbed to break Poisson symmetry.

Method:
1. Build a Beurling system with perturbed primes: p'_k = p_k + ε·sin(k·π/10)
2. Adapt the cascade computation using Beurling "von Mangoldt" weights
3. Accumulate residual contributions from Beurling primes (analogue of h(n))
4. Measure: does the residual telescope to 50+ orders (like ℤ), or halt early?

Output: raw residual values, order-of-magnitude drops, and comparison to ℤ case.
"""
import sys
import mpmath as mp
from fractions import Fraction
import json

def get_real_primes(limit):
    """Return list of primes up to limit."""
    def sieve(n):
        is_prime = [True] * (n + 1)
        is_prime[0] = is_prime[1] = False
        for i in range(2, int(n**0.5) + 1):
            if is_prime[i]:
                for j in range(i*i, n+1, i):
                    is_prime[j] = False
        return [i for i in range(2, n+1) if is_prime[i]]
    return sieve(limit)

def build_beurling_system(perturb_mode="sine", perturb_strength=0.1, max_val=100):
    """
    Build a Beurling generalized-number system.

    perturb_mode: "sine", "scale", "random"
    - "sine": p'_k = p_k + ε·sin(k·π/N) where N=π(max_val)
    - "scale": p'_k = α·p_k (density scaling)
    - "random": p'_k = p_k·(1 + ε·pseudo-random), breaks Poisson

    Returns: {
        "real_primes": [...],
        "beurling_primes": [...],
        "perturbation": {...},
        "description": str
    }
    """
    real_primes = get_real_primes(max_val)
    N = len(real_primes)

    if perturb_mode == "sine":
        # Perturb by sine wave — breaks Poisson summing structure
        beurling_primes = []
        perturbations = {}
        for k, p in enumerate(real_primes, 1):
            delta = perturb_strength * mp.sin(mp.mpf(k) * mp.pi / mp.mpf(N))
            p_beurling = mp.mpf(p) + delta
            beurling_primes.append(p_beurling)
            perturbations[p] = float(delta)
        desc = f"Sine perturb: p'_k = p_k + {perturb_strength}·sin(k·π/{N})"

    elif perturb_mode == "scale":
        # Scale all primes by α ≠ 1 — changes density, should break Poisson
        alpha = 1.0 + perturb_strength  # e.g., 1.1 for +10% scaling
        beurling_primes = [mp.mpf(p) * mp.mpf(alpha) for p in real_primes]
        perturbations = {p: float(mp.mpf(p) * (mp.mpf(alpha) - 1)) for p in real_primes}
        desc = f"Scale perturb: p'_k = {alpha:.3f}·p_k (density +{(alpha-1)*100:.1f}%)"

    else:  # random
        # Pseudo-random perturbation (seeded for reproducibility)
        import hashlib
        beurling_primes = []
        perturbations = {}
        for p in real_primes:
            seed = int(hashlib.md5(str(p).encode()).hexdigest(), 16) % 1000000
            mp.seed(seed)
            delta = perturb_strength * mp.mpf(p) * (mp.mpf(2) * mp.rand() - mp.mpf(1))
            p_beurling = mp.mpf(p) + delta
            beurling_primes.append(p_beurling)
            perturbations[p] = float(delta)
        desc = f"Random perturb: p'_k = p_k·(1 ± {perturb_strength}) [seeded]"

    return {
        "real_primes": real_primes,
        "beurling_primes": beurling_primes,
        "perturbations": perturbations,
        "perturb_mode": perturb_mode,
        "perturb_strength": perturb_strength,
        "description": desc,
        "num_primes": len(real_primes)
    }

def build_beurling_von_mangoldt(beurling_system, max_prime_power=None):
    """
    Build a "von Mangoldt function" for Beurling primes.

    For standard case: Λ(n) = log(p) if n = p^k, else 0.
    For Beurling: Λ_B(q_k) = log(q_k) where q_k are Beurling primes.
    Also include prime powers: q_k^2, q_k^3, etc.

    Returns: {
        "indices": [q_k, q_k^2, ...],  # as Decimal-like objects
        "weights": [log(q_k), 2·log(q_k), 3·log(q_k), ...]
    }
    """
    bp = beurling_system["beurling_primes"]
    if max_prime_power is None:
        max_prime_power = max(5, len(bp) // 3)  # include prime powers up to p^k, k~3-5

    indices = []
    weights = []

    for i, q_k in enumerate(bp):
        log_q = mp.log(q_k)
        for k in range(1, max_prime_power + 1):
            power_val = q_k ** k
            indices.append(power_val)
            weights.append(k * log_q)

    return {
        "indices": indices,
        "weights": weights,
        "num_terms": len(indices),
        "max_prime_power": max_prime_power
    }

def simple_cascade_analogue(beurling_system, dps=100, num_terms=None):
    """
    Simple analogue of the cascade: sum up Beurling "von Mangoldt" contributions
    and measure how fast the residual decays.

    This is a toy version: we don't build the full Connes operator, but instead
    just accumulate contributions from Beurling primes in order, similar to how
    h(n) accumulates in the real cascade.

    Returns: {
        "residual_log10": [log10 values at each step],
        "indices_processed": [corresponding Beurling prime indices],
        "num_orders_dropped": int,
        "halts_after": int (first index where residual plateaus)
    }
    """
    mp.mp.dps = dps
    bp = beurling_system["beurling_primes"]

    if num_terms is None:
        num_terms = len(bp)

    # Initialize residual: start with an "arch" expectation (toy: log of largest Beurling prime)
    residual = mp.log(bp[-1])

    residual_log10_vals = [float(mp.log10(residual)) if residual != 0 else float('-inf')]
    indices_processed = []

    # Accumulate contributions from Beurling primes in order
    prev_residual_log10 = residual_log10_vals[0]
    plateau_count = 0
    halts_after = None

    for i, q_k in enumerate(bp[:num_terms]):
        # Toy contribution: log(q_k) / sqrt(index) (analogue of Λ(n)/√n)
        contribution = mp.log(q_k) / mp.sqrt(mp.mpf(i + 2))
        residual = residual - contribution

        if residual != 0:
            res_log10 = float(mp.log10(abs(residual)))
        else:
            res_log10 = float('-inf')

        residual_log10_vals.append(res_log10)
        indices_processed.append(float(q_k))

        # Check for plateau (residual stops decreasing significantly)
        if res_log10 >= prev_residual_log10 - 0.5:  # drop < 0.5 orders → stalling
            plateau_count += 1
            if plateau_count >= 2 and halts_after is None:
                halts_after = i
        else:
            plateau_count = 0

        prev_residual_log10 = res_log10

    num_orders_dropped = residual_log10_vals[0] - residual_log10_vals[-1]

    return {
        "residual_log10": residual_log10_vals,
        "indices_processed": indices_processed,
        "num_orders_dropped": float(num_orders_dropped),
        "halts_after": halts_after,
        "total_steps": len(residual_log10_vals) - 1,
        "final_residual_log10": residual_log10_vals[-1]
    }

def run_cascade_tests():
    """Run cascade tests for several Beurling systems and compare to ℤ."""
    print("="*80)
    print("BEURLING CASCADE TEST — Does the 58-order cascade survive Beurling perturbations?")
    print("="*80)

    tests = [
        ("sine_mild", "sine", 0.05),
        ("sine_moderate", "sine", 0.1),
        ("sine_strong", "sine", 0.2),
        ("scale_105", "scale", 0.05),
        ("scale_110", "scale", 0.10),
    ]

    results = {}

    for test_name, perturb_mode, perturb_strength in tests:
        print(f"\n{'-'*80}")
        print(f"Test: {test_name} ({perturb_mode}, strength={perturb_strength})")
        print(f"{'-'*80}")

        # Build Beurling system
        beurling = build_beurling_system(
            perturb_mode=perturb_mode,
            perturb_strength=perturb_strength,
            max_val=100
        )
        print(f"Description: {beurling['description']}")
        print(f"Number of primes: {beurling['num_primes']}")
        print(f"First 5 real primes: {beurling['real_primes'][:5]}")
        print(f"First 5 Beurling primes: {[float(x) for x in beurling['beurling_primes'][:5]]}")

        # Run cascade analogue
        cascade_result = simple_cascade_analogue(beurling, dps=100, num_terms=min(20, len(beurling['beurling_primes'])))
        print(f"\nCascade residual (log10 scale):")
        for i in range(0, len(cascade_result['residual_log10']), max(1, len(cascade_result['residual_log10']) // 10)):
            print(f"  Step {i:2d}: log10|R| = {cascade_result['residual_log10'][i]:8.3f}")

        print(f"\nSummary:")
        print(f"  Total steps: {cascade_result['total_steps']}")
        print(f"  Orders dropped: {cascade_result['num_orders_dropped']:.3f}")
        print(f"  Final residual (log10): {cascade_result['final_residual_log10']:.3f}")
        print(f"  Halts after (plateau): {cascade_result['halts_after']} steps")

        results[test_name] = {
            "beurling": beurling,
            "cascade": cascade_result
        }

    # Summary comparison
    print(f"\n{'='*80}")
    print("SUMMARY: Cascade behavior across all Beurling systems")
    print(f"{'='*80}")
    print(f"{'Test Name':<20} {'Perturbation':<15} {'Orders Dropped':<18} {'Halts After':<15}")
    print(f"{'-'*68}")
    for test_name, data in results.items():
        halts = data["cascade"]["halts_after"] if data["cascade"]["halts_after"] else "Never"
        print(f"{test_name:<20} {data['beurling']['description'][:14]:<15} {data['cascade']['num_orders_dropped']:>6.2f} orders    {str(halts):>12}")

    # Comparison to ℤ case
    print(f"\n{'='*80}")
    print("COMPARISON TO ℤ (Standard primes)")
    print(f"{'='*80}")
    print("ℤ cascade (c=13, c=12): ~58 orders dropped, telescopes ALL the way (no halt)")
    print("If Beurling systems show HALT at low orders → cascade IS ℤ-specific (Poisson-dependent)")
    print("If Beurling systems show 50+ orders → cascade is more general (not ℤ-specific)")

    # Save results
    output_file = "/Users/scones/projects/riemann/operators/track_crossblock/beurling_cascade_results.json"
    json.dump({
        "test_timestamp": str(mp.mpf(0)),  # placeholder
        "num_tests": len(results),
        "per_test_results": {
            name: {
                "description": data["beurling"]["description"],
                "num_orders_dropped": data["cascade"]["num_orders_dropped"],
                "halts_after": data["cascade"]["halts_after"],
                "total_steps": data["cascade"]["total_steps"],
                "final_residual_log10": data["cascade"]["final_residual_log10"]
            }
            for name, data in results.items()
        }
    }, open(output_file, 'w'), indent=2)
    print(f"\nResults saved to: {output_file}")

    return results

if __name__ == "__main__":
    run_cascade_tests()
