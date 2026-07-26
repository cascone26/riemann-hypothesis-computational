#!/usr/bin/env python3
"""
BEURLING CASCADE MEASUREMENT - Test whether cascade halts for perturbed primes

Uses fast reduced-precision computation (dps=50, N=15, T=100) to test the hypothesis.

Command to run: python3 beurling_cascade_measure.py
"""

import sys
sys.path.insert(0, "/Users/scones/projects/riemann/.venv_cvs/lib/python3.14/site-packages")

import mpmath as mp
import connes_cvs.operator as O
import json
import time

def get_real_primes(limit):
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit+1, i):
                is_prime[j] = False
    return [i for i in range(2, limit+1) if is_prime[i]]

def build_beurling_prime_data(max_val=13, perturb_mode="sine", strength=0.1):
    """
    Build custom prime-power data for Beurling primes to pass to psi_prime.

    Returns: list of (prime_power_index, log_prime, multiplicity) tuples
    """
    real_primes = get_real_primes(max_val)
    N = len(real_primes)
    beurling_primes = []

    for k, p in enumerate(real_primes, 1):
        if perturb_mode == "sine":
            delta = strength * mp.sin(mp.mpf(k) * mp.pi / mp.mpf(N))
        elif perturb_mode == "none":
            delta = 0
        elif perturb_mode == "scale":
            delta = strength * mp.mpf(p)
        else:  # cosine
            delta = strength * mp.cos(mp.mpf(k) * mp.pi / mp.mpf(N))

        p_b = mp.mpf(p) + delta
        beurling_primes.append(p_b)

    # Build prime-power data
    pdata = []
    for q in beurling_primes:
        log_q = mp.log(q)
        for k in range(1, 3):  # up to q^2
            power_val = q ** k
            pdata.append((float(power_val), log_q, k))

    return beurling_primes, pdata

def run_cascade_for_primes(pdata_list, c=13, N=15, dps=50, T=100, description=""):
    """
    Run cascade computation for a given set of prime-power data.

    pdata_list: list of (prime_power_val, log_prime, multiplicity) for primes
    """
    mp.mp.dps = dps
    L = mp.log(mp.mpf(c))
    idx = list(range(-N, N+1))
    DIM = 2*N + 1

    print(f"    Computing cascade for {description}...", flush=True)
    start_time = time.time()

    # Arch/pole
    arcpol = {}
    arcpol_d = {}
    for n in idx:
        x = mp.mpf(n)
        O._hplus_cache_clear()
        arcpol[n] = O.psi_pole(x, L) + O.psi_arch(x, L, T, dps)
        arcpol_d[n] = O.psi_pole_deriv(x, L) + O.psi_arch_deriv(x, L, T, dps)
        O._hplus_cache_clear()

    # Per-prime contributions
    pp = {}
    for (ppval, logp, mult) in pdata_list:
        val = {}
        der = {}
        # Note: we're using the Connes library function with our custom pdata
        for n in idx:
            val[n] = O.psi_prime(mp.mpf(n), L, [(ppval, logp, mult)])
            der[n] = O.psi_prime_deriv(mp.mpf(n), L, [(ppval, logp, mult)])
        pp[ppval] = (val, der)

    # Build matrix
    psi_full = {}
    psid_full = {}
    for n in idx:
        pv = sum((pp[ppval][0][n] for ppval, _, _ in pdata_list), mp.mpf(0))
        pd = sum((pp[ppval][1][n] for ppval, _, _ in pdata_list), mp.mpf(0))
        psi_full[n] = pv + arcpol[n]
        psid_full[n] = pd + arcpol_d[n]

    Q = mp.matrix(DIM, DIM)
    for i in range(DIM):
        m = i - N
        for j in range(DIM):
            nn = j - N
            Q[i, j] = psid_full[nn] if m == nn else (psi_full[m] - psi_full[nn]) / (m - nn)

    for i in range(DIM):
        for j in range(i+1, DIM):
            a = (Q[i, j] + Q[j, i]) / 2
            Q[i, j] = a
            Q[j, i] = a

    # Even-sector projection
    V = mp.matrix(DIM, N+1)
    V[N, 0] = mp.mpf(1)
    r2 = 1 / mp.sqrt(2)
    for k in range(1, N+1):
        V[N+k, k] = r2
        V[N-k, k] = r2

    Qe = V.T * Q * V
    eigs, vecs = mp.eigsy(Qe)
    k0 = min(range(N+1), key=lambda i: eigs[i])
    psi_eigvec = mp.matrix([vecs[r, k0] for r in range(N+1)])

    # Cascade measurement
    def assemble_from(valmap, dermap):
        Q = mp.matrix(DIM, DIM)
        for i in range(DIM):
            m = i - N
            for j in range(DIM):
                nn = j - N
                Q[i, j] = dermap[nn] if m == nn else (valmap[m] - valmap[nn]) / (m - nn)
        for i in range(DIM):
            for j in range(i+1, DIM):
                a = (Q[i, j] + Q[j, i]) / 2
                Q[i, j] = a
                Q[j, i] = a
        return V.T * Q * V

    Qarch = assemble_from(arcpol, arcpol_d)
    residual = psi_eigvec.T * (Qarch * psi_eigvec)
    residual_val = residual[0]
    cascade_log10 = [float(mp.log10(abs(residual_val))) if residual_val != 0 else 0]

    # Accumulate per-prime
    for ppval, _, _ in pdata_list[:10]:  # First 10
        Vn = assemble_from(pp[ppval][0], pp[ppval][1])
        hn = psi_eigvec.T * (Vn * psi_eigvec)
        residual_val = residual_val + hn[0]
        cascade_log10.append(float(mp.log10(abs(residual_val))) if residual_val != 0 else 0)

    elapsed = time.time() - start_time

    # Analysis
    orders_dropped = cascade_log10[0] - cascade_log10[-1]
    halts_at = None
    for i in range(1, len(cascade_log10)):
        if cascade_log10[i] >= cascade_log10[i-1] - 0.3:
            halts_at = i
            break

    return {
        'cascade_steps': cascade_log10,
        'orders_dropped': float(orders_dropped),
        'halts_at': halts_at,
        'num_terms': len(cascade_log10) - 1,
        'elapsed_seconds': elapsed
    }

def main():
    print("="*80)
    print("BEURLING CASCADE MEASUREMENT")
    print("Testing: Does the cascade halt for Beurling generalized-number systems?")
    print("="*80)

    # Get real prime data (baseline)
    print("\n1. BASELINE: Real primes (ℤ)")
    print("-"*80)
    pdata_real, _ = O.prime_powers_up_to(13)
    baseline = run_cascade_for_primes(pdata_real, c=13, N=15, dps=50, T=100, description="ℤ (real primes)")
    print(f"    Orders dropped: {baseline['orders_dropped']:.3f}")
    print(f"    Halts at step: {baseline['halts_at']}")
    print(f"    Cascade pattern: {baseline['cascade_steps'][:6]}")

    # Test Beurling systems
    tests = [
        ("sine_weak", "sine", 0.05),
        ("sine_moderate", "sine", 0.1),
        ("sine_strong", "sine", 0.2),
        ("scale_110", "scale", 0.1),
        ("cosine_weak", "cosine", 0.05),
    ]

    results = {"baseline_real_primes": baseline}

    print("\n2. BEURLING SYSTEMS: Perturbed primes")
    print("-"*80)

    for test_name, perturb_mode, strength in tests:
        print(f"\n  Test: {test_name} (mode={perturb_mode}, strength={strength})")
        beurling_primes, pdata_beurling = build_beurling_prime_data(
            max_val=13,
            perturb_mode=perturb_mode,
            strength=strength
        )
        print(f"    Beurling primes (first 3): {[float(mp.nstr(p, 5)) for p in beurling_primes[:3]]}")

        try:
            result = run_cascade_for_primes(
                pdata_beurling,
                c=13,
                N=15,
                dps=50,
                T=100,
                description=f"Beurling ({test_name})"
            )
            print(f"    Orders dropped: {result['orders_dropped']:.3f}")
            print(f"    Halts at step: {result['halts_at']}")
            print(f"    Cascade pattern: {result['cascade_steps'][:6]}")
            print(f"    Ratio vs baseline: {result['orders_dropped'] / baseline['orders_dropped']:.3f}")

            results[test_name] = result
        except Exception as e:
            print(f"    ERROR: {e}")
            results[test_name] = {"error": str(e)}

    # Summary
    print("\n" + "="*80)
    print("SUMMARY & INTERPRETATION")
    print("="*80)

    print("\nRESULTS:")
    print(f"{'Test Name':<20} {'Orders Dropped':<18} {'Halts At':<15} {'vs Baseline':<15}")
    print("-"*68)

    for name, result in results.items():
        if "error" not in result:
            halts_str = str(result['halts_at']) if result['halts_at'] else "Never"
            ratio = result['orders_dropped'] / baseline['orders_dropped'] if baseline else 1.0
            print(f"{name:<20} {result['orders_dropped']:>10.3f} orders  {halts_str:>12}     {ratio:>6.3f}x")
        else:
            print(f"{name:<20} ERROR")

    print("\nINTERPRETATION:")
    print("-"*80)
    print("Hypothesis: Cascade HALTS for Beurling systems (not ℤ-specific)")
    print("Null hypothesis: Cascade SURVIVES for Beurling (robust, not ℤ-specific)")

    beurling_results = {k: v for k, v in results.items() if k != "baseline_real_primes" and "error" not in v}

    if beurling_results:
        avg_ratio = sum(r['orders_dropped'] for r in beurling_results.values()) / len(beurling_results) / baseline['orders_dropped']
        avg_halt = sum(r['halts_at'] if r['halts_at'] else 999 for r in beurling_results.values()) / len(beurling_results)

        print(f"\nBeurling cascade behavior:")
        print(f"  Average orders dropped (ratio to baseline): {avg_ratio:.3f}")
        print(f"  Average halt step: {avg_halt:.1f}")

        if avg_ratio > 0.8 and avg_halt > 5:
            print(f"\n  VERDICT: Cascade is ROBUST to Beurling perturbations")
            print(f"           ⟹ Cascading is NOT ℤ-specific (Poisson-independent)")
            print(f"           ⟹ Cascading mechanism is more general than Poisson summation")
        elif avg_halt < 3:
            print(f"\n  VERDICT: Cascade HALTS early for Beurling primes")
            print(f"           ⟹ Cascading IS ℤ-specific (Poisson-dependent)")
            print(f"           ⟹ Hypothesis CONFIRMED: cascade breaks without integer structure")
        else:
            print(f"\n  VERDICT: Inconclusive - cascade partially degrades")
            print(f"           ⟹ Cascading depends on integer structure but not completely")

    # Save results
    output_file = "/Users/scones/projects/riemann/operators/track_crossblock/beurling_cascade_measurement.json"
    with open(output_file, 'w') as f:
        # Clean up cascade_steps for JSON
        clean_results = {}
        for k, v in results.items():
            if "error" not in v:
                clean_results[k] = {
                    'orders_dropped': v['orders_dropped'],
                    'halts_at': v['halts_at'],
                    'num_terms': v['num_terms'],
                    'elapsed_seconds': v['elapsed_seconds']
                }
            else:
                clean_results[k] = v

        json.dump({"results": clean_results}, f, indent=2)

    print(f"\n{'='*80}")
    print(f"Results saved to: {output_file}")

if __name__ == "__main__":
    main()
