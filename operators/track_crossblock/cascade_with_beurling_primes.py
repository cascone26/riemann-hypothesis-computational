#!/usr/bin/env python3
"""
CASCADE WITH BEURLING PRIMES — Use real connes_cvs machinery but with Beurling primes

Strategy: The connes_cvs library accepts per-prime data as (prime_power_index, log_prime, weight).
Instead of using real primes, we'll:
1. Create Beurling primes (perturbed or scaled versions of real primes)
2. Build the prime_power data structure using Beurling primes
3. Feed it to psi_prime / psi_prime_deriv functions
4. Run the full cascade computation
5. Measure: does the cascade telescope to 50+ orders, or halt?

This is the MOST FAITHFUL test since it uses the actual connes_cvs library.
"""

import sys
sys.path.insert(0, "/Users/scones/projects/riemann/.venv_cvs/lib/python3.14/site-packages")

import mpmath as mp
import connes_cvs.operator as O
import time
import json

def get_real_primes(limit):
    """Sieve for primes up to limit."""
    is_prime = [True] * (limit + 1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, int(limit**0.5) + 1):
        if is_prime[i]:
            for j in range(i*i, limit+1, i):
                is_prime[j] = False
    return [i for i in range(2, limit+1) if is_prime[i]]

def build_beurling_primes(max_val=50, perturbation="sine", strength=0.1):
    """Build Beurling primes and return as (real_primes, beurling_primes)."""
    real_primes = get_real_primes(max_val)
    N = len(real_primes)
    beurling_primes = []

    for k, p in enumerate(real_primes, 1):
        if perturbation == "sine":
            delta = strength * mp.sin(mp.mpf(k) * mp.pi / mp.mpf(N))
        elif perturbation == "scale":
            delta = strength * mp.mpf(p)
        elif perturbation == "none":
            delta = 0
        else:  # cosine
            delta = strength * mp.cos(mp.mpf(k) * mp.pi / mp.mpf(N))

        beurling_primes.append(mp.mpf(p) + delta)

    return real_primes, beurling_primes

def build_prime_power_data(beurling_primes, max_power=3):
    """
    Build prime power data tuples for connes_cvs.psi_prime.

    Returns: list of (prime_power_val, log_prime, weight_index) tuples

    This mimics the output of O.prime_powers_up_to but using Beurling primes.
    """
    pdata = []
    for i, q in enumerate(beurling_primes):
        log_q = mp.log(q)
        for k in range(1, max_power + 1):
            power_val = q ** k
            pdata.append((float(power_val), log_q, k))
    return pdata

def run_cascade_with_beurling(max_val=50, perturbation="sine", strength=0.1, dps=100, N=20, T=200):
    """
    Run the cascade computation using Beurling primes with the real connes_cvs library.

    Returns: cascade_log10 values, num_orders_dropped, halt_info
    """
    mp.mp.dps = dps

    print(f"  Building Beurling primes: perturb={perturbation}, strength={strength}")
    real_primes, beurling_primes = build_beurling_primes(max_val, perturbation, strength)
    print(f"    Real primes (first 5): {real_primes[:5]}")
    print(f"    Beurling primes (first 5): {[float(mp.nstr(p, 6)) for p in beurling_primes[:5]]}")

    # Build the cascade computation
    # We'll use the logic from cvs_gate.py but with Beurling prime data

    c = max_val
    L = mp.log(mp.mpf(c))
    idx = list(range(-N, N+1))
    DIM = 2*N + 1

    print(f"  Computing arch/pole components (c={c}, N={N}, dps={dps}, T={T})...")
    t0 = time.time()

    # Compute archimedean components (same for any prime set — only depends on L)
    arcpol = {}
    arcpol_d = {}
    for n in idx:
        x = mp.mpf(n)
        O._hplus_cache_clear()
        arcpol[n] = O.psi_pole(x, L) + O.psi_arch(x, L, T, dps)
        arcpol_d[n] = O.psi_pole_deriv(x, L) + O.psi_arch_deriv(x, L, T, dps)
        O._hplus_cache_clear()

    t_arch = time.time() - t0
    print(f"    Arch components done in {t_arch:.1f}s")

    # Compute per-prime contributions using Beurling primes
    print(f"  Computing per-prime contributions for {len(beurling_primes)} Beurling primes...")
    t0 = time.time()

    pp = {}
    for i, q in enumerate(beurling_primes):
        log_q = mp.log(q)

        # Per-prime data tuple: (power_val, log_q, k) for k=1,2,3
        pdata = []
        for k in range(1, 4):  # prime powers up to q^3
            power_val = q ** k
            pdata.append((float(power_val), log_q, k))

        val = {}
        der = {}
        for n in idx:
            val[n] = O.psi_prime(mp.mpf(n), L, pdata)
            der[n] = O.psi_prime_deriv(mp.mpf(n), L, pdata)

        # Store using the Beurling prime as key (as a float for JSON-ability)
        pp[float(q)] = (val, der)

    t_prime = time.time() - t0
    print(f"    Per-prime components done in {t_prime:.1f}s")

    # Now assemble and compute eigenvector + cascade
    print(f"  Building operator matrix and extracting near-null eigenvector...")
    t0 = time.time()

    # Assemble full matrix (arch + all primes)
    psi_full = {}
    psid_full = {}
    for n in idx:
        pv = sum((pp[float(q)][0][n] for q in beurling_primes), mp.mpf(0))
        pd = sum((pp[float(q)][1][n] for q in beurling_primes), mp.mpf(0))
        psi_full[n] = pv + arcpol[n]
        psid_full[n] = pd + arcpol_d[n]

    # Build the difference-quotient matrix
    Q = mp.matrix(DIM, DIM)
    for i in range(DIM):
        m = i - N
        for j in range(DIM):
            nn = j - N
            Q[i, j] = psid_full[nn] if m == nn else (psi_full[m] - psi_full[nn]) / (m - nn)

    # Symmetrize
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

    # Eigenvalue decomposition
    eigs, vecs = mp.eigsy(Qe)
    k0 = min(range(N+1), key=lambda i: eigs[i])
    lam0 = eigs[k0]
    psi_eigvec = mp.matrix([vecs[r, k0] for r in range(N+1)])

    t_eig = time.time() - t0
    print(f"    Matrix + eigenvector done in {t_eig:.1f}s")
    print(f"    Near-null eigenvalue λ_min = {mp.nstr(lam0, 6)} (log10={float(mp.log10(abs(lam0))):.2f})")

    # Now compute the cascade: per-prime-power residual
    print(f"  Computing cascade residual...")

    # Assemble matrices for each Beurling prime
    def assemble_from(valmap, dermap):
        """Assemble difference-quotient matrix from value/deriv maps."""
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

    # Arch contribution
    Qarch = assemble_from(arcpol, arcpol_d)
    h_arch = psi_eigvec.T * (Qarch * psi_eigvec)
    h_arch_val = h_arch[0]

    cascade_log10 = []
    residual = h_arch_val

    if residual != 0:
        log10_val = float(mp.log10(abs(residual)))
    else:
        log10_val = float('-inf')

    cascade_log10.append(log10_val)
    print(f"    arch: log10|R| = {log10_val:8.3f}")

    # Per-prime contributions
    for i, q in enumerate(beurling_primes):
        Vn = assemble_from(pp[float(q)][0], pp[float(q)][1])
        hn = psi_eigvec.T * (Vn * psi_eigvec)
        hn_val = hn[0]
        residual = residual + hn_val

        if residual != 0:
            log10_val = float(mp.log10(abs(residual)))
        else:
            log10_val = float('-inf')

        cascade_log10.append(log10_val)
        print(f"    +q_{i}: log10|R| = {log10_val:8.3f}  (|h_n|=10^{float(mp.log10(abs(hn_val))) if hn_val!=0 else 0:.2f})")

    # Analysis
    num_orders = cascade_log10[0] - cascade_log10[-1]
    halts_at = None
    for i in range(1, len(cascade_log10)-1):
        if cascade_log10[i] >= cascade_log10[i-1] - 0.5:  # plateau
            halts_at = i
            break

    return {
        "cascade_log10": cascade_log10,
        "num_orders_dropped": float(num_orders),
        "halts_at": halts_at,
        "lambda_min": float(mp.nstr(lam0, 10)),
        "num_steps": len(cascade_log10) - 1,
        "final_residual": cascade_log10[-1] if cascade_log10 else None,
        "time_arch": t_arch,
        "time_prime": t_prime,
        "time_eig": t_eig
    }

def main():
    print("="*80)
    print("CASCADE WITH BEURLING PRIMES — Real connes_cvs machinery")
    print("="*80)

    # Test configurations
    tests = [
        ("regular_int", "none", 0.0),
        ("sine_weak", "sine", 0.05),
        ("sine_moderate", "sine", 0.1),
        ("scale_105", "scale", 0.05),
    ]

    results = {}

    for test_name, perturb, strength in tests:
        print(f"\n{'-'*80}")
        print(f"Test: {test_name}")
        print(f"{'-'*80}")

        try:
            result = run_cascade_with_beurling(
                max_val=50,
                perturbation=perturb,
                strength=strength,
                dps=100,
                N=20,
                T=200
            )
            results[test_name] = result

            print(f"\nSummary:")
            print(f"  Orders dropped: {result['num_orders_dropped']:.3f}")
            print(f"  Final residual (log10): {result['final_residual']:.3f}")
            print(f"  Halts after: {result['halts_at']} steps")
            print(f"  Total time: {result['time_arch'] + result['time_prime'] + result['time_eig']:.1f}s")

        except Exception as e:
            print(f"ERROR in {test_name}: {e}", flush=True)
            import traceback
            traceback.print_exc()
            results[test_name] = {"error": str(e)}

    # Summary
    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")

    for name, result in results.items():
        if "error" not in result:
            print(f"{name:<20}: {result['num_orders_dropped']:>6.2f} orders, halts_at={result['halts_at']}")
        else:
            print(f"{name:<20}: ERROR - {result['error']}")

    return results

if __name__ == "__main__":
    results = main()

    # Save results
    output = {
        "num_tests": len(results),
        "results": {
            name: {k: v for k, v in result.items() if k != "cascade_log10"}
            for name, result in results.items()
        }
    }

    with open("/Users/scones/projects/riemann/operators/track_crossblock/cascade_beurling_real_results.json", "w") as f:
        json.dump(output, f, indent=2)

    print(f"\nResults saved to cascade_beurling_real_results.json")
