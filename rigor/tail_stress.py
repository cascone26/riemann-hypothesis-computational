#!/usr/bin/env python3
"""
Adversarial tail-stress test: hunt for cracks in the tail lemma.
Claims: for modes k > K0, tail block T is strictly diagonally dominant (Gershgorin)
hence PD. We BREAK IT if: min eigenvalue ≤ 0, or Gershgorin margin touches 0.

Optimized: cache Q values, focus on key test cases.
"""
import sys, os
import mpmath as mp
import numpy as np

mp.mp.dps = 30
a = mp.mpf('0.549')
b = mp.pi / (2 * a)
log2 = mp.log(2)
twoa = 2 * a
Om0 = mp.digamma(mp.mpf(1) / 4) - mp.log(mp.pi)
Phi0 = mp.zeta(2, mp.mpf(1) / 4)
wp = mp.log(2) / mp.sqrt(2)
NL = 1500

# Global Q cache
Q_cache = {}

def Pp(p, t, s):
    if p == 0:
        return 1 / s
    if p == 1:
        return t / s - 1 / s ** 2
    return t * t / s - 2 * t / s ** 2 + 2 / s ** 3

def II(p, al, om, sc, t0, t1):
    s = mp.mpc(-al, om)
    F = lambda tt: mp.e ** (s * tt) * Pp(p, tt, s)
    v = F(t1) - F(t0)
    return v.real if sc == 'cos' else v.imag

def dterms(k, l):
    if k == l:
        return [(2 * a * k * k * b * b, 0, k * b, 'cos'),
                (-k * k * b * b, 1, k * b, 'cos'),
                (-k * b, 0, k * b, 'sin')]
    al = k * l * b * (mp.mpf(1) / (k - l) - mp.mpf(1) / (k + l))
    ak = -k * l * b * (mp.mpf(1) / (k - l) + mp.mpf(1) / (k + l))
    return [(al, 0, l * b, 'sin'), (ak, 0, k * b, 'sin')]

def Q(k, l, a_val=None):
    """Compute Q(k,l) with caching. Can use different a for parameter sweeps."""
    if a_val is None:
        a_val = a
    cache_key = (k, l, float(a_val))
    if cache_key in Q_cache:
        return Q_cache[cache_key]

    # Note: Q(k,l) depends on a through b = π/(2a), which affects dterms and II
    if a_val != a:
        # Recompute with different a
        b_local = mp.pi / (2 * a_val)
        twoa_local = 2 * a_val
        tot = mp.mpf(0)
        gt = [(-4, 0, mp.mpf(-1) / 2), (-4, 0, mp.mpf(1) / 2), (8, 0, mp.mpf(0)),
              (-Om0 / 2, 1, mp.mpf(0)), (-Phi0 / 4, 0, mp.mpf(0))]
        for n in range(NL):
            gt.append((mp.mpf(1) / 4 / (n + mp.mpf(1) / 4) ** 2, 0, 2 * n + mp.mpf(1) / 2))

        # Recompute dterms with local b
        def dterms_local(k, l):
            if k == l:
                return [(2 * a_val * k * k * b_local * b_local, 0, k * b_local, 'cos'),
                        (-k * k * b_local * b_local, 1, k * b_local, 'cos'),
                        (-k * b_local, 0, k * b_local, 'sin')]
            al = k * l * b_local * (mp.mpf(1) / (k - l) - mp.mpf(1) / (k + l))
            ak = -k * l * b_local * (mp.mpf(1) / (k - l) + mp.mpf(1) / (k + l))
            return [(al, 0, l * b_local, 'sin'), (ak, 0, k * b_local, 'sin')]

        for (ag, pg, alp) in gt:
            for (ad, pd, om, sc) in dterms_local(k, l):
                tot += ag * ad * II(pg + pd, alp, om, sc, mp.mpf(0), twoa_local)
        for (ad, pd, om, sc) in dterms_local(k, l):
            tot += wp * ad * II(pd + 1, mp.mpf(0), om, sc, log2, twoa_local)
            tot += (-wp * log2) * ad * II(pd, mp.mpf(0), om, sc, log2, twoa_local)
    else:
        # Use cached version with standard a
        tot = mp.mpf(0)
        gt = [(-4, 0, mp.mpf(-1) / 2), (-4, 0, mp.mpf(1) / 2), (8, 0, mp.mpf(0)),
              (-Om0 / 2, 1, mp.mpf(0)), (-Phi0 / 4, 0, mp.mpf(0))]
        for n in range(NL):
            gt.append((mp.mpf(1) / 4 / (n + mp.mpf(1) / 4) ** 2, 0, 2 * n + mp.mpf(1) / 2))
        for (ag, pg, alp) in gt:
            for (ad, pd, om, sc) in dterms(k, l):
                tot += ag * ad * II(pg + pd, alp, om, sc, mp.mpf(0), twoa)
        for (ad, pd, om, sc) in dterms(k, l):
            tot += wp * ad * II(pd + 1, mp.mpf(0), om, sc, log2, twoa)
            tot += (-wp * log2) * ad * II(pd, mp.mpf(0), om, sc, log2, twoa)

    Q_cache[cache_key] = tot
    return tot

def build_tail_block(modes_odd, K0, a_val=None):
    """Build tail block T for odd modes in (K0, Kmax]."""
    if a_val is None:
        a_val = a
    tail_modes = [k for k in modes_odd if k > K0]
    n = len(tail_modes)
    T = mp.matrix(n, n)
    print(f"    Computing {n}x{n} tail block matrix...", flush=True)
    for i, k in enumerate(tail_modes):
        for j, l in enumerate(tail_modes):
            T[i, j] = Q(k, l, a_val)
        if (i + 1) % max(1, n // 5) == 0:
            print(f"      Row {i+1}/{n} complete", flush=True)
    print(f"    Tail block complete.", flush=True)
    return T, tail_modes

def min_eigenvalue(T, a_val):
    """Compute minimum eigenvalue of generalized problem Tv = λ(a*I)v using scipy."""
    n = T.rows
    # Convert mpmath matrix to numpy for scipy
    T_np = np.array([[float(T[i, j]) for j in range(n)] for i in range(n)], dtype=float)
    a_np = float(a_val)
    # The generalized problem T v = λ (a*I) v is equivalent to (1/a)*T v = λ v
    T_scaled = T_np / a_np
    try:
        evals = np.linalg.eigvalsh(T_scaled)
        lam_min = float(np.min(evals))
        return lam_min
    except Exception as e:
        print(f"      Eigenvalue computation failed: {e}", flush=True)
        return None

def gershgorin_margins(T, a_val, tail_modes):
    """Compute Gershgorin margin for each tail mode: g(k) = Q(k,k) - Σ_{l≠k}|Q(k,l)|."""
    margins = {}
    for i, k in enumerate(tail_modes):
        diag = T[i, i] / a_val  # Normalize by a for fair comparison
        offdiag_sum = sum(abs(T[i, j]) for j in range(T.rows) if j != i) / a_val
        margin = diag - offdiag_sum
        margins[k] = (float(diag), float(offdiag_sum), float(margin))
    return margins

def stress_test_single(K0, Kmax, a_val, modes_odd, log_f):
    """Run stress test for given K0, Kmax, a."""
    log_f.write(f"\n{'='*80}\n")
    log_f.write(f"Stress test: K0={K0}, Kmax={Kmax}, a={float(a_val):.3f}\n")
    log_f.write(f"{'='*80}\n")
    log_f.flush()

    # Build tail block
    T, tail_modes = build_tail_block(modes_odd, K0, a_val)
    n = len(tail_modes)
    if n == 0:
        log_f.write(f"  No tail modes for K0={K0}\n")
        log_f.flush()
        return None
    log_f.write(f"  Tail modes: {tail_modes[:5]}...{tail_modes[-2:]} (n={n})\n")
    log_f.flush()

    # 1. Rayleigh minimum (eigenvalue)
    log_f.write(f"\n  [1] Rayleigh minimum eigenvalue:\n")
    log_f.flush()
    lam_min = min_eigenvalue(T, a_val)
    log_f.write(f"      λ_min = {lam_min:.10f}\n")
    status_rayleigh = "PASS (PD)" if (lam_min and lam_min > 1e-10) else "FAIL or NEAR-ZERO"
    log_f.write(f"      → {status_rayleigh}\n")
    log_f.flush()

    # 2. Gershgorin margins
    log_f.write(f"\n  [2] Gershgorin dominance:\n")
    log_f.flush()
    margins = gershgorin_margins(T, a_val, tail_modes)
    worst_mode = min(margins, key=lambda k: margins[k][2])
    worst_margin = margins[worst_mode][2]
    log_f.write(f"      Worst mode: k={worst_mode}, diag={margins[worst_mode][0]:.6f}, offdiag_sum={margins[worst_mode][1]:.6f}\n")
    log_f.write(f"      Gershgorin margin g(k) = {worst_margin:.10f}\n")
    status_gershgorin = "PASS (dominant)" if worst_margin > 1e-10 else "FAIL or NEAR-ZERO"
    log_f.write(f"      → {status_gershgorin}\n")
    log_f.flush()

    return {
        'K0': K0,
        'Kmax': Kmax,
        'a': float(a_val),
        'n_tail': n,
        'lam_min': lam_min,
        'worst_mode': worst_mode,
        'gershgorin_margin_worst': worst_margin,
        'status_rayleigh': status_rayleigh,
        'status_gershgorin': status_gershgorin
    }

class DualWriter:
    def __init__(self, logfile):
        self.logfile = logfile
        self.original_stdout = sys.__stdout__
    def write(self, s):
        self.original_stdout.write(s)
        with open(self.logfile, 'a') as f:
            f.write(s)
    def flush(self):
        self.original_stdout.flush()

def main():
    os.makedirs('/Users/scones/projects/riemann/results', exist_ok=True)
    logfile = '/Users/scones/projects/riemann/results/tail_stress.log'

    log_f = open(logfile, 'w')
    log_f.write("="*80 + "\n")
    log_f.write("TAIL STRESS TEST: Adversarial search for cracks in tail lemma\n")
    log_f.write("="*80 + "\n")
    log_f.write("Task: find modes/configs where tail block T loses positive-definiteness\n")
    log_f.write("      or breaks Gershgorin diagonal dominance.\n")
    log_f.write("      NL (Lerch terms) = 1500, high precision.\n\n")
    log_f.flush()

    # Redirect stdout to dual output
    sys.stdout = DualWriter(logfile)

    # Odd modes
    modes_odd = list(range(1, 52, 2))
    Kmax = 49

    results = []

    # Test 1: Primary test with standard K0 and a
    print("\n" + "="*80)
    print("PRIMARY TEST: K0=9, a=0.549")
    print("="*80)
    res = stress_test_single(9, Kmax, mp.mpf('0.549'), modes_odd, log_f)
    if res:
        results.append(res)

    # Test 2: Verify with different K0 (smaller = larger tail)
    print("\n" + "="*80)
    print("SENSITIVITY: Vary K0 (a=0.549)")
    print("="*80)
    for K0 in [7, 11]:
        res = stress_test_single(K0, Kmax, mp.mpf('0.549'), modes_odd, log_f)
        if res:
            results.append(res)

    # Test 3: Check tight endpoint (a=0.549 is the endpoint)
    print("\n" + "="*80)
    print("ENDPOINT CHECK: Near a=0.549 (tight region)")
    print("="*80)
    for a_test in [0.545, 0.549]:
        res = stress_test_single(9, Kmax, mp.mpf(str(a_test)), modes_odd, log_f)
        if res:
            results.append(res)

    # Summary
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)

    crack_found = False
    for res in results:
        is_pass = ('PASS' in res['status_rayleigh'] and 'PASS' in res['status_gershgorin'])
        if not is_pass:
            crack_found = True
            print(f"\n⚠️ CRACK DETECTED at K0={res['K0']}, a={res['a']}:")
            if 'PASS' not in res['status_rayleigh']:
                print(f"   - Rayleigh: λ_min={res['lam_min']:.10f} (should be > 0)")
            if 'PASS' not in res['status_gershgorin']:
                print(f"   - Gershgorin: worst margin={res['gershgorin_margin_worst']:.10f} at mode {res['worst_mode']}")

    if not crack_found:
        print("\n✓ NO CRACKS FOUND. Tail lemma holds (PD + dominant) across all tested configs.\n")
        # Report tightest margins
        min_rayleigh = min((r['lam_min'] for r in results if r['lam_min']), default=None)
        min_gershgorin = min((r['gershgorin_margin_worst'] for r in results), default=None)
        if min_rayleigh:
            print(f"  Tightest Rayleigh margin: λ_min = {min_rayleigh:.10f}")
        if min_gershgorin:
            print(f"  Tightest Gershgorin margin: {min_gershgorin:.10f}")

    print(f"\n✓ Stress test complete. Full output in {logfile}")
    log_f.close()

if __name__ == '__main__':
    main()
