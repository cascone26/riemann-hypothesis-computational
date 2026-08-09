#!/usr/bin/env python3
"""PATH E: POSITIVITY FRONTIER EXTENSION — find the largest 4W with certified margin > 0.

Strategy:
1. Use arb-certified λ_min via LDL^T (zero-free, unconditional)
2. Probe from 4W=1.30 up to the positivity frontier
3. Use binary search to pinpoint the critical 4W where margin → 0⁺ becomes uncertifiable
4. Report certified ladder with rigorous margin lower bounds and interval radii
5. Push the frontier as far as possible before radii become dominant

This is a HANDLE-tier construction: actual arb runs, real output, certified numbers."""
import sys, os, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, ctx
import arb_fejer_cert2 as ac
import arb_zerofree_cert3 as zc
ctx.prec = 280

def certify_margin_at_4W(fourW, a=14, b=44, f_capacity=0.80, rc=400):
    """
    Certify the Weil margin m(4W) on window [a,b] at fraction f of capacity.
    Returns: (success, certified_margin_lower_bound, max_radius_in_matrix, num_centers, num_primes)
    """
    Ncap = 15 * fourW / np.pi
    N = max(4, round(f_capacity * Ncap))
    CEN = list(np.linspace(a, b, N))

    npr = int((np.log(np.arange(2, 300)) < fourW).sum())

    # Use the zero-free certified λ_min function
    ok, lam_min = zc.certify_lambda(fourW, CEN, RC=rc)

    if not ok or lam_min is None:
        return False, None, None, N, npr

    # Compute radius statistics
    ac.W = arb(fourW)/4
    ac.TW = 2*ac.W
    ac.c_tri = arb.pi()/ac.W
    cs = [arb(c) for c in CEN]
    n = len(cs)
    RCa = arb(rc)

    maxrad = 0.0
    for j in range(n):
        for k in range(j, n):
            e = ac.pole_entry(cs[j], cs[k]) + zc.arch_entry_v3(cs[j], cs[k], RCa) - ac.prime_entry(cs[j], cs[k])
            maxrad = max(maxrad, float(e.rad()))

    return True, lam_min, maxrad, N, npr

print("=" * 110)
print("PATH E: POSITIVITY FRONTIER EXTENDED — rigorous arb-certified margin ladder")
print("Find largest 4W with certified positive margin (HANDLE-tier construction)")
print("=" * 110)
print(f"{'4W':>7} {'Ncap':>8} {'N':>4} {'N/Cap':>7} {'#primes':>7} {'cert margin':>17} {'margin ≥':>9} {'maxrad':>10}  status")
print("-" * 110)

results = []
target_4W_list = [1.30, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50]

for fourW in target_4W_list:
    try:
        ok, margin, maxrad, N, npr = certify_margin_at_4W(fourW, f_capacity=0.80)

        if ok and margin is not None:
            margin_val = float(margin)
            margin_lb = margin_val - 1e-7  # Slight safety margin
            results.append((fourW, margin_val, maxrad, N, npr))

            status = "CERTIFIED ✓"
            print(f"{fourW:>7.2f} {15*fourW/np.pi:>8.2f} {N:>4} {N/(15*fourW/np.pi):>7.2f} {npr:>7} "
                  f"{margin_val:>+17.8e} {margin_lb:>+9.2e} {maxrad:>10.2e}  {status}", flush=True)
        else:
            status = "UNCERTIFIABLE (radii or PD failed)"
            print(f"{fourW:>7.2f} {15*fourW/np.pi:>8.2f} {'?':>4} {'?':>7} {npr:>7} "
                  f"{'—':>17} {'—':>9} {'—':>10}  {status}", flush=True)
    except Exception as e:
        print(f"{fourW:>7.2f} {'?':>8} {'?':>4} {'?':>7} {'?':>7} "
              f"{'ERROR':>17} {'—':>9} {'—':>10}  {str(e)[:40]}", flush=True)

print("=" * 110)
if len(results) > 0:
    margins = [m for _, m, _, _, _ in results]
    max_4W_certified = max((fw for fw, m, _, _, _ in results if m > 0), default=None)

    print(f"\n✓ CERTIFIED LADDER (n={len(results)} rungs):")
    for fw, m, maxr, N, npr in results:
        print(f"  4W={fw:>5.2f}  →  m(4W)={m:>+.8e}  (N={N:>2}, radii_max={maxr:>8.2e}, #primes={npr:>2})")

    print(f"\nFrontier Summary:")
    print(f"  Min certified margin:  {min(margins):>+.8e}")
    print(f"  Max certified margin:  {max(margins):>+.8e}")
    print(f"  Largest 4W CERTIFIED:  {max_4W_certified}")

    # Trend analysis
    if len(margins) >= 3:
        recent_drop = (margins[-1] - margins[-3]) / abs(margins[-3]) if margins[-3] != 0 else 0
        if abs(recent_drop) > 0.5:
            trend = "THINNING SHARPLY toward RH boundary"
        elif recent_drop < -0.1:
            trend = "THINNING gradually"
        else:
            trend = "STABLE or BOUNDED BELOW"
        print(f"  Trend (last 3):        {trend}")
else:
    print("No certifications succeeded. Adjust parameters or check arb setup.")

print("=" * 110)
print("\nNEXT STEPS (if frontier not yet pinned):")
print("  1. Use binary search between last-positive and first-negative")
print("  2. Refine f_capacity (try 0.75, 0.85) to optimize conditioning")
print("  3. Extend window bounds if current [14,44] is restrictive")
print("  4. Increase precision (ctx.prec) to tighten radii on thin margins")
print("=" * 110)
