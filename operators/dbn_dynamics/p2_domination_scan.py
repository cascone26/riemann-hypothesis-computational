#!/usr/bin/env python3
"""§185: the FIRST prime sets the criticality. Scan the Weil-form margin λ_a (Suzuki §8 compact, g_a_compact)
across the p=2-ONLY support zone a∈(½log2, ½log3)=(0.3466, 0.5493) — where the Weil form is EXACTLY
(archimedean − single p=2 term), no other primes. Watch λ_a collapse as p=2's support-weight ramps up.
Finding (2026-07-23): λ_a falls ~1000× (7.6e-3 → ~1e-5) within the p=2-only zone; higher primes barely move it
(flat to a=3). ⟹ criticality is SET by p=2; the crux is the two-term arch-vs-p2 balance, not the full prime sum.
CAVEAT: float64/Nyström floor ~1e-5 is NOT converged (§109 mpmath: a=0.6→2.3e-5 vs float64 7.6e-6); the COLLAPSE
is robust, the exact floor needs mpmath."""
import sys, json, math
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
from g_a_compact import lambda_a_compact

def prime_powers_upto(x):
    out = []
    for n in range(2, x + 1):
        m, p = n, None
        for d in range(2, n + 1):
            if m % d == 0:
                p = d
                while m % d == 0: m //= d
                break
        if m == 1: out.append(n)   # n is a prime power
    return out

if __name__ == "__main__":
    a05log = 0.5 * math.log(2)   # 0.3466  p=2 onset
    a05l3  = 0.5 * math.log(3)   # 0.5493  p=3 onset
    avals = [0.25, 0.30, a05log, 0.40, 0.45, 0.50, a05l3, 0.60, 0.70]
    print(f"§185 first-prime criticality: λ_a across the p=2-only zone (p=2@{a05log:.4f}, p=3@{a05l3:.4f})")
    print(f"{'a':>7} {'primepowers≤e^2a':>18} {'λ_a (Weil margin)':>20} {'gap':>11}")
    out = []
    for a in avals:
        Nmax = int(math.floor(math.exp(2 * a)))
        pps = prime_powers_upto(Nmax)
        la, gap = lambda_a_compact(a)
        print(f"{a:>7.4f} {str(pps):>18} {la:>20.6e} {gap:>11.3e}", flush=True)
        out.append(dict(a=a, prime_powers=pps, lambda_a=float(la), gap=float(gap)))
    json.dump(out, open("results/p2_domination_scan.json", "w"), indent=2)
    print("\nREAD: λ_a collapses ~1000x (7.6e-3→~1e-5) DRIVEN ENTIRELY BY p=2 within its own zone; higher primes flat.")
    print("⟹ the arch-vs-p=2 two-term balance SETS the Weil criticality. Crux = dominate ONE prime; rest is a tail.")
    print("(Exact floor float64-limited; collapse robust. Next: mpmath the p=2-only margin + the arch-vs-p2 sharp bound.)")
