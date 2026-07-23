#!/usr/bin/env python3
"""§189: is Weil positivity a NORM-BUDGET RELAY? §188: p=2 place = involution of exact norm w=Λ(2)/√2; arch
negativity |α₁(a)| deepens with support a; p=3 enters at a=½log3 as the p=2 budget nears exhaustion.
TEST across the prime thresholds (½log2=.347, ½log3=.549, ½log4=.693, ½log5=.805): does the arch-negativity
depth |α₁(a)| stay COVERED by the cumulative prime budget B(a)=Σ_{n≤e^2a}Λ(n)/√n, with each prime entering as
|α₁| approaches the running budget? Report arch λ_min, full λ_min, B(a), and |arch λ_min|/B(a). If |α₁|≈B (ratio
→1) and full stays ~+2.5e-5, the relay is the ANATOMY of A+B≥P: primes relay to track arch-negativity growth."""
import sys, math
import mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym

def prime_budget(a):
    """B(a) = Σ_{prime powers n ≤ e^{2a}} Λ(n)/√n."""
    Nmax = int(math.floor(math.exp(2 * float(a))))
    tot = mp.mpf(0); contribs = []
    for n in range(2, Nmax + 1):
        m, p = n, None
        for d in range(2, n + 1):
            if m % d == 0:
                p = d
                while m % d == 0: m //= d
                break
        if m == 1:                      # n = p^k prime power, Λ(n)=log p
            wn = mp.log(p) / mp.sqrt(n); tot += wn; contribs.append((n, float(wn)))
    return tot, contribs

if __name__ == "__main__":
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    K   = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    M   = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    avals = [mp.mpf(x) for x in (sys.argv[4].split(',') if len(sys.argv) > 4 else
             ['0.45','0.549','0.60','0.693','0.75','0.805','0.85'])]
    mp.mp.dps = dps
    print(f"§189 norm-budget relay test | dps={dps} K={K} M={M}")
    print(f"{'a':>6} {'primepow≤e2a':>16} {'arch λ_min':>13} {'full λ_min':>13} {'budget B(a)':>13} {'|archλmin|/B':>13}")
    for a in avals:
        Qf, Qa, Mass, Vv, xs, ww = weil_matrices(a, K, M)
        Ea, _ = gen_eig_sym(Qa, Mass)
        Ef, _ = gen_eig_sym(Qf, Mass)
        B, contribs = prime_budget(a)
        npp = [n for n, _ in contribs]
        ratio = abs(Ea[0]) / B if B > 0 else mp.mpf('nan')
        print(f"{mp.nstr(a,4):>6} {str(npp):>16} {mp.nstr(Ea[0],6):>13} {mp.nstr(Ef[0],5):>13} {mp.nstr(B,6):>13} {mp.nstr(ratio,5):>13}", flush=True)
    print("\nREAD: if |arch λ_min|/B(a) → 1 (budget just covers negativity) with full λ_min ~+2.5e-5 flat across every")
    print("prime threshold ⟹ the RELAY is exact: Σ prime norms tracks arch-negativity growth = anatomy of A+B≥P.")
    print("If the ratio wanders / full λ_min drifts, the relay is only heuristic — record which.")
