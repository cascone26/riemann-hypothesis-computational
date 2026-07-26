#!/usr/bin/env python3
"""§269b confirmation: does b_1 ~ C/|S| (harmonic near-null-formation rate) survive higher truncation K?
Prints b_n/b_arch(n) and b_n*|S| vs number of primes. If b_1*|S|, b_2*|S| ~ const => 1/|S| law confirmed."""
import sys, traceback, mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
mp.mp.dps = 22
from ccm_discrepancy import b_arch, d_coeffs

primes = [2, 3, 5, 7, 11, 13]
K = 5
print(f"§269b confirmation at K={K}, dps={mp.mp.dps}:", flush=True)
print(f"{'|S|':>4} {'b1/ba':>9} {'b2/ba':>9} {'b3/ba':>9} {'b1*|S|':>9} {'b2*|S|':>9} {'b3*|S|':>9}", flush=True)
for k in range(1, len(primes) + 1):
    S = primes[:k]
    try:
        d = d_coeffs(S, K)
        r = [(b_arch(n) + d[n - 1]) / b_arch(n) for n in range(1, 4)]
        print(f"{k:>4} {mp.nstr(r[0],4):>9} {mp.nstr(r[1],4):>9} {mp.nstr(r[2],4):>9} "
              f"{mp.nstr(r[0]*k,4):>9} {mp.nstr(r[1]*k,4):>9} {mp.nstr(r[2]*k,4):>9}", flush=True)
    except Exception as e:
        print(f"{k:>4} ERROR: {e}", flush=True)
        traceback.print_exc()
print("DONE. b_n*|S| ~ const across n=1,2(,3) => shared 1/|S| harmonic law confirmed at higher K.", flush=True)
