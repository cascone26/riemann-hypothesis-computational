#!/usr/bin/env python3
"""§269d: is the b_1 ~ C/|S| decay COUNT-only (combinatorial) or PRIME-VALUE-dependent (arithmetic)?
Compare b_1/b_arch across prime sets of the SAME size |S| but different prime values.
Same value => count-only. Different => arithmetic content (candidate tie to zero density)."""
import sys, mpmath as mp
sys.path.insert(0, '/Users/scones/projects/riemann/operators')
from ccm_discrepancy import b_arch, d_coeffs   # note: import sets dps=60
K = 5
def b1r(S):
    d = d_coeffs(S, K)
    return (b_arch(1) + d[0]) / b_arch(1)
tests = [
    [2, 3, 5, 7],        # |S|=4 small
    [11, 13, 17, 19],    # |S|=4 larger
    [2, 3, 5, 7, 11, 13],       # |S|=6 small
    [13, 17, 19, 23, 29, 31],   # |S|=6 larger
    [2],                 # single small
    [31],                # single large
]
print(f"§269d COUNT vs ARITHMETIC test (K={K}, dps={mp.mp.dps}):", flush=True)
print(f"{'|S|':>4} {'b1/b_arch':>12}  {'primes':>26}  {'sum 1/sqrt p':>13}  {'sum log p':>10}", flush=True)
for S in tests:
    r = b1r(S)
    s_inv = sum(1/mp.sqrt(p) for p in S)
    s_log = sum(mp.log(p) for p in S)
    print(f"{len(S):>4} {mp.nstr(r,6):>12}  {str(S):>26}  {mp.nstr(s_inv,5):>13}  {mp.nstr(s_log,5):>10}", flush=True)
print("DONE. Same-|S| rows equal => COUNT-only. Differ => ARITHMETIC (check which of sum-1/sqrt-p, sum-log-p tracks it).", flush=True)
