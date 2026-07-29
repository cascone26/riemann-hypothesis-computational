"""Tao: does the de Bruijn-Newman reality-margin die as t or sqrt(t)?
Destabilizing heat operator q_s = e^{+s d^2/dx^2} p = sum_k s^k/k! p^{(2k)}.
An isolated real pair at gap delta goes complex at s* = delta^2/8 (analytic).
Test with FULL multi-zero configs: measure s*(delta) for a near-pair embedded in
a realistic spacing, confirm s* ~ delta^2 (beta=2) and read the constant."""
import numpy as np
np.set_printoptions(suppress=True)

def heat_forward_destabilize(coeffs, s, terms=40):
    # coeffs: polynomial highest-first (np.poly1d convention). q_s = sum_k s^k/k! p^{(2k)}
    p = np.poly1d(coeffs); q = np.poly1d(p.coeffs.astype(float)); acc = np.poly1d(q.coeffs.copy())
    d2 = p; fact = 1.0
    for k in range(1, terms):
        d2 = d2.deriv().deriv()
        if len(d2.coeffs)==0: break
        fact *= k
        acc = acc + (s**k/fact)*d2
    return acc

def first_complex_s(roots, s_hi=5.0, tol=1e-9):
    base = np.poly1d(roots, r=True).coeffs
    lo, hi = 0.0, s_hi
    def all_real(s):
        r = heat_forward_destabilize(base, s).roots
        return np.max(np.abs(r.imag)) < 1e-6
    if all_real(hi):  # never goes complex in range
        return None
    while hi-lo > tol:
        mid=(lo+hi)/2
        if all_real(mid): lo=mid
        else: hi=mid
    return (lo+hi)/2

print("ISOLATED pair (analytic s*=delta^2/8):")
for delta in [0.4,0.2,0.1,0.05]:
    s=first_complex_s([-delta/2, delta/2])
    print(f"  delta={delta:.3f}  s*={s:.6e}  s*/delta^2={s/delta**2:.4f}  (analytic 0.125)")

print("\nEMBEDDED pair in a spacing-~1 comb of 8 zeros (neighbors present):")
for delta in [0.30,0.20,0.10,0.05]:
    roots=[-3.0,-2.0,-1.0,-delta/2,delta/2,1.0,2.0,3.0]
    s=first_complex_s(roots)
    print(f"  delta={delta:.3f}  s*={s:.6e}  s*/delta^2={s/delta**2:.4f}")

# calibrate the UPPER-BOUND wall X ~ exp(c/sqrt(t)) from the two record points
import math
def c_from(Lam,X): return math.log(X)*math.sqrt(Lam)
c1=c_from(0.22, 1e10); c2=c_from(0.20, 3e12)
print(f"\nUpper-bound wall X~exp(c/sqrt(Lambda)):  c(0.22@1e10)={c1:.1f}  c(0.20@3e12)={c2:.1f}")
c=c2
for L in [0.10,0.05,0.01]:
    print(f"  Lambda<= {L:.2f}  needs  log X={c/math.sqrt(L):.1f}  X~1e{c/math.sqrt(L)/math.log(10):.1f}")
