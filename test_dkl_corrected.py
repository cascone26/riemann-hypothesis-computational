#!/usr/bin/env python3
"""
CORRECTED: Test D_kl = C_kl + C_lk (both directions).
The closed form computes D_kl(t) which includes BOTH:
  C_kl(t) = ∫ v_k'(x) v_l'(x-t) dx
  C_lk(t) = ∫ v_l'(x) v_k'(x-t) dx
My first test only computed C_kl, not the sum. Fixed now.
"""
import mpmath as mp
mp.mp.dps = 30

a = mp.mpf('0.549')
beta = mp.pi / (2*a)
twoa = 2*a

print(f"a = {mp.nstr(a, 8)}, beta = {mp.nstr(beta, 8)}\n")

def vprime(k, x):
    """v_k'(x) with x ∈ [0, 2a]"""
    return k * beta * mp.cos(k * beta * x)

def C_kl(k, l, t_val):
    """C_kl(t) = ∫_0^{2a} v_k'(x) v_l'(x-t) dx, using x ∈ [max(0,t), 2a]"""
    if t_val < 0:
        t_val = -t_val
    def integrand(xp):
        y_shifted = xp - t_val
        return vprime(k, xp) * vprime(l, y_shifted)
    lower = max(mp.mpf(0), t_val)
    upper = twoa
    if lower >= upper:
        return mp.mpf(0)
    return mp.quad(integrand, [lower, upper])

def D_kl_definition(k, l, t_val):
    """D_kl(t) = C_kl(t) + C_lk(t) from the definition"""
    return C_kl(k, l, t_val) + C_kl(l, k, t_val)

def D_kl_closedform(k, l, t_val):
    """Closed-form D_kl(t) from dkl_closedform.py"""
    if k == l:
        return k*k*beta*beta*(twoa - t_val)*mp.cos(k*beta*t_val) - k*beta*mp.sin(k*beta*t_val)
    else:
        return k*l*beta * (
            (mp.sin(l*beta*t_val) - mp.sin(k*beta*t_val))/(k-l)
            - (mp.sin(l*beta*t_val) + mp.sin(k*beta*t_val))/(k+l)
        )

print("Testing D_kl = C_kl + C_lk: definition vs closed form\n")
print(f"{'k,l':<7} {'t':<8} {'definition':<18} {'closed-form':<18} {'error':<10}")
print("-" * 70)

test_pairs = [(1,1), (3,3), (5,5), (7,7), (1,3), (3,5), (5,7)]
test_t = [mp.mpf('0.05'), mp.mpf('0.25'), mp.mpf('0.5'), mp.mpf('1.0'), mp.mpf('1.5')]

max_err = mp.mpf(0)
for k, l in test_pairs:
    for t_val in test_t:
        d_def = D_kl_definition(k, l, t_val)
        d_cf = D_kl_closedform(k, l, t_val)
        err = abs(d_def - d_cf)
        max_err = max(max_err, err)
        print(f"{str(k)+','+str(l):<7} {mp.nstr(t_val,4):<8} {mp.nstr(d_def, 16):<18} "
              f"{mp.nstr(d_cf, 16):<18} {mp.nstr(err, 3):<10}")

print(f"\nMax error: {mp.nstr(max_err, 3)}")
print(f"\nResult:")
if max_err < 1e-8:
    print("✓ D_kl closed form VERIFIED - matches definition to >1e-8")
else:
    print(f"✗ D_kl closed form has error: {mp.nstr(max_err, 3)}")
