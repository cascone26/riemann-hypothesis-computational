#!/usr/bin/env python3
"""
Test D_kl within its valid domain: t ∈ [0, 2a].
"""
import mpmath as mp
mp.mp.dps = 30

a = mp.mpf('0.549')
beta = mp.pi / (2*a)
twoa = 2*a

print(f"a = {mp.nstr(a, 8)}, 2a = {mp.nstr(twoa, 6)}")
print(f"Valid domain: t ∈ [0, {mp.nstr(twoa, 6)}]\n")

def vprime(k, x):
    """v_k'(x) with x ∈ [0, 2a]"""
    return k * beta * mp.cos(k * beta * x)

def C_kl(k, l, t_val):
    """C_kl(t) = ∫ v_k'(x) v_l'(x-t) dx"""
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
    """D_kl(t) = C_kl(t) + C_lk(t)"""
    return C_kl(k, l, t_val) + C_kl(l, k, t_val)

def D_kl_closedform(k, l, t_val):
    """Closed-form D_kl(t)"""
    if k == l:
        return k*k*beta*beta*(twoa - t_val)*mp.cos(k*beta*t_val) - k*beta*mp.sin(k*beta*t_val)
    else:
        return k*l*beta * (
            (mp.sin(l*beta*t_val) - mp.sin(k*beta*t_val))/(k-l)
            - (mp.sin(l*beta*t_val) + mp.sin(k*beta*t_val))/(k+l)
        )

print("Testing D_kl within valid range t ∈ [0, 2a]\n")
print(f"{'k,l':<7} {'t':<8} {'definition':<18} {'closed-form':<18} {'rel.error':<10}")
print("-" * 70)

test_pairs = [(1,1), (3,3), (5,5), (7,7), (1,3), (3,5), (5,7)]
# Use fractions of 2a instead of fixed values
test_t_fracs = [mp.mpf(i)/10 for i in range(1, 11)]  # 0.1, 0.2, ..., 1.0 of 2a
test_t = [f * twoa for f in test_t_fracs]

max_rel_err = mp.mpf(0)
max_abs_err = mp.mpf(0)
for k, l in test_pairs:
    for t_val in test_t:
        d_def = D_kl_definition(k, l, t_val)
        d_cf = D_kl_closedform(k, l, t_val)
        abs_err = abs(d_def - d_cf)
        if abs(d_def) > 1e-20:
            rel_err = abs_err / abs(d_def)
        else:
            rel_err = abs_err  # absolute error for near-zero values

        max_abs_err = max(max_abs_err, abs_err)
        max_rel_err = max(max_rel_err, rel_err)

        print(f"{str(k)+','+str(l):<7} {mp.nstr(t_val, 5):<8} {mp.nstr(d_def, 16):<18} "
              f"{mp.nstr(d_cf, 16):<18} {mp.nstr(rel_err, 3):<10}")

print(f"\nMax absolute error: {mp.nstr(max_abs_err, 3)}")
print(f"Max relative error: {mp.nstr(max_rel_err, 3)}")
print(f"\nResult:")
if max_abs_err < 1e-12 and max_rel_err < 1e-10:
    print("✓ D_kl CLOSED FORM FULLY VERIFIED")
    print("  - Matches definition to high precision (>1e-10 relative)")
    print("  - Reduction D_kl = C_kl + C_lk is correct")
else:
    print(f"✗ D_kl closed form has errors")
    print(f"  Max abs: {mp.nstr(max_abs_err, 3)}, Max rel: {mp.nstr(max_rel_err, 3)}")
