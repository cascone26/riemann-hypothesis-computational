#!/usr/bin/env python3
"""
MINIMAL TEST: Verify D_kl closed form only (avoid expensive Q integrals).
"""
import mpmath as mp
mp.mp.dps = 30

a = mp.mpf('0.549')
beta = mp.pi / (2*a)
twoa = 2*a

print(f"a = {mp.nstr(a, 8)}, beta = {mp.nstr(beta, 8)}\n")

def vprime(k, x):
    return k * beta * mp.cos(k * beta * x)

def D_kl_definition(k, l, t_val):
    """D_kl(t) by direct 1D integration."""
    if t_val < 0:
        t_val = -t_val
    def integrand(xp):
        y_shifted = xp - t_val
        return vprime(k, xp) * vprime(l, y_shifted)
    lower = max(mp.mpf(0), t_val)
    upper = twoa
    if lower >= upper:
        return mp.mpf(0)
    result = mp.quad(integrand, [lower, upper])
    return result

def D_kl_closedform(k, l, t_val):
    """Closed-form D_kl(t)."""
    if k == l:
        return k*k*beta*beta*(twoa - t_val)*mp.cos(k*beta*t_val) - k*beta*mp.sin(k*beta*t_val)
    else:
        return k*l*beta * (
            (mp.sin(l*beta*t_val) - mp.sin(k*beta*t_val))/(k-l)
            - (mp.sin(l*beta*t_val) + mp.sin(k*beta*t_val))/(k+l)
        )

print("Testing D_kl: definition (1D integral) vs closed form\n")
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
if max_err < 1e-10:
    print("✓ D_kl closed form VERIFIED - matches 1D integral from definition to >1e-10")
else:
    print(f"✗ D_kl closed form has error: {mp.nstr(max_err, 3)}")
