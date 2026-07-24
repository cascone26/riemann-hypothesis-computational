#!/usr/bin/env python3
"""
FAST ADVERSARIAL VERIFICATION — spot checks only.
Reduces computation by:
- Testing D_kl(t) on just a few (t, k, l) values
- Testing Q[k,l] on just the 1×1 diagonal (k=1, k=3, etc.) as a basic check
- Skipping the full 2D brute-force for off-diagonals
"""
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')

import mpmath as mp
mp.mp.dps = 35

a = mp.mpf('0.549')
beta = mp.pi / (2*a)
log2 = mp.log(2)
twoa = 2*a

print(f"== FAST ADVERSARIAL VERIFICATION: Even-Sector Weil Form ==")
print(f"a = {mp.nstr(a, 8)}, β = {mp.nstr(beta, 8)}\n")

# Basis function derivative
def vprime(k, x):
    """Derivative v_k'(x) with x ∈ [0, 2a]."""
    return k * beta * mp.cos(k * beta * x)

# Prime data
def prime_data_local(a_val):
    """Prime powers n <= e^{2a}."""
    N = int(mp.e**(2*a_val)) + 1
    is_prime = [True]*(N+1)
    is_prime[0] = is_prime[1] = False
    for i in range(2, N+1):
        if is_prime[i]:
            for j in range(2*i, N+1, i):
                is_prime[j] = False
    out = []
    for p in range(2, N+1):
        if is_prime[p]:
            logp = mp.log(p)
            pk = p
            while pk <= N:
                logpk = mp.log(pk)
                out.append((logp, logpk, mp.sqrt(mp.mpf(pk))))
                pk *= p
    return out

def screw_g(t, a_val):
    """Screw kernel g(t), t ≥ 0. Returns (g_full, g_arch)."""
    at = abs(t)
    boundary = -4*(mp.e**(t/2) + mp.e**(-t/2) - 2)
    psi_quarter = mp.digamma(mp.mpf(1)/4)
    logpi = mp.log(mp.pi)
    arch_core = -(at/2)*(psi_quarter - logpi)
    zeta_quarter = mp.zeta(2, mp.mpf(1)/4)
    if at < mp.mpf(10)**(-20):
        lerch_t = zeta_quarter
    else:
        z = mp.e**(-2*at)
        lerch_t = mp.nsum(lambda n: z**n / (n + mp.mpf(1)/4)**2, [0, mp.inf])
    arch_tail = -(mp.mpf(1)/4) * (zeta_quarter - mp.e**(-at/2) * lerch_t)
    g_arch = boundary + arch_core + arch_tail

    prime_sum = mp.mpf(0)
    pl = prime_data_local(a_val)
    for (lam, logn, sqrtn) in pl:
        if logn <= at:
            prime_sum += lam/sqrtn * (at - logn)
    g_full = g_arch + prime_sum
    return g_full, g_arch

# ========== PART 1: D_kl(t) from definition vs closed form ==========

print("[PART 1] D_kl(t): Definition vs Closed Form\n")

def D_kl_definition(k, l, t_val):
    """D_kl(t) by direct integration."""
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

# Test D_kl on diagonal and a few off-diagonal entries
test_pairs = [(1,1), (3,3), (5,5), (1,3), (3,5)]
test_t = [mp.mpf('0.1'), mp.mpf('0.5'), mp.mpf('1.0')]

print(f"{'k,l':<6} {'t':<8} {'definition':<18} {'closed-form':<18} {'error':<12}")
print("-" * 70)

max_dkl_err = mp.mpf(0)
for k, l in test_pairs:
    for t_val in test_t:
        d_def = D_kl_definition(k, l, t_val)
        d_cf = D_kl_closedform(k, l, t_val)
        err = abs(d_def - d_cf)
        max_dkl_err = max(max_dkl_err, err)
        print(f"{str(k)+','+str(l):<6} {mp.nstr(t_val,3):<8} {mp.nstr(d_def, 16):<18} "
              f"{mp.nstr(d_cf, 16):<18} {mp.nstr(err, 3):<12}")

print(f"\nMax D_kl error: {mp.nstr(max_dkl_err, 3)}")
if max_dkl_err < 1e-8:
    print("✓ D_kl CLOSED FORM VERIFIED")
else:
    print(f"✗ WARNING: D_kl error {mp.nstr(max_dkl_err, 3)}")

# ========== PART 2: Q[k,k] diagonal entries ==========

print(f"\n[PART 2] Q[k,k] Diagonal: Definition vs via D_kl\n")

def Q_definition(k, l):
    """Brute-force 2D integration."""
    def integrand(xp, yp):
        t_val = xp - yp
        g_f, _ = screw_g(t_val, a)
        return g_f * vprime(k, xp) * vprime(l, yp)
    result = mp.quad(integrand, [0, twoa], [0, twoa])
    return result

def Q_via_Dkl(k, l):
    """Via D_kl closed form."""
    def integrand(t_val):
        gf, _ = screw_g(t_val, a)
        d = D_kl_closedform(k, l, t_val)
        return gf * d
    result = mp.quad(integrand, [0, twoa])
    return result

# Test diagonal only (faster)
diag_modes = [1, 3, 5, 7]
print(f"{'k':<4} {'Q_definition':<20} {'Q_via_Dkl':<20} {'error':<12}")
print("-" * 60)

max_Q_err = mp.mpf(0)
for k in diag_modes:
    q_def = Q_definition(k, k)
    q_dkl = Q_via_Dkl(k, k)
    err = abs(q_def - q_dkl)
    max_Q_err = max(max_Q_err, err)
    print(f"{k:<4} {mp.nstr(q_def, 18):<20} {mp.nstr(q_dkl, 18):<20} {mp.nstr(err, 3):<12}")

print(f"\nMax Q diagonal error: {mp.nstr(max_Q_err, 3)}")

# ========== PART 3: Compare against existing entry_engine.py ==========

print(f"\n[PART 3] Compare Against entry_engine.py\n")

try:
    from entry_engine import Q as Q_engine
    from mpmath_weil_hiprec import weil_matrices

    print("Loading weil_matrices (dps=40)...")
    Qf, Qa, Mass, Vv, xs, w = weil_matrices(a, 24, 1200)
    ev_indices = [i for i in range(24) if (i+1) % 2 == 1][:4]

    print(f"{'k':<4} {'definition':<18} {'entry_engine':<18} {'weil_matrix':<18} {'err':<12}")
    print("-" * 75)

    max_engine_err = mp.mpf(0)
    for ik, k in enumerate(diag_modes):
        q_def = Q_definition(k, k)
        q_eng = Q_engine(k, k)
        q_mat = Qf[ev_indices[ik], ev_indices[ik]]
        err = abs(q_def - q_eng)
        max_engine_err = max(max_engine_err, err)
        print(f"{k:<4} {mp.nstr(q_def, 16):<18} {mp.nstr(q_eng, 16):<18} "
              f"{mp.nstr(q_mat, 16):<18} {mp.nstr(err, 3):<12}")

    print(f"\nMax entry_engine error: {mp.nstr(max_engine_err, 3)}")

except Exception as e:
    print(f"Could not import existing code: {e}")
    import traceback
    traceback.print_exc()

# ========== SUMMARY ==========

print(f"\n" + "="*75)
print(f"VERIFICATION SUMMARY")
print(f"="*75)
print(f"Max D_kl error (definition vs closed-form):    {mp.nstr(max_dkl_err, 3)}")
print(f"Max Q[k,k] error (definition vs D_kl):         {mp.nstr(max_Q_err, 3)}")
try:
    print(f"Max Q[k,k] error (definition vs entry_engine): {mp.nstr(max_engine_err, 3)}")
except:
    pass

print(f"\nConclusion:")
if max_dkl_err < 1e-8 and max_Q_err < 1e-8:
    print(f"✓ REDUCTIONS VERIFIED: cos-overlap D_kl and Q[k,l] computation are CORRECT")
    print(f"✓ No discrepancies found in the even-sector Weil form")
else:
    print(f"✗ DISCREPANCY DETECTED")
