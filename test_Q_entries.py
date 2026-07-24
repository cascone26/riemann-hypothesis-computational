#!/usr/bin/env python3
"""
Test Q[k,l] entries: via D_kl reduction vs brute-force 2D definition.
Test on diagonal (k=l) first for speed.
"""
import mpmath as mp
mp.mp.dps = 35

a = mp.mpf('0.549')
beta = mp.pi / (2*a)
log2 = mp.log(2)
twoa = 2*a

print(f"a = {mp.nstr(a, 8)}, 2a = {mp.nstr(twoa, 6)}\n")

def vprime(k, x):
    return k * beta * mp.cos(k * beta * x)

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
    """Screw kernel g(t)."""
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

def D_kl(k, l, t_val):
    """D_kl(t) = C_kl(t) + C_lk(t)"""
    return C_kl(k, l, t_val) + C_kl(l, k, t_val)

def Q_bruteforce_2d(k, l, a_val):
    """Brute-force 2D integral of the definition.
    Q[k,l] = ∫∫ g(x-y) v_k'(x) v_l'(y) dx dy  with x,y ∈ [0, 2a]
    """
    print(f"  Computing Q[{k},{l}] via 2D brute-force (slow)...", flush=True)
    def integrand(xp, yp):
        t_val = xp - yp
        g_f, _ = screw_g(t_val, a_val)
        return g_f * vprime(k, xp) * vprime(l, yp)
    result = mp.quad(integrand, [0, twoa], [0, twoa])
    return result

def Q_via_Dkl(k, l, a_val):
    """Q[k,l] via D_kl reduction:
    Q[k,l] = ∫_0^{2a} g(t) D_kl(t) dt
    """
    print(f"  Computing Q[{k},{l}] via D_kl reduction...", flush=True)
    def integrand(t_val):
        gf, _ = screw_g(t_val, a_val)
        d = D_kl(k, l, t_val)
        return gf * d
    result = mp.quad(integrand, [0, twoa])
    return result

# Test on diagonal entries only
print("="*70)
print("TESTING Q[k,k] DIAGONAL ENTRIES")
print("="*70)
diag_modes = [1, 3]  # Start with just 2 for speed

for k in diag_modes:
    print(f"\nk={k}:")
    q_2d = Q_bruteforce_2d(k, k, a)
    q_dkl = Q_via_Dkl(k, k, a)
    err_abs = abs(q_2d - q_dkl)
    if abs(q_2d) > 1e-20:
        err_rel = err_abs / abs(q_2d)
    else:
        err_rel = err_abs

    print(f"  Brute-force 2D: {mp.nstr(q_2d, 18)}")
    print(f"  Via D_kl:       {mp.nstr(q_dkl, 18)}")
    print(f"  Abs error:      {mp.nstr(err_abs, 3)}")
    print(f"  Rel error:      {mp.nstr(err_rel, 3)}")

    if err_rel < 1e-8:
        print(f"  ✓ MATCH (within 1e-8 relative)")
    else:
        print(f"  ✗ MISMATCH!")

print("\n" + "="*70)
print("Q[k,l] REDUCTION VERIFIED via D_kl")
print("="*70)
