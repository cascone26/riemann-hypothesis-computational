#!/usr/bin/env python3
"""
ADVERSARIAL VERIFICATION — independent re-derivation of the even-sector Weil form.
Lane: verify the cos-overlap reduction D_kl(t) and entry engine against first principles.

Task:
(1) From DEFINITION Q[k,l] = ∫∫ g(x-y) v_k'(x) v_l'(y) dx dy, INDEPENDENTLY re-derive:
    - the cos-overlap D_kl(t) reduction
    - the D_kl closed form (handed from dkl_closedform.py)
    - the elementary entry integrals
(2) Numerically verify against existing code:
    - dkl_closedform.py (D_kl closed form)
    - entry_engine.py (Q[k,l] entries via closed form)
    - brute-force 2D definition integral
(3) Report any DISCREPANCY.
"""
import sys
sys.path.insert(0, '/Users/scones/projects/riemann/operators/dbn_dynamics')
sys.path.insert(0, '/Users/scones/projects/riemann/flint_venv/lib/python3.11/site-packages')

import mpmath as mp
mp.mp.dps = 40

# ============================================================================
# PART 0: Setup — parameters, basis functions, screw kernel
# ============================================================================

a = mp.mpf('0.549')
beta = mp.pi / (2*a)  # β = π/(2a)
log2 = mp.log(2)
twoa = 2*a

print(f"== ADVERSARIAL VERIFICATION: Even-Sector Weil Form ==")
print(f"a = {mp.nstr(a, 8)}, β = π/(2a) = {mp.nstr(beta, 8)}")
print(f"log 2 = {mp.nstr(log2, 8)}, 2a = {mp.nstr(twoa, 8)}")

# Basis: v_k(x) = sin(kπ(x+a)/(2a)) for x ∈ [-a, a]
# Derivative: v_k'(x) = (kπ/(2a)) cos(kπ(x+a)/(2a)) = kβ cos(kβ(x+a))
# Using x' = x+a ∈ [0, 2a]: v_k'(x') = kβ cos(kβ x')

def vprime(k, x):
    """Derivative v_k'(x) with x ∈ [0, 2a] (after shift by +a)."""
    return k * beta * mp.cos(k * beta * x)

# Screw kernel g(t) for t ≥ 0 (we exploit evenness)
def prime_data_local(a_val):
    """Prime powers n <= e^{2a}: return [(Lambda(n), log n, sqrt(n))]."""
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
    """Full screw kernel g(t), t ≥ 0. Returns (g_full, g_arch) where g_full includes primes."""
    at = abs(t)

    # Boundary: -4(e^{t/2} + e^{-t/2} - 2)
    boundary = -4*(mp.e**(t/2) + mp.e**(-t/2) - 2)

    # Arch core: -(|t|/2)(ψ(1/4) - log π)
    psi_quarter = mp.digamma(mp.mpf(1)/4)
    logpi = mp.log(mp.pi)
    arch_core = -(at/2)*(psi_quarter - logpi)

    # Arch tail: -(1/4)[ζ(2,1/4) - e^{-|t|/2} Φ(e^{-2|t|}, 2, 1/4)]
    zeta_quarter = mp.zeta(2, mp.mpf(1)/4)
    if at < mp.mpf(10)**(-20):
        lerch_t = zeta_quarter
    else:
        z = mp.e**(-2*at)
        # Φ(z, 2, a) = sum_{n=0}^∞ z^n / (n+a)^2
        lerch_t = mp.nsum(lambda n: z**n / (n + mp.mpf(1)/4)**2, [0, mp.inf])

    arch_tail = -(mp.mpf(1)/4) * (zeta_quarter - mp.e**(-at/2) * lerch_t)
    g_arch = boundary + arch_core + arch_tail

    # Prime sum: Σ_{log n <= |t|} Lambda(n)/sqrt(n) (|t| - log n)
    prime_sum = mp.mpf(0)
    pl = prime_data_local(a_val)
    for (lam, logn, sqrtn) in pl:
        if logn <= at:
            prime_sum += lam/sqrtn * (at - logn)

    g_full = g_arch + prime_sum
    return g_full, g_arch

# Test screw kernel
print(f"\n[KERNEL TEST] g(t) at sample points:")
for t_val in [0.1, 0.5, 1.0, mp.log(2), 1.5]:
    gf, ga = screw_g(t_val, a)
    print(f"  t={mp.nstr(t_val, 6):8s}: g_full={mp.nstr(gf, 10):16s} g_arch={mp.nstr(ga, 10):16s}")

# ============================================================================
# PART 1: INDEPENDENT DERIVATION — D_kl(t) by direct numerical integration
# ============================================================================

print(f"\n== PART 1: D_kl(t) REDUCTION — Direct 2D Integration ==")

def D_kl_definition(k, l, t_val):
    """
    Compute D_kl(t) by direct integration from the definition.
    D_kl(t) = ∫_0^{2a} v_k'(x') v_l'(x' - t) dx'
    where x' ∈ [0, 2a].

    The integral is over x' such that both x' and (x'-t) are in [0, 2a].
    For fixed t ≥ 0, that means x' ∈ [max(0, t), 2a].
    """
    if t_val < 0:
        t_val = -t_val  # D_kl is even in t

    def integrand(xp):
        """v_k'(x') * v_l'(x'-t)"""
        y_shifted = xp - t_val
        return vprime(k, xp) * vprime(l, y_shifted)

    # Bounds: x' ∈ [max(0,t), 2a]
    lower = max(mp.mpf(0), t_val)
    upper = twoa

    if lower >= upper:
        return mp.mpf(0)

    result = mp.quad(integrand, [lower, upper])
    return result

# Test D_kl for a few (k,l) and t values
print(f"Testing D_kl(t) by direct integration from definition:")
test_modes = [(1, 1), (1, 3), (3, 3), (3, 5), (5, 7), (7, 7)]
test_t = [mp.mpf('0.2'), mp.mpf('0.5'), mp.mpf('1.0')]

dkl_results = {}
for k, l in test_modes:
    dkl_results[(k, l)] = {}
    for t_val in test_t:
        d = D_kl_definition(k, l, t_val)
        dkl_results[(k, l)][t_val] = d
        print(f"  D[{k},{l}]({mp.nstr(t_val, 3)}): {mp.nstr(d, 12)}")

# ============================================================================
# PART 2: COMPARE against CLOSED-FORM formula
# ============================================================================

print(f"\n== PART 2: D_kl CLOSED FORM Verification ==")

def D_kl_closedform(k, l, t_val):
    """
    Closed-form D_kl(t) from dkl_closedform.py.
    k ≠ l: D_kl = klβ { [sin(lβt) - sin(kβt)]/(k-l) - [sin(lβt) + sin(kβt)]/(k+l) }
    k = l: D_kk = k²β²(2a - t) cos(kβt) - kβ sin(kβt)
    """
    if k == l:
        return k*k*beta*beta*(twoa - t_val)*mp.cos(k*beta*t_val) - k*beta*mp.sin(k*beta*t_val)
    else:
        return k*l*beta * (
            (mp.sin(l*beta*t_val) - mp.sin(k*beta*t_val))/(k-l)
            - (mp.sin(l*beta*t_val) + mp.sin(k*beta*t_val))/(k+l)
        )

print(f"Comparing D_kl: definition (direct integral) vs closed form:")
print(f"{'k,l':<8} {'t':<8} {'definition':<20} {'closed-form':<20} {'error':<12}")
print("-" * 75)

max_dkl_err = mp.mpf(0)
for k, l in test_modes:
    for t_val in test_t:
        d_def = dkl_results[(k, l)][t_val]
        d_cf = D_kl_closedform(k, l, t_val)
        err = abs(d_def - d_cf)
        max_dkl_err = max(max_dkl_err, err)

        print(f"{str(k)+','+str(l):<8} {mp.nstr(t_val, 3):<8} {mp.nstr(d_def, 18):<20} "
              f"{mp.nstr(d_cf, 18):<20} {mp.nstr(err, 3):<12}")

print(f"\nMax D_kl error: {mp.nstr(max_dkl_err, 3)}")
if max_dkl_err < 1e-10:
    print("✓ D_kl CLOSED FORM VERIFIED — matches definition to high precision")
else:
    print(f"✗ WARNING: D_kl closed form has significant error {mp.nstr(max_dkl_err, 3)}")

# ============================================================================
# PART 3: INDEPENDENT COMPUTATION OF Q[k,l] FROM DEFINITION
# ============================================================================

print(f"\n== PART 3: Q[k,l] ENTRY COMPUTATION ==")

def Q_definition_bruteforce(k, l, a_val):
    """
    Brute-force 2D integration of the definition.
    Q[k,l] = ∫∫_{[-a,a]²} g(x-y) v_k'(x) v_l'(y) dx dy

    Using change of variables: x' = x+a, y' = y+a, so x',y' ∈ [0, 2a].
    Then v_k'(x) becomes v_k'(x'-a) = kβ cos(kβx'), etc.
    And we integrate ∫∫_{[0,2a]²} g(x'-y') v_k'(x') v_l'(y') dx' dy'.
    """
    def integrand(xp, yp):
        t_val = xp - yp
        g_f, _ = screw_g(t_val, a_val)
        return g_f * vprime(k, xp) * vprime(l, yp)

    # 2D quadrature over [0, 2a] × [0, 2a]
    result = mp.quad(integrand, [0, twoa], [0, twoa])
    return result

print(f"Computing Q[k,l] by brute-force 2D integration (slow, high-precision reference):")
modes_even = [1, 3, 5, 7]  # EVEN sector (odd k modes)
Q_bruteforce = {}

for k in modes_even:
    for l in modes_even:
        if l < k:
            continue  # symmetric, so skip lower triangle for now
        q = Q_definition_bruteforce(k, l, a)
        Q_bruteforce[(k, l)] = q
        print(f"  Q[{k},{l}] = {mp.nstr(q, 15)}")

# ============================================================================
# PART 4: COMPARE against entry_engine.py (closed-form via D_kl)
# ============================================================================

print(f"\n== PART 4: Closed-Form Entry Engine (via D_kl) ==")

def Q_via_Dkl(k, l, a_val):
    """
    Compute Q[k,l] via D_kl reduction:
    Q[k,l] = ∫_0^{2a} g(t) D_kl(t) dt
    where we split: ∫ = ∫_0^{log2} + ∫_{log2}^{2a}
    """
    def integrand(t_val):
        gf, _ = screw_g(t_val, a_val)
        d = D_kl_closedform(k, l, t_val)
        return gf * d

    # Integrate over [0, 2a]
    result = mp.quad(integrand, [0, twoa])
    return result

print(f"Computing Q[k,l] using closed-form D_kl and numerical g(t) integration:")
Q_via_dkl = {}

for k in modes_even:
    for l in modes_even:
        if l < k:
            continue
        q = Q_via_Dkl(k, l, a)
        Q_via_dkl[(k, l)] = q
        print(f"  Q[{k},{l}] = {mp.nstr(q, 15)}")

print(f"\n[COMPARISON] Brute-force vs via D_kl:")
print(f"{'k,l':<8} {'brute-force':<20} {'via D_kl':<20} {'error':<12}")
print("-" * 65)

max_Q_dkl_err = mp.mpf(0)
for k in modes_even:
    for l in modes_even:
        if l < k:
            continue
        qb = Q_bruteforce[(k, l)]
        qd = Q_via_dkl[(k, l)]
        err = abs(qb - qd)
        max_Q_dkl_err = max(max_Q_dkl_err, err)

        print(f"{str(k)+','+str(l):<8} {mp.nstr(qb, 18):<20} {mp.nstr(qd, 18):<20} {mp.nstr(err, 3):<12}")

print(f"\nMax Q entry error (brute-force vs D_kl): {mp.nstr(max_Q_dkl_err, 3)}")

# ============================================================================
# PART 5: COMPARE against existing entry_engine.py and weil_matrices
# ============================================================================

print(f"\n== PART 5: Compare Against Existing Code ==")

# Import existing code
try:
    from mpmath_weil_hiprec import weil_matrices
    from entry_engine import Q as Q_engine

    print(f"Loading weil_matrices (mpmath_weil_hiprec.py, K=24, M=1200)...")
    Qf, Qa, Mass, Vv, xs, w = weil_matrices(a, 24, 1200)

    print(f"Extracting even-sector 4×4 block...")
    # Even modes are k=1,3,5,7, which have indices 0,2,4,6 in the full basis (odd indices)
    ev_indices = [i for i in range(24) if (i+1) % 2 == 1][:4]  # k=1,3,5,7

    print(f"[COMPARISON] Brute-force vs entry_engine vs weil_matrices:")
    print(f"{'k,l':<8} {'brute-force':<20} {'entry_engine':<20} {'weil_matrix':<20} {'err vs BF':<12}")
    print("-" * 88)

    max_engine_err = mp.mpf(0)
    max_matrix_err = mp.mpf(0)

    for ik, k in enumerate(modes_even):
        for il, l in enumerate(modes_even):
            if il < ik:
                continue

            qb = Q_bruteforce[(k, l)]
            qe = Q_engine(k, l)
            qm = Qf[ev_indices[ik], ev_indices[il]]

            err_engine = abs(qb - qe)
            err_matrix = abs(qb - qm)

            max_engine_err = max(max_engine_err, err_engine)
            max_matrix_err = max(max_matrix_err, err_matrix)

            print(f"{str(k)+','+str(l):<8} {mp.nstr(qb, 18):<20} {mp.nstr(qe, 18):<20} "
                  f"{mp.nstr(qm, 18):<20} {mp.nstr(err_engine, 3):<12}")

    print(f"\nMax error vs brute-force:")
    print(f"  entry_engine: {mp.nstr(max_engine_err, 3)}")
    print(f"  weil_matrices: {mp.nstr(max_matrix_err, 3)}")

except Exception as e:
    print(f"Could not import/run existing code: {e}")
    import traceback
    traceback.print_exc()

# ============================================================================
# SUMMARY
# ============================================================================

print(f"\n" + "="*80)
print(f"SUMMARY OF ADVERSARIAL VERIFICATION")
print(f"="*80)
print(f"D_kl closed form vs definition:    max error = {mp.nstr(max_dkl_err, 3)}")
print(f"Q via D_kl vs brute-force 2D:      max error = {mp.nstr(max_Q_dkl_err, 3)}")
print(f"\nConclusion:")
if max_dkl_err < 1e-10 and max_Q_dkl_err < 1e-10:
    print(f"✓ ALL REDUCTIONS VERIFIED: D_kl closed form and Q[k,l] computation are CORRECT")
    print(f"✓ No discrepancies found in the even-sector Weil form structure")
else:
    print(f"✗ DISCREPANCY DETECTED")
    if max_dkl_err > 1e-10:
        print(f"  - D_kl closed form error is significant: {mp.nstr(max_dkl_err, 3)}")
    if max_Q_dkl_err > 1e-10:
        print(f"  - Q[k,l] integration via D_kl has error: {mp.nstr(max_Q_dkl_err, 3)}")
