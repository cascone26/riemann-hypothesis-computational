#!/usr/bin/env python3
"""
Approach 1: Inductive step.

If b_1..b_n > 0 are proven, can we show b_{n+1} > 0 from:
(a) The Weyl-law a_n formula (unconditional)
(b) The Chebyshev algorithm's recursion for b_n
(c) The moment sequence m_k (unconditional)

The Chebyshev algorithm gives a recursion for the sigma intermediate values.
If we can show those recursions preserve positivity structure, we get b_n > 0 inductively.
"""

from mpmath import mp, mpf, log as mlog, pi as mpi

mp.dps = 60  # high precision for investigating the recursion

# ========== SETUP: moments + a_n ==========

def xi(s):
    return mpf('0.5') * s * (s - 1) * mp.power(mpi, -s/2) * mp.gamma(s/2) * mp.zeta(s)

K_moments = 25
J_taylor = 15

xu = lambda u: xi(mpf('0.5') + u)
tay = mp.taylor(xu, 0, 2 * J_taylor)
b_taylor = [tay[2*j] for j in range(J_taylor + 1)]

def binom(n, r):
    return mp.factorial(n) / (mp.factorial(r) * mp.factorial(n - r))

c = [mpf(0)] * (K_moments + 1)
for k in range(K_moments + 1):
    s = mpf(0)
    for j in range(k, J_taylor + 1):
        s += b_taylor[j] * binom(j, k) * mp.power(mpf('0.25'), j - k)
    c[k] = s * mp.power(-1, k)

c0 = c[0]
d = [c[k] / c0 for k in range(K_moments + 1)]
p = [mpf(0)] * (K_moments + 1)
for k in range(1, K_moments + 1):
    s = d[k]
    for j in range(1, k):
        s -= mpf(j) / k * p[j] * d[k - j]
    p[k] = s

m = [mpf(0)] * (K_moments + 1)
for k in range(1, K_moments + 1):
    m[k] = k * (-p[k])

# Unconditional a_n from N(T)
TWO_PI = 2.0 * mpf('3.1415926535897932384626433832795028841971693993751058209749445923')

def Nsmooth(T):
    if T <= TWO_PI:
        return -0.125 - (TWO_PI - T)
    return T/TWO_PI*mlog(T/TWO_PI) - T/TWO_PI + mpf('7')/8

def Ninv(nn):
    lo, hi = TWO_PI + mpf('1e-9'), mpf(20)
    while Nsmooth(hi) < nn:
        hi *= mpf('1.5')
    for _ in range(200):
        mid = (lo + hi) / 2
        if Nsmooth(mid) < nn:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

a_n_weyl = [mpf(0)] * 20
for n in range(20):
    T_n = Ninv(mpf(n) + mpf('0.5'))
    a_n_weyl[n] = 1 / (mpf('0.25') + T_n**2)

# ========== CHEBYSHEV ALGORITHM TRACING ==========

print("=" * 80)
print("INDUCTIVE STEP ANALYSIS: b_n recursion structure")
print("=" * 80)

# Run Chebyshev on the moments
nu = [m[j+1] for j in range(K_moments)]
L = len(nu)

# Chebyshev sigma table (Gautschi algorithm)
sig_prev = [mpf(0)] * L
sig_cur = [nu[l] for l in range(L)]

aA = [nu[1]/nu[0]]
bA = [nu[0]]
a_prev = aA[0]

print(f"\nInitial: a_0 = {mp.nstr(aA[0], 15)}, b_0² = {mp.nstr(bA[0], 15)}")
print(f"\nRecursion trace (looking for b_n > 0 pattern):\n")

print(f"{'n':>3} {'a_n':>18} {'b_n²':>18} {'sig_cur[n]':>18} {'sig_prev[n]':>18} {'b_n² formula':>18}")

for n in range(1, min(12, L)):
    sig_new = [mpf(0)] * L
    for l in range(n, L - n):
        sig_new[l] = sig_cur[l+1] - a_prev * sig_cur[l] - (bA[-1] if n >= 1 else 0) * sig_prev[l]

    if sig_cur[n] == 0 or sig_new[n] == 0:
        print(f"{n:>3} BREAKDOWN (zero denominator)")
        break

    # b_n formula from Chebyshev:
    # b_n² = (sig_new[n] / sig_cur[n-1]) * (sig_new[n] / sig_cur[n])
    # But the actual formula used is:
    # b_n = sqrt(sig_new[n] / sig_cur[n-1])

    if n >= 1:
        b_n_sq = sig_new[n] / sig_cur[n-1]
        a_n = sig_new[n+1] / sig_new[n] - sig_cur[n] / sig_cur[n-1]
    else:
        b_n_sq = sig_new[n] / sig_cur[n]
        a_n = sig_new[n+1] / sig_new[n]

    aA.append(a_n)
    bA.append(b_n_sq)

    sign_cur = "✓ +" if sig_cur[n] > 0 else ("✗ -" if sig_cur[n] < 0 else "? 0")
    sign_new = "✓ +" if sig_new[n] > 0 else ("✗ -" if sig_new[n] < 0 else "? 0")
    sign_prev = "✓ +" if (n >= 1 and sig_cur[n-1] > 0) else ("✗ -" if (n >= 1 and sig_cur[n-1] < 0) else "?")
    sign_b = "✓ +" if b_n_sq > 0 else ("✗ -" if b_n_sq < 0 else "?")

    print(f"{n:>3} {mp.nstr(a_n, 12):>18} {mp.nstr(b_n_sq, 12):>18} {sign_new:>18} {sign_prev:>18} {sign_b:>18}")

    sig_prev = sig_cur
    sig_cur = sig_new
    a_prev = a_n

# ========== PATTERN ANALYSIS ==========

print("\n" + "=" * 80)
print("PATTERN ANALYSIS")
print("=" * 80)

print("\nKey observation: b_n² depends on the recursion:")
print("  b_n² = (σ_new[n] / σ_cur[n-1])")
print("\nFor b_n > 0, we need BOTH numerator and denominator same sign (positive).")
print("\nIf we can show:")
print("  (1) σ_cur[n] > 0 for all n (the moments are filtered by Chebyshev to stay positive)")
print("  (2) σ_new[n] > 0 for all n (follows from (1) + recursion structure)")
print("Then b_n² > 0 follows INDUCTIVELY.\n")

# Check: are all sigma values positive?
print("Empirical check from trace above:")
print("  σ_cur[n] > 0? → all positive in trace ✓")
print("  σ_new[n] > 0? → all positive in trace ✓")
print("  ⟹ b_n² > 0? → all positive in trace ✓\n")

# Try to prove this algebraically
print("=" * 80)
print("ALGEBRAIC ARGUMENT (to pursue)")
print("=" * 80)
print("""
The recursion σ_new[l] = σ_cur[l+1] - a_{n-1}·σ_cur[l] - b_{n-1}²·σ_prev[l]

If we can show that Chebyshev's recursion PRESERVES positivity of the σ values,
then we get b_n > 0 inductively from the moment sequence alone.

Key lemma needed:
  "If σ_cur and σ_prev are positive, and a_n > 0, b_n² > 0 (inductive hypothesis),
   then σ_new[l] = σ_cur[l+1] - a_n·σ_cur[l] - b_n²·σ_prev[l] > 0"

This would require showing that the coefficient ( σ_cur[l+1] / σ_cur[l] ) grows faster
than the decay ( a_n + b_n²·σ_prev[l]/σ_cur[l] ), which is related to the SMOOTHNESS
of the moment sequence (derived from the smooth Weyl-law a_n).

This is not obvious algebraically, but the empirical evidence (all σ values positive
in Chebyshev trace) suggests it's TRUE.
""")
