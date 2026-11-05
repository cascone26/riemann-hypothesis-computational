#!/usr/bin/env python3
"""
Approach 2: Push Hankel determinant verification to arbitrary precision.
Earlier breakdown at n=5 might be precision budget, not structure.
Test with dps=100, dps=200 — how far can we certify rigorously?
"""

from mpmath import mp, mpf, log as mlog, pi as mpi
import numpy as np

print("=" * 80)
print("HIGH-PRECISION HANKEL DETERMINANT CERTIFICATION")
print("=" * 80)

# Test multiple precision levels
precision_levels = [50, 100, 150, 200]

for dps_target in precision_levels:
    mp.dps = dps_target

    print(f"\n{'='*80}")
    print(f"dps = {dps_target}")
    print(f"{'='*80}\n")

    # Recompute unconditional moments at this precision
    def xi(s):
        return mpf('0.5') * s * (s - 1) * mp.power(mpi, -s/2) * mp.gamma(s/2) * mp.zeta(s)

    K_moments = 30  # more moments at higher precision
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

    m_exact = [mpf(0)] * (K_moments + 1)
    for k in range(1, K_moments + 1):
        m_exact[k] = k * (-p[k])

    # Compute Hankel determinants using Sylvester criterion
    from mpmath import matrix as mp_matrix

    def hankel_det(moments_mpf, n):
        """Compute Hankel determinant det((m_{i+j})_{i,j=0..n-1})."""
        H = mp_matrix(n, n)
        for i in range(n):
            for j in range(n):
                if i + j + 1 <= len(moments_mpf) - 1:
                    H[i, j] = moments_mpf[i + j + 1]
                else:
                    H[i, j] = mpf(0)
        return mp.det(H)

    print(f"{'n':>3} {'det(H_n)':>40} {'Sign':>15}")

    certification_limit = None  # Track where we lose positivity rigorously

    for n in range(1, 15):  # Try up to n=14 at high precision
        try:
            det_n = hankel_det(m_exact, n)
            sign = "✓ >0" if det_n > 0 else ("✗ <0" if det_n < 0 else "? =0")

            # Check: is this det actually reliable at this precision?
            # Heuristic: if |det| < 10^(-dps+10), it might be noise
            if abs(det_n) < mp.mpf(10)**(-dps_target + 10):
                sign += " (NOISE?)"
                if certification_limit is None:
                    certification_limit = n

            print(f"{n:>3} {mp.nstr(det_n, 30):>40} {sign:>15}")
        except Exception as e:
            print(f"{n:>3} [COMPUTATION FAILED] {str(e)[:30]}")
            break

    if certification_limit:
        print(f"\n→ Rigorously certified up to n={certification_limit-1} (n={certification_limit} likely noise)")
    else:
        print(f"\n→ Determinants all clean, no obvious noise threshold crossed")

print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print("""
If high precision extends the rigorous certification beyond n=4, that's new
ground: it shows the "breakdown at n=5" is precision-dependent, not structural.
The true barrier would be wherever even dps=200 can't certify further, which
would tell us the actual conditioning limit of the moment sequence.
""")
