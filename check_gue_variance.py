"""
Check: what IS the actual GUE nearest-neighbor spacing variance?
The Wigner surmise for GUE gives p(s) = (32/pi^2) * s^2 * exp(-4s^2/pi).
Let's compute the exact variance from this distribution.
"""

import math


def wigner_gue_moments():
    """
    Compute moments of the GUE Wigner surmise distribution analytically.
    p(s) = (32/pi^2) * s^2 * exp(-4*s^2/pi)

    <s> = integral of s * p(s) ds from 0 to inf
    <s^2> = integral of s^2 * p(s) ds from 0 to inf
    Var(s) = <s^2> - <s>^2
    """
    # Using known results for integrals of s^n * exp(-a*s^2):
    # integral_0^inf s^n * exp(-a*s^2) ds = Gamma((n+1)/2) / (2 * a^((n+1)/2))

    a = 4 / math.pi
    c = 32 / (math.pi ** 2)

    # <s> = c * integral s^3 * exp(-a*s^2) ds = c * Gamma(2) / (2*a^2) = c / (2*a^2)
    mean = c * math.gamma(2) / (2 * a ** 2)

    # <s^2> = c * integral s^4 * exp(-a*s^2) ds = c * Gamma(5/2) / (2*a^(5/2))
    mean_sq = c * math.gamma(5 / 2) / (2 * a ** 2.5)

    var = mean_sq - mean ** 2

    print("GUE Wigner Surmise:")
    print(f"  <s>   = {mean:.10f}")
    print(f"  <s^2> = {mean_sq:.10f}")
    print(f"  Var   = {var:.10f}")
    print()

    # Also compute numerically via trapezoidal rule
    ds = 0.0001
    s_vals = [i * ds for i in range(1, 100000)]

    num_mean = sum(s * (32 / math.pi ** 2) * s ** 2 * math.exp(-4 * s ** 2 / math.pi) * ds for s in s_vals)
    num_mean_sq = sum(s ** 2 * (32 / math.pi ** 2) * s ** 2 * math.exp(-4 * s ** 2 / math.pi) * ds for s in s_vals)
    num_var = num_mean_sq - num_mean ** 2

    print("Numerical verification:")
    print(f"  <s>   = {num_mean:.10f}")
    print(f"  <s^2> = {num_mean_sq:.10f}")
    print(f"  Var   = {num_var:.10f}")
    print()

    # The exact GUE pair correlation spacing variance (not Wigner surmise)
    # is actually different. Wigner surmise is an approximation.
    # Exact GUE spacing variance = 1 - 4/pi^2 - ... (complicated)
    # Known value is approximately 0.178
    print("Note: Wigner surmise is an APPROXIMATION to exact GUE.")
    print("The exact GUE nearest-neighbor spacing variance is ~0.178")
    print("(See Mehta, Random Matrices, 3rd edition)")
    print()
    print("Our empirical variance: 0.1607")
    print("Exact GUE prediction:  ~0.178")
    print("Wigner surmise:         0.286 was WRONG (that's the GOE value!)")


if __name__ == "__main__":
    wigner_gue_moments()
