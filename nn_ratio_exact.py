"""
Verify the nearest neighbor ratio reference value for GUE.

The often-cited value 0.5307 (or 4 - 2*sqrt(3) ≈ 0.5359) comes from
the 2x2 Wigner surmise approximation. Let's compute what the actual
GUE prediction is from the Wigner surmise distribution and check
whether our empirical value is closer than we thought.

For GUE Wigner surmise: p(s) = (32/pi^2) * s^2 * exp(-4s^2/pi)
The ratio r = min(s_n, s_{n+1}) / max(s_n, s_{n+1}) has an analytically
computable distribution when consecutive spacings are independent
(they're not exactly independent, but it's a first check).
"""

import math
import random


def gue_wigner_sample(n=10_000_000):
    """
    Sample from GUE Wigner surmise and compute mean ratio.
    p(s) = (32/pi^2) * s^2 * exp(-4s^2/pi)

    This is a Gamma-like distribution. We can sample via rejection.
    Actually, s^2 * exp(-4s^2/pi) looks like it can be sampled:
    Let u = s^2, then u ~ Gamma(3/2, pi/4) -> s = sqrt(u).
    """
    # Gamma(3/2, scale=pi/8) gives us the right distribution for s^2
    # Actually: p(s) ∝ s^2 exp(-4s^2/pi)
    # Let u = 4s^2/pi, then s = sqrt(pi*u/4)
    # p(u) ∝ u * exp(-u) -> Gamma(2, 1)

    samples = []
    for _ in range(n):
        # Sample from Gamma(2, 1): u = -ln(r1*r2) where r1,r2 ~ Uniform(0,1)
        u = -math.log(random.random()) - math.log(random.random())
        s = math.sqrt(math.pi * u / 4)
        samples.append(s)

    # Compute mean ratio from consecutive pairs
    ratios = []
    for i in range(0, len(samples) - 1, 2):
        s1, s2 = samples[i], samples[i + 1]
        if max(s1, s2) > 0:
            ratios.append(min(s1, s2) / max(s1, s2))

    mean_r = sum(ratios) / len(ratios)
    return mean_r


def analytical_mean_ratio():
    """
    Compute <r> analytically for independent GUE Wigner surmise spacings.

    <r> = 2 * integral_0^inf integral_0^s2 (s1/s2) * p(s1) * p(s2) ds1 ds2

    where p(s) = (32/pi^2) * s^2 * exp(-4s^2/pi)
    """
    # Numerical integration
    ds = 0.002
    s_max = 5.0
    c = 32 / math.pi ** 2
    a = 4 / math.pi

    def p(s):
        return c * s ** 2 * math.exp(-a * s ** 2)

    total = 0.0
    for i in range(1, int(s_max / ds)):
        s2 = i * ds
        for j in range(1, i + 1):
            s1 = j * ds
            total += (s1 / s2) * p(s1) * p(s2) * ds * ds

    # Factor of 2 for symmetry
    mean_r = 2 * total
    return mean_r


def main():
    print("NEAREST NEIGHBOR RATIO — GUE REFERENCE VALUES")
    print("=" * 55)

    print("\nCommonly cited values:")
    print(f"  4 - 2*sqrt(3) = {4 - 2*math.sqrt(3):.8f} (2x2 GOE surmise)")
    print(f"  0.5307 (often cited for GUE, source unclear)")

    print("\nGUE Wigner surmise (independent spacings):")
    analytical = analytical_mean_ratio()
    print(f"  Analytical (numerical integration): {analytical:.8f}")

    print("\n  Monte Carlo (10M samples):")
    mc = gue_wigner_sample(10_000_000)
    print(f"  Monte Carlo:                        {mc:.8f}")

    print(f"\nOur empirical value (100K zeros):     0.61091678")
    print(f"Our empirical value (2M zeros):       [pending]")

    print("\nNote: consecutive spacings of zeta zeros are NOT independent.")
    print("The true GUE ratio accounts for correlations between neighbors.")
    print("Our empirical value should match the CORRELATED GUE prediction,")
    print("not the independent surmise prediction.")


if __name__ == "__main__":
    main()
