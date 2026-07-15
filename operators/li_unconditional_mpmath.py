"""
Compute Li coefficients lambda_n unconditionally using mpmath.
Extends the Stieltjes-constant approach beyond n=200.

Key: avoid the pole of zeta(s) at s=1 by using
f(z) = z * zeta(1+z) which is analytic at z=0 with f(0)=1.

xi(1+z) = (1/2)(1+z) * f(z) * Gamma((1+z)/2) * pi^{-(1+z)/2}
"""
import mpmath
from mpmath import mp, mpf, mpc, gamma, pi, log, zeta, fac, binomial, nsum, inf
import json
import os

# Set high precision
mp.dps = 60  # 60 decimal places

print(f"Precision: {mp.dps} decimal places")

def xi_mod(z):
    """
    Compute xi(1+z) = (1/2)(1+z)*z*Gamma((1+z)/2)*pi^{-(1+z)/2}*zeta(1+z)
    Works for small z ≠ 0.
    At z=0, xi(1) = 1/2.
    """
    s = mpf('1') + z
    if abs(z) < mpf('1e-40'):
        return mpf('1')/2
    # z * zeta(1+z) is analytic at z=0 (handles pole)
    fz = z * zeta(s)
    return mpf('1')/2 * s * fz * gamma(s/2) * pi**(-s/2)

def log_xi_mod(z):
    """log(xi(1+z))"""
    if abs(z) < mpf('1e-40'):
        return log(mpf('1')/2)
    return log(xi_mod(z))

print("\nComputing Li coefficients via Taylor expansion of log xi(1+z)...")
print("Method: mpmath numerical differentiation\n")

# Compute lambda_n for n=1..N_MAX
N_MAX = 500

lambda_vals = []

# Taylor coefficient of log xi(1+z) at z=0, order n
# = (1/n!) * d^n/dz^n [log xi(1+z)]|_{z=0}
# lambda_n = (-1)^{n+1} * n * Taylor_coeff_n

print(f"Computing n=1..{N_MAX}...")

for n in range(1, N_MAX+1):
    # Use mpmath.diff with higher-order derivative
    deriv = mpmath.diff(log_xi_mod, mpf('0'), n=n)
    taylor_coeff = deriv / fac(n)
    lam = (-1)**(n+1) * n * taylor_coeff
    lambda_vals.append(float(lam))

    if n % 50 == 0:
        print(f"  n={n:4d}: lambda_{n} = {float(lam):.10f} ({'POSITIVE' if lam > 0 else 'NEGATIVE!'}, min so far: {min(lambda_vals):.8f})")

print(f"\nDone!")
print(f"All positive for n=1..{N_MAX}: {all(v > 0 for v in lambda_vals)}")
print(f"Min: lambda_{lambda_vals.index(min(lambda_vals))+1} = {min(lambda_vals):.10f}")
print(f"Max: lambda_{lambda_vals.index(max(lambda_vals))+1} = {max(lambda_vals):.10f}")

# Compare with known values at n=1..10
known = {1: 0.023016666057459372, 2: 0.0920295635934749, 10: 2.2714350723887677}
print("\nComparison with known values:")
for n, kn in known.items():
    if n <= N_MAX:
        print(f"  n={n}: computed={lambda_vals[n-1]:.12f}, known={kn:.12f}, "
              f"diff={lambda_vals[n-1]-kn:.3e}")

# Save results
results = {
    'N_max': N_MAX,
    'precision_dps': mp.dps,
    'all_positive': all(v > 0 for v in lambda_vals),
    'min_lambda': min(lambda_vals),
    'min_at_n': lambda_vals.index(min(lambda_vals)) + 1,
    'lambda_vals': lambda_vals,
}
outfile = '/Users/scones/projects/riemann/results/li_unconditional_mpmath.json'
with open(outfile, 'w') as f:
    json.dump(results, f, indent=2)
print(f"\nSaved to {outfile}")
print(f"\n=== UNCONDITIONAL RESULT ===")
print(f"lambda_n > 0 for all n=1..{N_MAX}: {all(v > 0 for v in lambda_vals)}")
