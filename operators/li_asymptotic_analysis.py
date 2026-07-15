"""
Verify and analyze the asymptotic formula for Li coefficients.

Phi's formula (unconditional main term):
  lambda_n = (n/2)(log n - log(2*pi) - gamma_EM + 1) + R_n
  where R_n = sum_rho [1 - (1-1/rho)^n]

Under RH: R_n = O(n^{1/2} * log^2 n)
Off RH: R_n grows super-polynomially.

This script:
1. Loads computed lambda_n from N=2000 results
2. Computes the main asymptotic term
3. Measures R_n = lambda_n - main_term
4. Checks R_n / n^{1/2} stability (RH predicts bounded ratio)
5. Finds where main term dominates: i.e., when main_term > 2*|R_n|
"""

import json, math, os

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "results")

# Load N=2000 results
with open(os.path.join(RESULTS_DIR, 'li_unconditional_2000.json')) as f:
    data = json.load(f)

lambda_vals = data['lambda_vals']
N = len(lambda_vals)

gamma_EM = 0.5772156649015328606  # Euler-Mascheroni constant

def main_term(n):
    """Keiper-Li asymptotic main term (unconditional)."""
    # lambda_n ~ (n/2) * (log(n) - log(2*pi) - gamma_EM + 1)
    # From Bombieri-Lagarias: the 'trivial' contribution to lambda_n
    return (n / 2) * (math.log(n) - math.log(2 * math.pi) - gamma_EM + 1)

def main_term_v2(n):
    """Alternative main term from Keiper 1992 numerics."""
    # Keiper observed lambda_n / (n/2 * log(n)) -> 1 for large n
    return (n / 2) * math.log(n)

print("Asymptotic Analysis of Li Coefficients (from N=2000 data)")
print("="*70)
print()
print("Testing: lambda_n = (n/2)(log n - log(2*pi) - gamma + 1) + R_n")
print()

print(f"{'n':>8} | {'lambda_n':>14} | {'main_term':>14} | {'R_n':>12} | {'R_n/n^0.5':>12} | {'R_n/main%':>10}")
print("-"*80)

ns = list(range(1, 11)) + [20, 50, 100, 200, 500, 1000, 1500, 2000]
for n in ns:
    if n > N: continue
    lam = lambda_vals[n-1]
    mt = main_term(n)
    R_n = lam - mt
    rn_normalized = R_n / math.sqrt(n) if n > 0 else 0
    rn_pct = 100 * R_n / mt if abs(mt) > 0 else float('inf')
    print(f"{n:>8d} | {lam:>14.6f} | {mt:>14.6f} | {R_n:>12.6f} | {rn_normalized:>12.6f} | {rn_pct:>10.4f}%")

print()
print("R_n / n^{1/2} stability check:")
print("  (Under RH, this should be O(log^2 n) — slowly growing)")
print()

# Check if R_n / sqrt(n) is bounded and growing slowly
rn_data = []
for n in [10, 50, 100, 200, 500, 1000, 2000]:
    if n > N: continue
    lam = lambda_vals[n-1]
    mt = main_term(n)
    R_n = lam - mt
    rn_norm = R_n / math.sqrt(n)
    log2n = math.log(n)**2
    print(f"  n={n:5d}: R_n/sqrt(n) = {rn_norm:10.4f}, log^2(n) = {log2n:.2f}, ratio = {rn_norm/log2n:.4f}")
    rn_data.append((n, rn_norm))

print()
print("Main term dominance (when main_term > 10*|R_n|):")
for n in range(1, N+1):
    lam = lambda_vals[n-1]
    mt = main_term(n)
    R_n = lam - mt
    if abs(mt) > 10 * abs(R_n):
        print(f"  First n where main_term > 10*|R_n|: n = {n}")
        print(f"    lambda_{n} = {lam:.6f}, main_term = {mt:.6f}, R_n = {R_n:.6f}")
        break

for threshold in [2, 5, 10, 100]:
    for n in range(1, N+1):
        lam = lambda_vals[n-1]
        mt = main_term(n)
        R_n = lam - mt
        if mt > 0 and abs(mt) > threshold * abs(R_n):
            print(f"  n where main > {threshold:3d}*|R_n|: n = {n}")
            break

print()
print("Lambda_n / main_term ratio (approaching 1 from above):")
for n in [10, 50, 100, 200, 500, 1000, 1500, 2000]:
    if n > N: continue
    lam = lambda_vals[n-1]
    mt = main_term(n)
    print(f"  n={n:5d}: ratio = {lam/mt:.8f}")

print()
print("KEIPER ASYMPTOTIC FORMULA (simpler form):")
print("  lambda_n ~ (n/2) * log(n/(4*pi)) + (n/2)*(1 - gamma)")
print()
for n in [100, 500, 1000, 2000]:
    if n > N: continue
    lam = lambda_vals[n-1]
    keiper = (n/2) * math.log(n / (4 * math.pi)) + (n/2) * (1 - gamma_EM)
    simple = (n/2) * math.log(n)
    print(f"  n={n}: computed={lam:.4f}, keiper_form={keiper:.4f}, simple=(n/2)log(n)={simple:.4f}")

print()
print("CONCLUSION:")
main_dom_n = None
for n in range(2, N+1):
    lam = lambda_vals[n-1]
    mt = main_term(n)
    R_n = lam - mt
    if mt > 0 and mt > 100 * abs(R_n):
        main_dom_n = n
        break

if main_dom_n:
    print(f"  Main term dominates 100x at n >= {main_dom_n}")
    print(f"  For n >= {main_dom_n}, lambda_n > 0 follows from main term alone IF we trust it")
    print(f"  BUT: trusting the main term requires knowing |R_n| < main_term/100")
    print(f"       which requires bounding the zero sum R_n unconditionally => equivalent to RH")
else:
    print("  Main term never dominates 100x within n=1..2000 (or data insufficient)")
    # Check max dominance
    max_ratio = max(main_term(n) / abs(lambda_vals[n-1] - main_term(n))
                   for n in range(2, N+1)
                   if abs(lambda_vals[n-1] - main_term(n)) > 0)
    print(f"  Maximum main/|R_n| ratio in data: {max_ratio:.1f}x at n ~", end="")
    for n in range(N, 0, -1):
        mt = main_term(n)
        R_n = lambda_vals[n-1] - mt
        if abs(R_n) > 0 and mt / abs(R_n) > max_ratio * 0.99:
            print(f" {n}")
            break
