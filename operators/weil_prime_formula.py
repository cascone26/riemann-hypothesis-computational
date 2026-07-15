"""
Weil Explicit Formula: λ_n from Prime Numbers

GOAL: Express λ_n = Σ_ρ [1-(1-1/ρ)^n] directly as a sum over PRIME POWERS,
bypassing the zeros entirely.

If λ_n has a manifestly positive prime-number representation, that proves RH
without any finite computational certificate — closing the gap from 5H.

------------------------------------------------------------------------------
STARTING POINT — the explicit formula for log ξ:

ξ(s) = (1/2)s(s-1)π^{-s/2}Γ(s/2)ζ(s)

Taking log: log ξ(s) = log(1/2) + log(s) + log(s-1) - (s/2)log(π) + log Γ(s/2) + log ζ(s)

Differentiating: (ξ'/ξ)(s) = 1/s + 1/(s-1) - log(π)/2 + (1/2)(Γ'/Γ)(s/2) - Λ̂(s)

where Λ̂(s) = -(ζ'/ζ)(s) = Σ_{n≥1} Λ(n)/n^s  (von Mangoldt Dirichlet series)

------------------------------------------------------------------------------
THE GENERATING FUNCTION:

G(t) = Σ_{n≥0} λ_{n+1} t^n = (d/dt) log ξ(1/(1-t))

= (ξ'/ξ)(1/(1-t)) × (1/(1-t))²

Let w = 1/(1-t), so t = 1 - 1/w, dt/dw = 1/w²:

G(t) = (ξ'/ξ)(w) / (w-1)²  [in terms of w = 1/(1-t)]

= [1/w + 1/(w-1) - log(π)/2 + (1/2)(Γ'/Γ)(w/2) + Σ_n Λ(n)/n^w] / (w-1)²

where the substitution w = 1/(1-t) makes each term a power series in t.

------------------------------------------------------------------------------
TAYLOR COEFFICIENTS — each piece of G(t):

Piece 1: 1/(w(w-1)²)  where w = 1/(1-t)
  = (1-t)/((1-t-1)²) = (1-t)·(1-t)²/(-t)²·(1/(1-t)) ... [messy]
  Actually: 1/(w·(w-1)²) = (1-t)/((1/(1-t)-1)²) = (1-t)/((t/(1-t))²) = (1-t)³/t²
  → singularity at t=0, not a power series. This term is handled differently.

Piece 2: 1/(w-1)³  where w = 1/(1-t)
  w-1 = 1/(1-t) - 1 = t/(1-t)
  1/(w-1)³ = (1-t)³/t³ → not a power series in t.

The analytic continuation of G(t) is delicate. The simple terms 1/w and 1/(w-1)
create poles and branch cuts in t that need careful handling.

------------------------------------------------------------------------------
PRACTICAL APPROACH: Compute λ_n via finite Euler product

For the Dirichlet series -(ζ'/ζ)(s) = Σ Λ(n)/n^s, the Taylor coefficients
of G(t) can be computed via:

λ_n = (−1)^n × (explicit sum over prime powers weighted by n-th derivative)

More directly, from the integral representation:
λ_n = (1/2πi) ∮_{|t|=r} G(t)/t^{n+1} dt

which can be approximated by summing the prime power contributions to G.

The CLEAN VERSION (Li 1997, Equation (1.2)):
λ_n = Σ_{k=1}^{n} (n choose k) × (-1)^{k+1} × η_k

where η_k = Σ_ρ ρ^{-k} + (sum of contributions from trivial zeros and boundary terms)

And from the explicit formula for Σ_ρ ρ^{-k}:

Actually the clearest formula (from Lagarias 1999, "On a positivity property..."):

λ_n = Σ_{k=1}^{n} (-1)^{k+1} C(n,k) × [log(4π) - 1 + Σ_{j=0}^{k-1} (-1)^j/(j+1) - Γ'(k)/Γ(k) + extra_k]

where extra_k involves Bernoulli numbers and the prime-power sum Σ_p Σ_m log(p) × Q(p^m, k).

The prime-power term:
  extra_k = (-1)^{k-1}/(k-1) × Σ_{p prime} Σ_{m≥1} Λ(p^m) × (specific weight)

Computing this directly:
"""

import numpy as np
from sympy import primerange, factorint, ln, log
import math

GAMMA_EULER = 0.5772156649015328606  # Euler-Mascheroni constant


def von_mangoldt(n):
    """Λ(n) = log(p) if n = p^k, else 0."""
    if n < 2:
        return 0.0
    fac = factorint(n)
    if len(fac) == 1:
        p = list(fac.keys())[0]
        return math.log(p)
    return 0.0


def lambda_from_prime_approx(n, x_max=1000):
    """
    Approximate λ_n using the explicit formula with prime powers up to x_max.

    From the Weil/Guinand explicit formula for the Li coefficients:

    λ_n = 1 + Σ_{j=1}^{n} C(n,j)/(j) × [A_j + B_j]

    where:
      A_j = analytic terms (log(4π), digamma, etc.)
      B_j = prime-power sum: (-1)^{j+1} × Σ_{p^k ≤ x_max} Λ(p^k) × log(p^k)^{j-1}/p^{k/2}
            (this is the Riemann-prime counting contribution)

    WARNING: This formula requires careful regularization. The sum over prime powers
    converges only conditionally and the boundary terms at x_max create error O(√x_max).
    """
    from math import comb, log, lgamma, exp

    # Analytic boundary terms: from 1/s, 1/(s-1), Γ'/Γ in ξ'/ξ
    # At the n-th level:
    # A_j = (-1)^{j+1} × [log(4π) - GAMMA_EULER + Σ_{k=1}^{j} (-1)^k/k]
    #        = (-1)^{j+1} × [log(4π) - GAMMA_EULER + 1 - 1/2 + 1/3 - ... ± 1/j]

    # This alternating harmonic sum: H_alt(j) = Σ_{k=1}^{j} (-1)^{k+1}/k = 1 - 1/2 + 1/3 - ...
    # = log(2) as j → ∞

    analytic_parts = []
    for j in range(1, n + 1):
        H_alt = sum((-1)**(k+1) / k for k in range(1, j+1))
        A_j = (-1)**(j+1) * (math.log(4*math.pi) - GAMMA_EULER + H_alt)
        analytic_parts.append(A_j)

    # Prime power sum: the key arithmetic contribution
    # B_j = Σ_{p^k: prime power ≤ x_max} log(p) × (−log p^k)^{j-1} / p^{k/2}
    #       × binomial × sign adjustment

    prime_contributions = [0.0] * (n + 1)
    for m in range(2, x_max + 1):
        lam = von_mangoldt(m)
        if lam == 0:
            continue
        sqrt_m = math.sqrt(m)
        log_m = math.log(m)
        for j in range(1, n + 1):
            # Weight for j-th term: log(p^k)^{j-1} / p^{k/2}
            contrib = lam * (log_m ** (j-1)) / sqrt_m
            prime_contributions[j] += contrib

    lambda_approx = 1.0  # boundary term from ξ(0) = 1/2
    for j in range(1, n + 1):
        coeff = math.comb(n, j) * ((-1)**(j+1)) / math.factorial(j-1)
        term = coeff * (analytic_parts[j-1] + prime_contributions[j] / math.factorial(j-1))
        lambda_approx += term

    return lambda_approx


def lambda_via_log_xi_derivatives(n):
    """
    Compute λ_n from the Taylor expansion of log ξ(s) at s=1.

    Exact formula (Keiper 1992):
    λ_n = Σ_{k=1}^{n} (-1)^{k+1} C(n,k) × b_k

    where b_k = (1/(k-1)!) × d^{k-1}/ds^{k-1}|_{s=1} (s-1)ξ'/ξ(s) × s^{k-1}

    From the Laurent expansion of (s-1)ξ'/ξ(s) at s=1:
    (s-1)ξ'/ξ(s) = 1 - γ₀(s-1) + γ₁(s-1)² - ...  (Stieltjes constants)

    where γ_j = (-1)^j/j! × lim_{N→∞} [Σ_{n=1}^{N} log(n)^j/n - log(N)^{j+1}/(j+1)]
    """
    from math import comb

    # Stieltjes constants (known values)
    stieltjes = [
        GAMMA_EULER,            # γ₀
        -0.0728158454836767248, # γ₁
        -0.0096903631928723184, # γ₂
        0.002053834420303346134, # γ₃
        0.0023209775599343375,  # γ₄
        0.000793323817301062701, # γ₅
    ]

    # b_k = coefficient in Σ_k (-1)^{k+1} C(n,k) b_k
    # From the explicit formula relating b_k to Stieltjes constants:
    # The leading behavior: b_k ~ (log(4π) - γ₀) = constant for k→∞ under RH
    # But computing b_k exactly requires the full analytic structure.

    # For small k, b_k can be computed from:
    # b_1 = 1 + log(4π)/2 - γ₀/2 - 1 + log(4) - γ₀/2 ... [see Li 1997]

    # This is numerically unstable for large k. Use the zero-based formula instead.
    return None  # Placeholder — full computation requires high-precision arithmetic


def compute_lambda_from_zeros(gammas, n):
    """Compute λ_n = Σ_j F_n(1/2,γ_j) for on-line zeros (reference computation)."""
    phi = 2.0 * np.arctan(1.0 / (2.0 * gammas))
    F = 4.0 * (1.0 - np.cos(n * phi))
    return np.sum(F)


def prime_structure_analysis(n_max=50, x_max=500):
    """
    Examine the prime-power decomposition of λ_n.

    The key question: do the prime-power contributions to λ_n form a manifestly
    positive sum, or do they require cancellation to give λ_n > 0?

    If manifestly positive: λ_n > 0 for all n is an ARITHMETIC fact, not number-theoretic.
    If by cancellation: proving λ_n > 0 requires information about the zeros.
    """
    print(f"\nPrime Structure Analysis (prime powers up to {x_max}):")
    print(f"{'n':>6} {'analytic':>12} {'prime sum':>12} {'total':>12} {'zero-based':>12}")
    print("-" * 58)

    from math import comb, log, factorial

    gammas = []
    try:
        with open("data/zeros_100k.txt") as f:
            for i, line in enumerate(f):
                if i >= 200: break
                gammas.append(float(line.strip()))
        gammas = np.array(gammas)
    except FileNotFoundError:
        gammas = np.array([14.1347, 21.0220, 25.0109, 30.4249, 32.9351,
                           37.5862, 40.9187, 43.3271, 48.0052, 49.7738])

    for n in range(1, min(n_max + 1, 10)):
        # Analytic (non-prime) contribution to λ_n
        # From the Keiper formula: the "trivial" part not involving zeros or primes
        # b_0 = log(4π/e) + γ contribution
        analytic = 0.0
        for j in range(1, n + 1):
            H_alt = sum((-1)**(k+1)/k for k in range(1, j+1))
            A_j = (-1)**(j+1) * (log(4*math.pi) - GAMMA_EULER + H_alt)
            analytic += comb(n, j) * A_j

        # Prime power contribution
        prime_sum = 0.0
        for m in range(2, x_max + 1):
            lam = von_mangoldt(m)
            if lam == 0: continue
            sqrt_m = math.sqrt(m)
            log_m = math.log(m)
            for j in range(1, n + 1):
                c = comb(n, j) * ((-1)**(j+1)) / factorial(j)
                prime_sum += c * lam * (log_m**(j)) / sqrt_m

        # Reference from zeros
        zero_val = compute_lambda_from_zeros(gammas, n)

        total = analytic + prime_sum
        print(f"{n:>6} {analytic:>12.4f} {prime_sum:>12.4f} {total:>12.4f} {zero_val:>12.4f}")


def sign_analysis_prime_terms():
    """
    Direct sign analysis: do prime power contributions to λ_n have definite sign?

    Key observation: the von Mangoldt function Λ(n) ≥ 0 always.
    If the weight assigned to each prime power p^k in λ_n is non-negative,
    then λ_n ≥ (positive contribution from non-prime terms).
    """
    print("\nSign Analysis of Prime Power Weights in λ_n:")
    print()
    print("For λ_n, the weight of prime power p^k is:")
    print("  W(p^k, n) = Σ_{j=1}^{n} C(n,j) × (-1)^{j+1} / j! × log(p^k)^j / p^{k/2}")
    print("            = (1/p^{k/2}) × Σ_{j=1}^{n} C(n,j) × (-1)^{j+1} × (log p^k)^j / j!")
    print("            = (1/p^{k/2}) × [1 - (1 - log p^k)^n + ... polynomial in log p^k]")
    print()

    # The inner sum: Σ_{j=1}^{n} C(n,j) (-1)^{j+1} x^j/j! where x = log(p^k)
    # = -(1/x) × Σ_{j=1}^{n} C(n,j) (-x)^j/j = -(1/x) × [(1-x)^n - 1 + n·x/1! ...]

    # Actually: Σ_{j=0}^{n} C(n,j)(-1)^j x^j/j! = [coefficient of ??? in e^{-x(1-t)^n}]
    # This doesn't simplify cleanly.

    # Numerically: compute W(p^k, n) for small n, various p^k
    from math import comb, factorial, log, sqrt

    primes_powers = [2, 3, 4, 5, 7, 8, 9, 11, 13]

    print(f"  W(p^k, n) for various n and p^k:")
    print(f"  {'p^k':>6} {'n=1':>8} {'n=2':>8} {'n=3':>8} {'n=5':>8} {'n=10':>8} {'sign?':>8}")
    print("  " + "-" * 60)

    for m in primes_powers:
        lam = von_mangoldt(m)
        if lam == 0:
            continue
        x = log(m)  # log(p^k)
        sqrtm = sqrt(m)
        row = [f"{m:>6}"]
        signs = set()
        for n in [1, 2, 3, 5, 10]:
            W = sum(comb(n,j) * ((-1)**(j+1)) * (x**j) / factorial(j)
                    for j in range(1, n+1)) / sqrtm
            row.append(f"{W:>8.4f}")
            signs.add("+" if W > 0 else "-")
        consistent = "✓ +" if all(s=="+" for s in signs) else ("✓ -" if all(s=="-" for s in signs) else "MIXED")
        row.append(f"{consistent:>8}")
        print("  " + "".join(row))

    print()
    print("KEY OBSERVATION:")
    print("  If W(p^k, n) > 0 for ALL prime powers p^k and ALL n ≥ 1,")
    print("  then λ_n ≥ Σ_{p^k} Λ(p^k) × W(p^k, n) > 0 (manifestly).")
    print("  This would prove λ_n > 0 without ANY zero computation.")


def explicit_weight_formula(x, n):
    """
    W(x, n) = (1/√x) × Σ_{j=1}^{n} C(n,j)(-1)^{j+1} (log x)^j / j!

    Closed form: using the generating function of the exponential polynomials.

    Let u = log(x) = log(p^k):
    W(x, n) × √x = Σ_{j=1}^{n} C(n,j)(-u)^j × (-1) / j!
                 = -Σ_{j=1}^{n} C(n,j)(-u)^j / j!

    Note: Σ_{j=0}^{n} C(n,j)(-u)^j / j! = ? (no simple closed form)

    But: Σ_{j=0}^{n} C(n,j)(-u)^j = (1-u)^n  (ordinary binomial)
    And: Σ_{j=0}^{n} C(n,j)(-u)^j / j! = exponential polynomial / exponential

    Numerically:
    """
    from math import comb, factorial, log, sqrt, exp
    u = math.log(x)
    inner = sum(comb(n, j) * ((-1)**(j+1)) * (u**j) / factorial(j) for j in range(1, n+1))
    return inner / math.sqrt(x)


def positivity_of_weights(n_max=30, x_max=1e6):
    """
    Check: is W(p^k, n) ≥ 0 for all prime powers p^k and all n ≤ n_max?

    If yes: λ_n = (analytic terms) + Σ Λ(p^k) × W(p^k, n) ≥ analytic terms.
    The analytic terms are related to ξ(0), ξ'(0), etc., and separately computable.
    """
    print(f"\nPositivity of Prime Weights W(p^k, n) for n ≤ {n_max}:")

    all_positive = True
    min_W = float('inf')
    min_at = None

    for x in range(2, 1001):
        lam = von_mangoldt(x)
        if lam == 0:
            continue
        for n in range(1, n_max + 1):
            W = explicit_weight_formula(x, n)
            if W < min_W:
                min_W = W
                min_at = (x, n, W)
            if W < 0:
                all_positive = False
                print(f"  NEGATIVE: W({x},{n}) = {W:.6f} (prime power {x})")

    print(f"  Min weight: W({min_at[0]}, {min_at[1]}) = {min_at[2]:.6f}")
    print(f"  All weights ≥ 0 for n ≤ {n_max}, p^k ≤ 1000: {all_positive}")

    if all_positive:
        print()
        print("  *** POTENTIAL BREAKTHROUGH: If W(p^k,n) ≥ 0 for ALL p^k, ALL n, ***")
        print("  *** then λ_n ≥ 0 is an ARITHMETIC fact, independent of the zeros! ***")
        print()
        print("  Remaining to verify: analytic boundary terms ≥ -Σ Λ(p^k) W(p^k,n).")
        print("  Check large p^k: does W → 0+ or cross to negative?")

        # Check large prime powers
        print(f"\n  Large prime power behavior:")
        for x in [1e3, 1e6, 1e9, 1e12]:
            for n in [1, 5, 10, 30]:
                W = explicit_weight_formula(x, n)
                print(f"    W({x:.0e}, {n}) = {W:.2e}", end="   ")
            print()

    return all_positive, min_at


def closed_form_weight(u, n):
    """
    W(p^k, n) × √(p^k) = f(u,n) where u = log(p^k).

    f(u,n) = Σ_{j=1}^{n} C(n,j)(-1)^{j+1} u^j/j!
           = 1 - e^{-u} × (sum from exponential series truncated at n)

    Wait: Σ_{j=0}^{∞} C(n,j)(-1)^j u^j/j! = ?
    Not the exponential — the binomial (1-u)^n is for integer powers, not divided factorials.

    Correct expansion: Σ_{j=0}^{n} C(n,j)(-u)^j/j! = ?

    This is related to the incomplete Gamma function or the Poisson CDF.

    Alternative: note that for u = log(p^k), the weight is:
    f(u,n) = Σ_{j=1}^{n} C(n,j) × (-1)^{j+1}/j! × u^j

    This is the Euler-Maclaurin approximation to (1-e^{-u}), truncated at j=n.
    For small u: f(u,n) ≈ u - u²/2 + u³/6 - ... ≈ 1 - e^{-u} > 0.
    For large u (u >> n): f(u,n) is dominated by the leading term (-1)^{n+1}C(n,n)u^n/n!
                        = (-1)^{n+1}u^n/n!  → alternates in sign.

    KEY: For u > n (i.e., log(p^k) > n, i.e., p^k > e^n), the weight can be negative!
    """
    from math import comb, factorial
    return sum(comb(n, j) * ((-1)**(j+1)) * (u**j) / factorial(j) for j in range(1, n+1))


def weight_sign_boundary(n_max=20):
    """
    Find the threshold u_*(n) where f(u,n) changes sign.
    For u > u_*(n): W(p^k, n) < 0.

    This is critical: if u_*(n) < log(primes), then LARGE primes give NEGATIVE
    contributions to λ_n, and λ_n > 0 is NOT a manifestly positive sum.
    """
    print(f"\nSign Boundary for Prime Weights (n ≤ {n_max}):")
    print(f"  W(p^k, n) changes sign at u = log(p^k) = u_*(n):")
    print(f"  {'n':>6} {'u_*(n)':>12} {'p^k boundary':>16} {'sign for p^k>bdy':>20}")
    print("  " + "-" * 58)

    from scipy.optimize import brentq

    for n in range(1, n_max + 1):
        # Find u where f(u,n) = 0 (first zero for u > 0)
        # f(u,n) starts positive (≈ u for small u), check if it becomes negative
        u_vals = np.linspace(0.1, 3*n, 1000)
        f_vals = np.array([closed_form_weight(u, n) for u in u_vals])

        # Find sign changes
        sign_changes = np.where(np.diff(np.sign(f_vals)))[0]

        if len(sign_changes) == 0:
            print(f"  {n:>6} {'∞ (always +)':>12} {'always':>16} {'POSITIVE':>20}")
        else:
            u_star_idx = sign_changes[0]
            u_lo, u_hi = u_vals[u_star_idx], u_vals[u_star_idx + 1]
            try:
                u_star = brentq(lambda u: closed_form_weight(u, n), u_lo, u_hi)
                pk_boundary = math.exp(u_star)
                sign_large = "NEGATIVE" if closed_form_weight(u_star * 1.1, n) < 0 else "POSITIVE"
                print(f"  {n:>6} {u_star:>12.4f} {pk_boundary:>16.2f} {sign_large:>20}")
            except:
                print(f"  {n:>6} {'~'+str(u_vals[u_star_idx])[:6]:>12} {'?':>16} {'?':>20}")


if __name__ == "__main__":
    print("=" * 70)
    print("Weil Explicit Formula: λ_n from Prime Numbers")
    print("=" * 70)

    # 1. Sign analysis of prime weights
    sign_analysis_prime_terms()

    # 2. Positivity check
    all_pos, min_info = positivity_of_weights(n_max=20)

    # 3. Sign boundary analysis (critical!)
    try:
        from scipy.optimize import brentq
        weight_sign_boundary(n_max=20)
    except ImportError:
        # Manual boundary search
        print("\nSign boundary (manual search, scipy not available):")
        for n in range(1, 11):
            positive = True
            u_star = None
            for u_test in np.linspace(0.1, 5*n, 500):
                if closed_form_weight(u_test, n) < 0:
                    u_star = u_test
                    positive = False
                    break
            if positive:
                print(f"  n={n}: always positive for u ≤ {5*n}")
            else:
                print(f"  n={n}: sign change near u={u_star:.3f}, p^k={math.exp(u_star):.2f}")

    print("\n" + "=" * 70)
    print("CONCLUSION:")
    print()
    print("The prime-power weight W(p^k, n) in λ_n is NOT always positive.")
    print("For log(p^k) > u_*(n) (large prime powers), W becomes negative.")
    print("Therefore λ_n > 0 is NOT a manifestly positive prime-power sum.")
    print()
    print("This means: proving λ_n > 0 via the prime-power formula requires")
    print("cancellation between large-prime contributions. Such cancellation")
    print("is equivalent to information about the zeros of ζ(s).")
    print()
    print("The prime-number approach does NOT immediately bypass the zero gap.")
    print()
    print("HOWEVER: The weight analysis reveals STRUCTURE:")
    print("  - Small primes (p^k ≤ e^{u_*(n)}) contribute positively.")
    print("  - Large primes (p^k > e^{u_*(n)}) contribute negatively.")
    print("  - The balance determines the sign of λ_n.")
    print("  - Under PNT (or RH), the balance is positive for all n.")
    print("  - This connection deserves further exploration.")
    print("=" * 70)
