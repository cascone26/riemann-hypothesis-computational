"""
CLEAN RIGOROUS PROOF: F_n(0, gamma) > 0 for all n <= 85, gamma >= gamma_1

The correct approach — not sinh < |sin|, but the exact boundary analysis.

THEOREM: cosh(x)*cos(y) < 1 for 0 < x < y < 2*pi - arctan(sinh(x)).

PROOF STRUCTURE:

Step A (Log inequality): log(r) < arctan(1/gamma) for all gamma > 0.
    Equivalently: (1/2)*log(1+t^2) < arctan(t) for all t > 0.
    Proof: Let f(t) = arctan(t) - (1/2)*log(1+t^2). Then f(0)=0 and
           f'(t) = 1/(1+t^2) - t/(1+t^2) = (1-t)/(1+t^2).
    Wait, that's wrong. Let me redo:
           f'(t) = 1/(1+t^2) - t/(1+t^2) = (1-t)/(1+t^2).
    Hmm f'(t) > 0 for t < 1, negative for t > 1? That can't be right since arctan(t) ~ pi/2 and log grows.

    Actually: d/dt[(1/2)*log(1+t^2)] = t/(1+t^2), and d/dt[arctan(t)] = 1/(1+t^2).
    So f'(t) = 1/(1+t^2) - t/(1+t^2) = (1-t)/(1+t^2).

    At t=1: f(1) = pi/4 - (1/2)*log(2) = 0.7854 - 0.3466 = 0.439 > 0. f is increasing up to t=1.
    At t=2: f(2) = arctan(2) - (1/2)*log(5) = 1.1071 - 0.8047 = 0.302 > 0.

    But f'(t) < 0 for t > 1, so f is decreasing. Does it ever hit 0?
    As t -> inf: arctan(t) -> pi/2 ~ 1.5708, (1/2)*log(1+t^2) ~ log(t) -> inf.

    So actually (1/2)*log(1+t^2) > arctan(t) for large t! The inequality reverses!

    Let me recheck: for t = 10: arctan(10) = 1.4711, (1/2)*log(101) = 2.307. So log > arctan.
    For t = 0.1: arctan(0.1) = 0.0997, (1/2)*log(1.01) = 0.00499. arctan > log.

    For t = 1/gamma where gamma >= gamma_1 = 14.13: t <= 1/14.13 = 0.0708.
    At t = 0.0708: arctan(0.0708) = 0.07063, (1/2)*log(1.005) = 0.002496.
    arctan >> log here.

    More carefully: for t in (0,1), arctan(t) > (1/2)*log(1+t^2)?
    f'(t) = (1-t)/(1+t^2) > 0 for t < 1.
    f(0) = 0.
    So f(t) > 0 for t in (0,1). YES! The inequality holds for t in (0,1). ✓

    For gamma >= gamma_1 > 14: t = 1/gamma < 1/14 < 1. So the inequality holds!
    arctan(1/gamma) > (1/2)*log(1+1/gamma^2) = log(r).

    => x = n*log(r) < n*arctan(1/gamma) = y. ✓

Step B (arctan-sinh inequality): arctan(sinh(x)) < x for all x > 0.
    Let g(x) = arctan(sinh(x)) - x.
    g(0) = 0.
    g'(x) = cosh(x)/(1+sinh^2(x)) - 1 = cosh(x)/cosh^2(x) - 1 = sech(x) - 1.
    For x > 0: sech(x) < 1, so g'(x) < 0.
    Therefore g(x) < g(0) = 0 for all x > 0. ✓

Step C (Characterization of cosh(x)*cos(y) < 1):
    For x > 0, the equation cosh(x)*cos(y) = 1 has solutions in (0, 2*pi):
    y_* = arctan(sinh(x)) and y_** = 2*pi - arctan(sinh(x)).
    [Verification: cos(y_*) = cos(arctan(sinh(x))) = 1/sqrt(1+sinh^2(x)) = 1/cosh(x) = sech(x).
     So cosh(x)*cos(y_*) = cosh(x)*sech(x) = 1. ✓
     Similarly cos(y_**) = cos(2*pi - y_*) = cos(y_*) = sech(x). ✓]

    For y in (y_*, y_**): cos(y) < cos(y_*) = sech(x), so cosh(x)*cos(y) < 1. ✓
    [cos is decreasing on (0, pi) and symmetric about pi, so cos(y) < sech(x) = cos(y_*) for y in (y_*, 2*pi - y_*).]
    For y in (0, y_*): cos(y) > sech(x), so cosh(x)*cos(y) > 1.
    For y in (y_**, 2*pi): cos(y) > sech(x), so cosh(x)*cos(y) > 1.

Step D (Applying to our problem):
    Lower bound: From Steps A and B: y = n*arctan(1/gamma) > n*log(r) = x > arctan(sinh(x)) = y_*.
    So y > y_*. ✓

    Upper bound: We need y < y_** = 2*pi - arctan(sinh(x)).
    Equivalently: y + arctan(sinh(x)) < 2*pi.
    i.e., n*arctan(1/gamma) + arctan(sinh(n*log(r))) < 2*pi.

    Since both arctan(1/gamma) and log(r) are decreasing in gamma, the LHS is maximized at gamma = gamma_1.

    So we need: n*arctan(1/gamma_1) + arctan(sinh(n*log(r_1))) < 2*pi for n = 1..85.
    This is a FINITE numerical check at gamma = gamma_1!

Step E (Monotonicity in gamma):
    Let h(gamma) = n*arctan(1/gamma) + arctan(sinh(n*log(r(gamma)))).
    d(arctan(1/gamma))/dgamma = -1/(gamma^2*(1+1/gamma^2)) = -1/(gamma^2+1) < 0.
    d(log(r))/dgamma = d[(1/2)*log(1+1/gamma^2)]/dgamma = -1/(gamma^3+gamma) < 0.
    So both terms decrease as gamma increases. h(gamma) is decreasing in gamma. ✓
    max h = h(gamma_1).
"""

import numpy as np
from mpmath import mp, mpf, sqrt, log, cos, cosh, sin, sinh, atan, pi, nstr, acos, acosh, asin, asinh

mp.dps = 50

GAMMA_1 = mpf('14.134725141734693790')

print("=" * 70)
print("CLEAN RIGOROUS PROOF: cosh(x)*cos(y) < 1 for n=1..85")
print("=" * 70)

# =============================================================================
# Verify Step A: log(r) < arctan(1/gamma) for t = 1/gamma in (0,1)
# =============================================================================
print("\n--- Step A: log(r) < arctan(1/gamma) for all gamma >= gamma_1 ---")
print()
print("  (1/2)*log(1+t^2) < arctan(t) for all t in (0,1)")
print("  Proof: f(t) = arctan(t) - (1/2)*log(1+t^2) satisfies f(0)=0, f'(t)=(1-t)/(1+t^2) > 0 for t<1.")
print("  => f(t) > 0 for t in (0,1). Since 1/gamma < 1 for gamma > 1, log(r) < arctan(1/gamma). ✓")
print()

# Numerical verification
t = mpf('1') / GAMMA_1
log_r = (mpf('1')/2) * log(1 + t**2)
arctan_t = atan(t)
print(f"  At gamma_1: t = 1/gamma_1 = {float(t):.8f}")
print(f"  log(r) = {float(log_r):.8f}")
print(f"  arctan(1/gamma) = {float(arctan_t):.8f}")
print(f"  log(r) < arctan(1/gamma): {float(log_r) < float(arctan_t)}")

# =============================================================================
# Verify Step B: arctan(sinh(x)) < x for x > 0
# =============================================================================
print("\n--- Step B: arctan(sinh(x)) < x for all x > 0 ---")
print()
print("  g(x) = arctan(sinh(x)) - x satisfies g(0)=0, g'(x) = sech(x) - 1 < 0.")
print("  => g(x) < 0 for all x > 0. ✓")
print()

# Numerical verification for several x values
for x_val in [0.01, 0.1, 0.2122, 0.5, 1.0]:
    x = mpf(str(x_val))
    lhs = atan(sinh(x))
    print(f"  x = {x_val:.4f}: arctan(sinh(x)) = {float(lhs):.8f} < x = {x_val:.8f}: {float(lhs) < x_val}")

# =============================================================================
# Step C: Exact characterization
# =============================================================================
print("\n--- Step C: cosh(x)*cos(y) < 1 iff y in (y_*, 2*pi - y_*) ---")
print("  where y_* = arctan(sinh(x))")
print()
print("  Verification: at y_* = arctan(sinh(x)), cos(y_*) = sech(x).")
for x_val in [0.1, 0.2122]:
    x = mpf(str(x_val))
    y_star = atan(sinh(x))
    cos_y_star = cos(y_star)
    sech_x = 1/cosh(x)
    prod = float(cosh(x) * cos(y_star))
    print(f"  x={x_val:.4f}: y_*={float(y_star):.6f}, cos(y_*)={float(cos_y_star):.8f}, sech(x)={float(sech_x):.8f}, cosh*cos={prod:.8f}")

# =============================================================================
# Step D: The key condition — monotone in gamma, tightest at gamma_1
# =============================================================================
print("\n--- Step D: Key condition: n*arctan(1/gamma) + arctan(sinh(n*log(r))) < 2*pi ---")
print("  Tightest at gamma = gamma_1. Verify for each n = 1..85:")
print()

all_pass = True
max_lhs = mpf('0')
max_n = 0
min_gap = mpf('100')

print(f"  {'n':>4}  {'y = n*arctan':>14}  {'arctan(sinh(x))':>16}  {'LHS':>14}  {'2*pi - LHS':>14}  {'OK?':>6}")
print(f"  {'-'*75}")

for n in range(1, 86):  # Only n=1..85
    gamma = GAMMA_1
    t = 1/gamma
    x = n * log(sqrt(1 + t**2))   # n * log(r)
    y = n * atan(t)                 # n * arctan(1/gamma)
    y_star = atan(sinh(x))          # arctan(sinh(x))

    lhs = y + y_star
    gap = 2*pi - lhs

    ok = gap > 0

    if lhs > max_lhs:
        max_lhs = lhs
        max_n = n
    if gap < min_gap:
        min_gap = gap

    if not ok:
        all_pass = False

    if n <= 5 or n >= 83:
        status = "OK" if ok else "FAIL"
        print(f"  {n:>4}  {float(y):>14.8f}  {float(y_star):>16.8f}  {float(lhs):>14.8f}  {float(gap):>14.8f}  {status:>6}")

print(f"  {'...':>4}")
print()
print(f"  Maximum LHS (over n=1..85) = {float(max_lhs):.10f} at n={max_n}")
print(f"  2*pi = {float(2*pi):.10f}")
print(f"  Minimum gap = {float(min_gap):.10f}  > 0: {'YES' if min_gap > 0 else 'NO'}")
print(f"  All n=1..85 pass: {'YES' if all_pass else 'NO'}")

# =============================================================================
# Step E: At n=86, the condition fails
# =============================================================================
print("\n--- Step E: Threshold confirmation at n=86 ---")
n = 86
gamma = GAMMA_1
t = 1/gamma
x = n * log(sqrt(1 + t**2))
y = n * atan(t)
y_star = atan(sinh(x))
lhs = y + y_star
gap = 2*pi - lhs
cosh_cos = cosh(x) * cos(y)

print(f"  n=86, gamma=gamma_1:")
print(f"    y = {float(y):.10f}")
print(f"    y_* = arctan(sinh(x)) = {float(y_star):.10f}")
print(f"    LHS = y + y_* = {float(lhs):.10f}")
print(f"    2*pi = {float(2*pi):.10f}")
print(f"    gap = 2*pi - LHS = {float(gap):.10f}  {'> 0 OK' if gap > 0 else '< 0 FAIL!'}")
print(f"    cosh(x)*cos(y) = {float(cosh_cos):.10f}  {'< 1' if cosh_cos < 1 else '> 1 (F_n < 0!)'}")

# =============================================================================
# Step F: Verify cosh(x)*cos(y) < 1 for all n=1..85 numerically
# =============================================================================
print("\n--- Step F: Direct verification that cosh(x)*cos(y) < 1 at gamma=gamma_1 ---")
print()

max_val = mpf('-1')
max_n2 = 0
for n in range(1, 86):
    gamma = GAMMA_1
    t = 1/gamma
    x = n * log(sqrt(1 + t**2))
    y = n * atan(t)
    val = cosh(x) * cos(y)
    if val > max_val:
        max_val = val
        max_n2 = n

print(f"  Max cosh(x)*cos(y) over n=1..85 at gamma=gamma_1: {float(max_val):.12f}")
print(f"  Achieved at n={max_n2}")
print(f"  This is < 1: {'YES' if max_val < 1 else 'NO'}")
print(f"  F_n = 4*(1 - {float(max_val):.12f}) = {4*(1-float(max_val)):.12f} > 0")

# =============================================================================
# COMPLETE THEOREM STATEMENT
# =============================================================================
print()
print("=" * 70)
print("COMPLETE PROOF STATEMENT")
print("=" * 70)
print(f"""
  THEOREM 3 (RIGOROUS PAIR POSITIVITY):
  For n = 1, 2, ..., 85 and ALL gamma >= gamma_1 = 14.134725...:

    F_n(0, gamma) = 4*(1 - cosh(n*log(r)) * cos(n*arctan(1/gamma))) > 0,

  where r = sqrt(1 + 1/gamma^2).

  PROOF:
  Let x = n*log(r) and y = n*arctan(1/gamma). We show cosh(x)*cos(y) < 1.

  Step 1 (Boundary analysis): For x > 0, define y_* = arctan(sinh(x)).
    Then: cosh(x)*cos(y) = 1 iff y = y_* or y = 2*pi - y_* (in [0, 2*pi]).
    [Verification: cos(y_*) = sech(x), so cosh(x)*cos(y_*) = 1.]
    Consequently: cosh(x)*cos(y) < 1 iff y in (y_*, 2*pi - y_*).

  Step 2 (Lower bound: y > y_*):
    (a) By the log inequality: for t = 1/gamma in (0,1),
        f(t) = arctan(t) - (1/2)*log(1+t^2) satisfies f(0)=0,
        f'(t) = (1-t)/(1+t^2) > 0, so f(t) > 0. Hence log(r) < arctan(1/gamma).
        => x = n*log(r) < n*arctan(1/gamma) = y.
    (b) By the arctan-sinh inequality: g(x) = arctan(sinh(x)) - x satisfies
        g(0) = 0, g'(x) = sech(x) - 1 < 0, so g(x) < 0 for x > 0.
        => y_* = arctan(sinh(x)) < x < y.

  Step 3 (Upper bound: y < 2*pi - y_*):
    Equivalently: y + y_* < 2*pi, i.e., n*arctan(1/gamma) + arctan(sinh(n*log(r))) < 2*pi.

    (a) Monotonicity: both arctan(1/gamma) and log(r) are decreasing in gamma,
        so the LHS is maximized at gamma = gamma_1 (the smallest allowed value).

    (b) Numerical verification at gamma = gamma_1 for n = 1, ..., 85:
        max LHS (at n=85) = {float(max_lhs):.10f} < 2*pi = {float(2*pi):.10f}.
        Min gap = 2*pi - max LHS = {float(min_gap):.10f} > 0. ✓

        This is a rigorous finite computation at one point (gamma_1) for each n.
        Verified using mpmath with {mp.dps}-digit arithmetic.

  Steps 2 and 3 together give y in (y_*, 2*pi - y_*), hence cosh(x)*cos(y) < 1.
  Therefore F_n(0, gamma) = 4*(1 - cosh(x)*cos(y)) > 0.

  THRESHOLD: (computed in Step E above)
  At n=86, gamma=gamma_1: y + y_* = 6.2872439286 > 2*pi = {float(2*pi):.10f}.
  So the condition FAILS at n=86, and F_86(0, gamma_1) < 0 as verified. ✓

  QED
""")

# Save the key numbers
import json
result = {
    "theorem": "Pair positivity: F_n(0,gamma) > 0 for n=1..85, all gamma >= gamma_1",
    "proof_method": "Boundary analysis: cosh(x)*cos(y) < 1 iff y in (arctan(sinh(x)), 2pi - arctan(sinh(x)))",
    "key_inequality": "log(r) < arctan(1/gamma) for gamma > 1 [f'(t) = (1-t)/(1+t^2) > 0 for t in (0,1)]",
    "arctan_sinh_inequality": "arctan(sinh(x)) < x for x > 0 [g'(x) = sech(x)-1 < 0]",
    "max_LHS_at_gamma1": float(max_lhs),
    "max_LHS_n": max_n,
    "two_pi": float(2*pi),
    "gap": float(2*pi - max_lhs),
    "n_threshold": 85,
    "n_fails": 86,
    "F86_gamma1": float(cosh(86 * log(sqrt(1 + (1/GAMMA_1)**2))) * cos(86 * atan(1/GAMMA_1))),
    "mpmath_precision_digits": mp.dps,
}
with open("/Users/scones/projects/riemann/results/clean_proof_verified.json", "w") as f:
    json.dump(result, f, indent=2)
print("  Results saved to results/clean_proof_verified.json")
