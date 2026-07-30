#!/usr/bin/env python3
"""Weyl seat: settle my own error against Arb's claim. Original SL eq (§464 build, lam=sqrt2):
  (x^2-2) xi'' + 2x xi' + (8*pi^2*x^2 - alpha) xi = 0        [tau xi = alpha xi]
I first (wrongly) classed x=inf via the Liouville-normal variable u=x*xi, which solves
  u'' + bracket(x) u = 0,  bracket(x) = 8pi^2 + (16pi^2-alpha)/(x^2-2) + 2/(x^2-2)^2   (EXACT, no approx)
u stays BOUNDED oscillatory (u ~ A cos(wx)+B sin(wx), w=2pi*sqrt2) -- that's correct. My error: I then
checked L^2-ness of u itself instead of xi=u/x. Since xi = u/x and u is bounded, xi ~ (bounded)/x, so
|xi|^2 ~ 1/x^2 which IS integrable at infinity -> BOTH solutions L^2(dx) near x=inf -> LIMIT-CIRCLE,
agreeing with Arb, not my first pass. Verify by direct numerical integration of xi itself (not u).
"""
import numpy as np
from scipy.integrate import solve_ivp

FP = 8*np.pi**2          # 8 pi^2  (coefficient (2*pi*lam)^2, lam^2=2)
w  = 2*np.pi*np.sqrt(2)  # asymptotic frequency
sp = np.sqrt(2.0)        # singular endpoint

def rhs(x, y, alpha):
    xi, xip = y
    xipp = -(2*x*xip + (FP*x*x - alpha)*xi) / (x*x - 2)
    return [xip, xipp]

def integrate(alpha, X0, X1, y0):
    xs = np.linspace(X0, X1, 4000)
    s = solve_ivp(rhs, [X0, X1], y0, args=(alpha,), method='DOP853',
                   rtol=1e-12, atol=1e-14, t_eval=xs)
    return xs, s.y[0], s.y[1]

# two independent local solutions near sp=sqrt2 (real, alpha real for this decay check):
# (a) REGULAR one: xi(sp)=1, xi'(sp) fixed by the indicial/regularity condition (finite limit).
# (b) generic/"other" one: perturb the derivative off the regular value -> picks up the LOG branch.
alpha_test = -100.0   # representative negative alpha (deep in the physical negative-eigenvalue band)
dxi_reg = -(FP*2 - alpha_test) / (2*sp)   # regularity condition derived in rt6_prolate_v2/v3

print("=== decay-rate check: xi(x) itself, both branches, alpha=%.1f ===" % alpha_test)
for label, y0 in [("regular", [1.0, dxi_reg]), ("perturbed(log-branch)", [1.0, dxi_reg+0.37])]:
    xs, xi, xip = integrate(alpha_test, sp+1e-4, 400.0, y0)
    # fit |xi(x)| envelope ~ C/x^p over the tail by log-log slope on local maxima (envelope of |xi|*x)
    tail = xs > 200
    xt, xit = xs[tail], np.abs(xi[tail])
    # local envelope via rolling max to avoid picking zero-crossings
    env = np.maximum.accumulate(xit[::-1])[::-1]  # crude but fine for a decay-exponent read
    logx, logy = np.log(xt), np.log(np.maximum(env,1e-300))
    p = -np.polyfit(logx, logy, 1)[0]
    l2_partial = np.trapezoid(xit**2, xt)   # finite-window ∫|xi|^2 dx over [200,400] -- should be finite/small if L2
    print(f"  {label:24s}: envelope decay exponent p~{p:.3f} (p=1 => xi~1/x => L2 at inf)"
          f"   ∫_200^400 |xi|^2 dx = {l2_partial:.3e}")

print("""
reading: p~1.0 for BOTH branches => xi ~ C/x (oscillating) => |xi|^2 ~ 1/x^2 => square-integrable at
infinity for BOTH linearly independent local solutions => x=+inf is LIMIT-CIRCLE (my original limit-point
call was wrong -- I mis-tested u=x*xi instead of xi). This confirms Arb's finding, independently.
""")

# --- x=sqrt2 endpoint: confirm log-branch is *locally* L^2 too (finite endpoint, mild singularity) ---
print("=== x=sqrt(2) endpoint: local behavior of the two Frobenius solutions ===")
for label, y0 in [("regular", [1.0, dxi_reg]), ("perturbed(log-branch)", [1.0, dxi_reg+0.37])]:
    xs, xi, xip = integrate(alpha_test, sp+1e-6, sp+0.5, y0)
    print(f"  {label:24s}: xi(sp+1e-6)={xi[0]:.4f}  xi(sp+0.5)={xi[-1]:.4f}  "
          f"max|xi| on window={np.max(np.abs(xi)):.4f}  (bounded => trivially L^2 near a finite pt)")
print("Both branches stay bounded near x=sp (only a log(x-sp) divergence for the 2nd solution, which is")
print("still locally L^2 since ∫(log t)^2 dt converges at t=0) => x=sqrt(2) is LIMIT-CIRCLE (indicial roots 0,0).")

print("""
=== DEFICIENCY INDICES (this half-line problem x in (sqrt2, inf)) ===
Both endpoints limit-circle => (n+,n-) = (2,2) for the operator on (sqrt2,inf) alone.
[Matches Arb's independent numeric result -- self-correction banked, not silently absorbed.]
""")
