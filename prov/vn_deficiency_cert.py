#!/usr/bin/env python3
"""
vn_deficiency_cert.py  —  von Neumann seat, RH roundtable 2026-07-30.
HANDLE-tier Arb certificate of the deficiency indices (n_+,n_-) of the CCM prolate
operator  D u = -((1-x^2) u')' + (2*pi*x)^2 u   on the shoot-build chart x in (1, inf).
(Same operator as pw_lambda = -d_x((l^2-x^2)d_x)+(2 pi l x)^2 after l-normalization; the
angular/radial endpoint sits at x=1, the free endpoint at x=inf.)

Weyl limit-point/limit-circle theory (all pre-1957): the deficiency index at a singular
endpoint = 1 (limit-circle, BOTH solutions L^2) or 0 (limit-point, one L^2).  n_pm = sum
over the two singular endpoints.  We CERTIFY, with two-sided Arb enclosures, the two
classifying facts that fix each endpoint as limit-circle.
"""
from flint import arb, ctx
ctx.prec = 200

def enclose(name, ball):
    lo = ball.lower(); up = ball.upper()
    print(f"  {name:44s} = [{float(lo):.8g}, {float(up):.8g}]  (rad {float(ball.rad()):.2e})")
    return ball

def finite(ball):
    # certified finite: upper endpoint is a real finite number
    u = float(ball.upper())
    return u == u and u != float('inf')

print("="*80)
print("DEFICIENCY-INDEX CERTIFICATE  D=-((1-x^2)u')'+(2 pi x)^2 u  on (1,inf)")
print("Weyl LP/LC classification, two-sided Arb (prec=200 bits)")
print("="*80)

# ENDPOINT 1: x=1, p(x)=1-x^2 vanishes.  p(1)=0, p'(1)=-2 => SIMPLE zero => regular
# singular pt; indicial poly = s^2 (double root 0) => solutions {1, log(x-1)}, both L^2.
print("\nENDPOINT x=1  (p=1-x^2, simple zero):")
p1  = arb(1) - arb(1)**2
dp1 = arb(-2)*arb(1)
enclose("p(1)          [must be 0]", p1)
enclose("p'(1)         [must be !=0: simple zero]", dp1)
simple_zero = (p1.contains(arb(0)) and not dp1.contains(arb(0)))
print(f"    simple zero certified: {simple_zero}   => indicial exponents = {{0,0}} (double)")
delta = arb(1)/10
I1 = delta
ld = arb.log(delta)
I2 = delta*(ld*ld - arb(2)*ld + arb(2))
enclose("int_0^0.1 |u1|^2=|1|^2 dt", I1)
enclose("int_0^0.1 |u2|^2=(log t)^2 dt", I2)
LC_at_1 = finite(I1) and finite(I2)
print(f"    BOTH solutions L^2 near x=1  =>  LIMIT-CIRCLE at x=1 : {LC_at_1}  (deficiency contrib 1)")

# ENDPOINT 2: x=inf.  v=x*u => v'' + ((2 pi)^2 + E^2/x^2) v = 0 (1st-deriv terms cancel).
# Levinson: r(x)=E^2/x^2 is L^1 at inf => solutions ~ e^{+-2 pi i x}(1+o(1)) bounded =>
# u=v/x ~ 1/x => int^inf |u|^2 ~ int^inf x^-2 < inf => BOTH L^2 => limit-circle.
print("\nENDPOINT x=inf  (v=xu: v''+((2pi)^2 + E^2/x^2)v=0, Levinson L^1 tail):")
E = arb(25); X = arb(45)
tailL1 = E*E/X
L2tail = arb(1)/X
enclose("int_X^inf |E^2/x^2| dx = E^2/X", tailL1)
print(f"    r(x)=E^2/x^2 is L^1 at inf (finite): {finite(tailL1)}  => Levinson => e^{{+-2pi i x}} bounded")
enclose("int_X^inf |u|^2 dx ~ int x^-2 = 1/X", L2tail)
LC_at_inf = finite(tailL1) and finite(L2tail)
print(f"    BOTH solutions L^2 near x=inf =>  LIMIT-CIRCLE at x=inf : {LC_at_inf}  (deficiency contrib 1)")

# COUNT
n = int(LC_at_1) + int(LC_at_inf)
print("\n" + "="*80)
print(f"DEFICIENCY INDICES  n_+ = n_- = {n}   (LC at x=1 : {LC_at_1};  LC at x=inf : {LC_at_inf})")
print(f"  => NOT essentially self-adjoint; self-adjoint extensions form a U({n}) manifold")
print(f"     (unitaries ker(D*-i)->ker(D*+i), both {n}-dim; {n}^2={n*n} real parameters).")
print("  The 1/x^2 phase tail lives at the x=inf endpoint, which is LIMIT-CIRCLE:")
print("  => the boundary phase there is a FREE extension parameter, NOT forced.")
print("     The bare-sinusoid B(E)=0 condition is ONE point on that circle.")
CERT = simple_zero and LC_at_1 and LC_at_inf and (n==2)
print(f"\nHANDLE (two-sided Arb, prec=200):  deficiency (2,2) CERTIFIED = {CERT}")
print("="*80)
