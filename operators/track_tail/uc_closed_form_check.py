#!/usr/bin/env python3
"""u_c closed form check (see tail_QFourier_loewner_note.md sec 6).
u(m) := (2m/A)*hhat(m), hhat(m) = int_0^2a H'(t)cos(mbt)dt computed EXACTLY via IBP on the
kernel's closed forms (no quadrature): hhat = (-1)^m H(2a) - H(0) + mb*int H sin(mbt),
with H(0)=0 exact (Lerch sum telescopes to Phi0/4).
Verified law (odd m): u(m) = -1 - (4wp/pi)*sin(m*b*log2) - 2/(pi^2 m) + 1.720/m + O(1/m^2)
  -> u_c = -1 EXACT, beta = 4wp/pi = 0.624052 EXACT; resid*m = +1.720 constant (4 digits)
     over m = 25..2001 at NL=100000 (m >= 5001 hits Lerch truncation resolution).
Derivation: b*2a=pi kills all sin(m pi) boundary terms; log-singularity (coef +1/2) gives
-(2/pi)*Si(m pi) -> -1 - 2/(pi^2 m) [odd m]; jump wp at log2 gives -(4wp/pi)sin(theta m);
smooth part gives the +c1/m with c1 = 1.720 measured (closed form = boundary values of the
smooth part's derivative, to be pinned in the rigorization)."""
import math, cmath
from scipy.special import digamma, zeta as hurwitz
A = 0.549; b = math.pi/(2*A); log2 = math.log(2); twoa = 2*A
wp = log2/math.sqrt(2); Om0 = digamma(0.25) - math.log(math.pi); Phi0 = hurwitz(2, 0.25)
NL = 100000
def Pp(p,t,s):
    if p==0: return 1/s
    if p==1: return t/s-1/s**2
    return t*t/s-2*t/s**2+2/s**3
def II(p,al,om,sc,t0,t1):
    s = complex(-al,om); F = lambda tt: cmath.exp(s*tt)*Pp(p,tt,s)
    v = F(t1)-F(t0)
    return v.real if sc=='cos' else v.imag
gt = [(-4,0,-0.5),(-4,0,0.5),(8,0,0.0),(-Om0/2,1,0.0),(-Phi0/4,0,0.0)]
for n in range(NL): gt.append((0.25/(n+0.25)**2,0,2*n+0.5))
H2a = sum(ag*(twoa**pg)*math.exp(-alp*twoa) for ag,pg,alp in gt) + wp*(twoa-log2)
def hhat(m):
    om = m*b
    s = sum(ag*II(pg,alp,om,'sin',0,twoa) for ag,pg,alp in gt)
    s += wp*(II(1,0,om,'sin',log2,twoa) - log2*II(0,0,om,'sin',log2,twoa))
    return ((-1)**m)*H2a + om*s
if __name__ == "__main__":
    theta = b*log2; beta = 4*wp/math.pi
    print("theta=%.6f beta=4wp/pi=%.6f" % (theta, beta))
    print("  m       u(m)        pred=-1-b*sin-2/(pi^2 m)   resid     resid*m")
    for m in [25, 51, 101, 251, 501, 1001, 2001]:
        u = (2*m/A)*hhat(m)
        pred = -1 - beta*math.sin(theta*m) - 2/(math.pi**2*m)
        print("%6d  %+.6f  %+.6f  %+.2e  %+.3f" % (m, u, pred, u-pred, (u-pred)*m))
    # prefactor identity spot check vs Isa's v2 / Loewner collapse is in the note (rel 3e-12)
