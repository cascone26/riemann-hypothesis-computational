#!/usr/bin/env python3
"""§308: analyze the EXPLICIT Polya kernel Phi(u) (Xi(t)=int Phi cos(tu)du) against classical real-zeros
sufficient conditions. Polya/de Bruijn: sufficient conditions for int Phi(u)cos(tu)du to have ONLY REAL zeros
include Phi>0 even + certain log-concavity/monotonicity. Riemann's Phi must NARROWLY MISS the simplest ones
(else RH proven). Locate which it meets vs misses: log-concavity of Phi, monotonicity of Phi'/Phi, and the
sign of the 'de Bruijn' quantity. Charisma: REPORT (structural diagnostics on the explicit kernel)."""
import mpmath as mp
mp.mp.dps=30
def Phi(u):
    s=mp.mpf(0); e2u=mp.e**(2*u)
    for n in range(1,80):
        term=(2*mp.pi**2*n**4*mp.e**(mp.mpf(9)*u/2)-3*mp.pi*n*n*mp.e**(mp.mpf(5)*u/2))*mp.e**(-mp.pi*n*n*e2u)
        s+=term
        if abs(term)<mp.mpf(10)**(-35) and n>3: break
    return s
# Phi is even in u (Phi(-u)=Phi(u)); analyze on u>=0.
grid=[mp.mpf(x)/10 for x in range(0,46,2)]  # u=0,0.2,...,4.4
print("u,  Phi,  Phi'/Phi,  (logPhi)''  [logPhi concave <=> (logPhi)'' < 0]:",flush=True)
lc_all=True; peak=None
prev=None
for u in grid:
    p=Phi(u)
    d1=mp.diff(Phi,u); d2=mp.diff(Phi,u,2)
    dlog=d1/p
    d2log=(d2*p-d1*d1)/p**2
    lc = d2log<0
    lc_all = lc_all and lc
    print(f"  u={mp.nstr(u,3):>5}: Phi={mp.nstr(p,7):>12}  Phi'/Phi={mp.nstr(dlog,6):>11}  (logPhi)''={mp.nstr(d2log,6):>11}  {'concave' if lc else 'CONVEX!'}",flush=True)
print(f"  => log Phi concave on grid: {lc_all}",flush=True)
print("",flush=True)
# de Bruijn / Newman: the relevant object is whether Phi is in the class guaranteeing real zeros.
# Newman: writing Phi ~ exp(-Q(u)), real-zeros-sufficiency ties to Q being of a growth class. Riemann Phi ~ exp(-pi e^{2u})
# grows DOUBLE-EXPONENTIALLY -> outside Polya's simplest sufficient class (which needs slower growth). Show the growth:
print("Phi decay (double-exponential ~ exp(-pi e^{2u}) => OUTSIDE Polya's simplest real-zeros class):",flush=True)
for u in [mp.mpf(1),mp.mpf(2),mp.mpf(3)]:
    p=Phi(u); pred=mp.e**(-mp.pi*mp.e**(2*u))
    print(f"  u={u}: Phi={mp.nstr(p,4)}  ~exp(-pi e^2u)={mp.nstr(pred,4)}  log(-logPhi)={mp.nstr(mp.log(-mp.log(p)),5)} (linear in u => double-exp)",flush=True)
print("DONE",flush=True)
