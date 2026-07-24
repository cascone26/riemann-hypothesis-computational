#!/usr/bin/env python3
"""§199.3a: CLOSED-FORM entry pieces. Q[k,l]=∫_0^{2a} g(t)D_kl(t)dt, D_kl=C_kl+C_lk (cos-overlap, sympy exact).
g=boundary + arch_core + arch_tail + prime; arch_tail=-¼Φ0 + ¼e^{-t/2}Φ(e^{-2t}). Closed-form (sympy) ALL parts
EXCEPT the Lerch ¼e^{-t/2}Φ term (numeric mpmath remainder for now, → rigorous quadrature §199.3b). Validate the
assembled 4×4 vs mpmath weil_matrices."""
import sys, mpmath as mp
import sympy as sp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, prime_data
mp.mp.dps=40
aval=mp.mpf('0.549'); A=sp.Rational('549','1000'); t=sp.symbols('t',positive=True)
L2=sp.log(2); TWOA=2*A
Om0=sp.digamma(sp.Rational(1,4))-sp.log(sp.pi)     # -5.3722
Phi0=sp.zeta(2,sp.Rational(1,4))                    # 17.1973
wprime=sp.log(2)/sp.sqrt(2)                          # Λ(2)/√2
x=sp.symbols('x',real=True)
def C(k,l):  # C_kl(t)=∫_{t-a}^a d_k(x)d_l(x-t)dx  (t in [0,2a])
    dk=(k*sp.pi/(2*A))*sp.cos(k*sp.pi*(x+A)/(2*A))
    dl=(l*sp.pi/(2*A))*sp.cos(l*sp.pi*(x-t+A)/(2*A))
    return sp.integrate(dk*dl,(x,t-A,A))
modes=[1,3,5,7]
print("§199.3a closed-form 4x4 (Lerch part numeric) vs mpmath  @ a=0.549")
Qf,Qa,Mass,Vv,xs,w=weil_matrices(aval,24,1200)
ev=[c for c in range(24) if (c+1)%2==1][:4]
maxerr=mp.mpf(0)
for ik,k in enumerate(modes):
    for il,l in enumerate(modes):
        if il<ik: continue
        D=sp.simplify(C(k,l)+C(l,k))
        gcf=Om0*(-t/2) + (-4*(sp.exp(t/2)+sp.exp(-t/2)-2)) + (-sp.Rational(1,4)*Phi0)   # boundary+core+(-Φ0/4)
        # prime part on [log2,2a]
        Icf = sp.integrate(gcf*D,(t,0,TWOA)) + sp.integrate(wprime*(t-L2)*D,(t,L2,TWOA))
        val_cf=mp.mpf(str(sp.N(Icf,40)))
        # Lerch remainder ¼∫_0^{2a} e^{-t/2}Φ(e^{-2t})D(t)dt  (numeric)
        Dnum=sp.lambdify(t,D,'mpmath')
        def lerch_int(tt):
            z=mp.e**(-2*tt); Phi=mp.nsum(lambda n:z**n/(n+mp.mpf(1)/4)**2,[0,mp.inf])
            return mp.e**(-tt/2)*Phi*Dnum(tt)
        val_L=mp.mpf('0.25')*mp.quad(lerch_int,[0,float(mp.log(2)),float(2*aval)])
        val=val_cf+val_L
        ref=Qf[ev[ik],ev[il]]
        err=abs(val-ref); maxerr=max(maxerr,err)
        print(f"  Q[{k},{l}]: closed+Lerch={mp.nstr(val,12)}  mpmath={mp.nstr(ref,12)}  err={mp.nstr(err,3)}")
print(f"\nmax |err| vs mpmath = {mp.nstr(maxerr,3)}  (small ⟹ closed-form decomposition CORRECT; only Lerch needs rigor)")
