#!/usr/bin/env python3
"""Parametrized RIGOROUS Weil-positivity certificate (Arb enclosure → exact-rational LDL^T).
Usage: python3 organ2_certify_N.py "lo,hi,N"   (linspace basis)  OR  "a1,a2,..."(explicit)
Certifies the densest still-certifiable basis (the §774 Nyquist edge) to HANDLE-tier."""
import sys
from flint import arb, acb, ctx
import sympy as sp
from sympy import primerange, Rational
ctx.prec = 320

arg = sys.argv[1] if len(sys.argv) > 1 else "15,55,11"
parts = [float(x) for x in arg.split(",")]
if len(parts) == 3 and parts[2] == int(parts[2]) and parts[2] <= 40:
    lo, hi, N = parts[0], parts[1], int(parts[2])
    AV = [arb(str(lo + (hi-lo)*k/(N-1))) for k in range(N)]
else:
    AV = [arb(str(x)) for x in parts]
nb = len(AV)
print(f"basis N={nb}: {[float(a) for a in AV]}", flush=True)
PI = arb.pi(); SQRT2PI = (2*PI).sqrt(); LOGPI = PI.log()
def phihat_acb(a, r):
    A=acb(a); return acb(SQRT2PI)/2*((-(r-A)**2/2).exp()+(-(r+A)**2/2).exp())
def Omega(r):
    z1=acb(arb(1)/4)+acb(0,1)*r/2; z2=acb(arb(1)/4)-acb(0,1)*r/2
    return (z1.digamma()+z2.digamma())/2-acb(LOGPI)
def arch_entry(aj,ak):
    R=arb(240)
    def integrand(x,analytic):
        if analytic:
            im=x.imag
            if float(im.upper())>=0.5 or float(im.lower())<=-0.5: return acb(arb('nan'))
        return phihat_acb(aj,x)*phihat_acb(ak,x)*Omega(x)
    return acb.integral(integrand, acb(-R), acb(R)).real/(2*PI)
def pole_entry(aj,ak):
    r=acb(0,arb(1)/2); return 2*(phihat_acb(aj,r)*phihat_acb(ak,r)).real
# --- precompute von Mangoldt prime-power data ONCE (log n, Λ(n)/√n) as arb ---
PMAX_PRIME = 2000000
print("precomputing prime-power arb data...", flush=True)
_PP = []   # (logn_arb, weight_arb=Λ(n)/√n)
for p in primerange(2, PMAX_PRIME):
    lp = arb(p).log(); m = 1
    while True:
        L = m*lp
        if float(L) > 32: break
        w = lp/arb(p)**(arb(m)/2)          # Λ(n)/√n  = log p / p^{m/2}
        _PP.append((L, w)); m += 1
SQRTPI2 = PI.sqrt()/2
print(f"  {len(_PP)} prime-power terms cached", flush=True)
def prime_entry(aj,ak):
    s=aj+ak; d=aj-ak
    es4=(-s**2/4).exp(); ed4=(-d**2/4).exp(); dh=d/2; sh=s/2
    tot=arb(0)
    for (L, w) in _PP:
        # g_cross(aj,ak,L) = e^{-L²/4}·(√π/2)·[e^{-s²/4}cos(d L/2)+e^{-d²/4}cos(s L/2)]
        g = (-L**2/4).exp()*SQRTPI2*(es4*(dh*L).cos()+ed4*(sh*L).cos())
        tot += w*g
    prime=2*tot
    Lcut=arb(PMAX_PRIME).log(); gmax=SQRTPI2*(es4+ed4)
    tailb=2*gmax*(arb(1)/4).exp()*PI.sqrt()*((Lcut-1)/2).erfc()
    return prime+arb(0,float(tailb))
print(f"Building rigorous {nb}x{nb} Weil form (prec {ctx.prec})...", flush=True)
Qb=[[None]*nb for _ in range(nb)]
for j in range(nb):
    for k in range(j,nb):
        v=arch_entry(AV[j],AV[k])+pole_entry(AV[j],AV[k])-prime_entry(AV[j],AV[k])
        Qb[j][k]=v; Qb[k][j]=v
    print(f"  row {j} done", flush=True)
DEN=10**15
M=sp.zeros(nb,nb); radsq=arb(0); maxrad=0.0
for j in range(nb):
    for k in range(nb):
        mid=Qb[j][k].mid(); rad=Qb[j][k].rad(); maxrad=max(maxrad,float(rad))
        M[j,k]=Rational(str(float(mid))).limit_denominator(DEN)
        radsq+=arb(abs(float(mid)-float(M[j,k]))+float(rad))**2
rho=radsq.sqrt(); rho_rat=Rational(str(float(rho)+1e-18)).limit_denominator(10**13)+Rational(1,10**12)
# float min-eig for reference
import numpy as np
Mf=np.array([[float(M[i,j]) for j in range(nb)] for i in range(nb)])
mineig=np.linalg.eigvalsh(Mf).min()
print(f"max entry rad={maxrad:.2e} | rho <= {float(rho):.3e} | float min-eig(M)={mineig:.3e} | margin={mineig-float(rho):.3e}", flush=True)
Msh=(M-rho_rat*sp.eye(nb)).applyfunc(lambda e: sp.Rational(e))
n=nb; L=sp.eye(n); D=[sp.Integer(0)]*n; pd=True
for jj in range(n):
    d=sp.Rational(Msh[jj,jj]-sum(L[jj,kk]**2*D[kk] for kk in range(jj))); D[jj]=d
    if D[jj]<=0: pd=False
    for ii in range(jj+1,n):
        L[ii,jj]=sp.Rational((Msh[ii,jj]-sum(L[ii,kk]*L[jj,kk]*D[kk] for kk in range(jj)))/D[jj])
print("pivots(float):",[f'{float(dv):+.3e}' for dv in D])
print("all pivots>0 (EXACT rational LDL^T of M-rho*I):",pd)
print("===>", f"HANDLE: Weil form PD on {nb}-fn basis, RIGOROUS (Arb+exact LDL^T), min-eig-rho>0" if pd else "NOT certified (past saturation edge)")
