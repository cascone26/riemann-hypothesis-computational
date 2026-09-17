#!/usr/bin/env python3
"""
CYCLE 293 — corrected full arch−prime Weil form, anchored to the CERTIFIED convention (weil_dirichlet_gram_arb.py).

Cycle-292 bug diagnosis (VERIFIED): (a) degree-2 gamma Γ(s)=Γ_ℝ(s)Γ_ℝ(s+1) (Legendre) needs the FULL
duplication density Re ψ(¼+ir/2)+Re ψ(¾+ir/2), NOT a single Re ψ(½+ir) (off by ~2× on the r-part);
(b) small-a basis is catastrophically ill-conditioned (cond~3e27, marginal) — use the certified large-a basis.
Certified convention (matched to weil_dirichlet_gram_arb.py, χ₄): Ω = log(cond) − deg·log π + Σ_j Re ψ(¼+μ_j/2+ir/2),
μ_j∈{0,1}; prime factor 2; digamma coeff 1. Cross-checked: Ω_ζ+Ω_χ4 = Ω_{ζ_K}, cond=|D|. deg-2 Epstein: Ω_Q=log(|D|/π²)+ψ(¼)+ψ(¾).

HARD GATE: h=1 must give n_neg=0 BOTH as direct degree-2 Epstein AND as ζ_K=ζ·L(χ₄)=(ζ-form)+(χ₄-form), and they must agree;
h=2,h=3 must give n_neg>=1. If not, the full-form route is reported fragile and the robust discriminator (off-line count + jump-weight sign) stands.
"""
import numpy as np, mpmath as mp
from math import log
from sympy import primerange
mp.mp.dps=30
SQ2PI=float(mp.sqrt(2*mp.pi)); PI=float(mp.pi)
def phihat(a,r): return SQ2PI/2*(mp.e**(-(r-a)**2/2)+mp.e**(-(r+a)**2/2))
def g(a,ap,y):
    s=a+ap; d=a-ap; sp=float(mp.sqrt(mp.pi))/2
    return mp.e**(-y*y/4)*sp*(mp.e**(-s*s/4)*mp.cos(d*y/2)+mp.e**(-d*d/4)*mp.cos(s*y/2))
def chi4(m): r=m%4; return 1 if r==1 else(-1 if r==3 else 0)
def psi(x,r): return float(mp.re(mp.digamma(mp.mpf(x)+1j*r/2)))
A=[10.,14.,18.,22.,26.]; n=len(A); NMAX=600

# ---- prime weight lists ----
def lam_zeta():
    L=[]
    for p in primerange(2,NMAX):
        pk=p; lp=log(p)
        while pk<NMAX: L.append((pk,lp)); pk*=p
    return L
def lam_chi4():
    L=[]
    for p in primerange(2,NMAX):
        pk=p; lp=log(p)
        while pk<NMAX:
            c=chi4(pk)
            if c: L.append((pk,c*lp))
            pk*=p
    return L
def reduced_principal_form(d):
    return (1,0,d//4) if d%4==0 else (1,1,(d+1)//4)
def lam_epstein(d):
    a2,b2,c2=reduced_principal_form(d); r=np.zeros(NMAX+1); B=int(2*np.sqrt(NMAX))+3
    for x in range(-B,B+1):
        for y in range(-B,B+1):
            v=a2*x*x+b2*x*y+c2*y*y
            if 1<=v<=NMAX: r[v]+=1
    a=r/r[1]; Lam=np.zeros(NMAX+1)
    for N in range(2,NMAX+1):
        s=a[N]*log(N)
        for e in range(1,N):
            if N%e==0: s-=Lam[e]*a[N//e]
        Lam[N]=s
    return [(N,Lam[N]) for N in range(2,NMAX+1) if abs(Lam[N])>1e-12]

def build(Omega, primes, pole):
    W=np.zeros((n,n))
    for j in range(n):
        for k in range(j,n):
            aj,ak=A[j],A[k]
            arch=float(mp.quad(lambda r: phihat(aj,r)*phihat(ak,r)*Omega(r),[-70,0,70]))/(2*PI)
            pr=2*sum(w/np.sqrt(N)*float(g(aj,ak,log(N))) for (N,w) in primes)
            pl=0.0
            if pole:
                pl=float(mp.re(phihat(aj,-0.5j)*phihat(ak,-0.5j)))+float(mp.re(phihat(aj,0.5j)*phihat(ak,0.5j)))
            W[j,k]=W[k,j]=pl+arch-pr
    return W

def rep(W): ev=np.linalg.eigvalsh(W); return int((ev<-1e-9).sum()), ev.min(), ev.max()/max(abs(ev.min()),1e-30)

# component forms (certified convention)
O_zeta = lambda r: log(1)  -1*log(PI) + psi('0.25',r)                       # cond 1, deg1, mu=0
O_chi4 = lambda r: log(4)  -1*log(PI) + psi('0.75',r)                       # cond 4, deg1, mu=1
Wz=build(O_zeta, lam_zeta(), pole=True)     # zeta has the pole
Wc=build(O_chi4, lam_chi4(), pole=False)
for nm,W in [("zeta (pole)",Wz),("L(chi4)",Wc),("zeta_K=zeta+chi4 [h=1 via product]",Wz+Wc)]:
    ng,mn,cd=rep(W); print(f"  {nm:38s} n_neg={ng} min-eig={mn:+.4e} cond={cd:.1e}")

print("\nDirect degree-2 Epstein (Ω_Q=log(|D|/π²)+ψ(¼)+ψ(¾)+pole), gate: h=1→0, h≥2→≥1")
known={4:0,20:8,23:3}
for d,h in [(4,1),(20,2),(23,3)]:
    O = (lambda dd: (lambda r: log(dd) -2*log(PI) + psi('0.25',r) + psi('0.75',r)))(d)
    W=build(O, lam_epstein(d), pole=True)
    ng,mn,cd=rep(W); print(f"  disc -{d} (h={h}): n_neg={ng} min-eig={mn:+.4e} cond={cd:.1e}  off-line={known[d]}  "
                           f"{'PASS' if (ng==0)==(known[d]==0) else 'MISMATCH'}")
