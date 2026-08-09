#!/usr/bin/env python3
"""§928 (autonomous, new): is the RH-corner razor UNIVERSAL (same for Dirichlet L) or ζ-specific? Build the Weil form
for L(s,χ) χ mod 4 (odd, conductor q=4): arch_χ(r)=Re ψ₀(3/4+ir/2)−ln π + ln 4 (odd shift + conductor floor-raise);
prime_χ=(1/π)Σ_{n odd,log n<4W} Λ(n)χ(n)/√n cos(r log n) with χ(n)=+1 if n≡1(4), −1 if n≡3(4); NO pole (non-principal).
Compare the corner (#neg, margin) to ζ. If same razor structure → universal (proof can't be ζ-arithmetic-specific)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vm(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=3200; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
CHI4=np.where(NS%2==0,0,np.where(NS%4==1,1,-1))  # χ mod 4 on n
def build(a,b,fourW,kind,frac=0.7,dr=0.004):
    W=fourW/4; L=b-a; Ncap=L*2*W/np.pi; N=max(4,int(round(frac*Ncap)))
    RG=np.arange(-(b+60),(b+60),dr); CEN=np.linspace(a,b,N)
    if kind=='zeta': OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi)
    else:            OM=np.real(digamma(0.75+0.5j*RG))-np.log(np.pi)+np.log(4.0)  # χ4 odd + conductor
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    m=LOGN<fourW; ln=LOGN[m]
    if kind=='zeta': lw=LAMW[m]
    else:            lw=(LAMW*CHI4)[m]   # χ-twisted prime weights
    cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N):
        prime[j,:]=2.0*(((P[j][None,:]*P)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    if kind=='zeta':
        def fc(c):
            y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
        fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    else: pole=np.zeros((N,N))  # no pole for non-principal L
    gw,gv=eigh(Gm); B=gv@np.diag(1/np.sqrt(gw))@gv.T
    ev=eigh(B@(arch-prime)@B,eigvals_only=True)
    margin=float(eigh(arch+pole-prime,Gm,eigvals_only=True)[0])
    floor=float(eigh(arch+pole,Gm,eigvals_only=True)[0])
    return N,int((ev<-1e-11).sum()),float(ev[0]),margin,floor
print("="*90)
print("§928 UNIVERSALITY: L(s,χ4) corner vs ζ corner. #neg(M_σ), min M_σ, full margin, floor")
print("="*90)
print(f"{'window':>12} {'4W':>4} | {'ζ: #neg minMσ margin floor':>34} | {'χ4: #neg minMσ margin floor':>34}")
for (a,b) in [(14,44),(20,50),(30,60)]:
    for fourW in [3.0,4.0]:
        z=build(a,b,fourW,'zeta'); c=build(a,b,fourW,'chi4')
        zs=f"{z[1]} {z[2]:+.2e} {z[3]:+.2e} {z[4]:.3f}"
        cs=f"{c[1]} {c[2]:+.2e} {c[3]:+.2e} {c[4]:.3f}"
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} | {zs:>34} | {cs:>34}")
print("="*90)
print("If χ4 ALSO has #neg=1 razor (margin~1e-6+) → razor is UNIVERSAL (Selberg-class), NOT ζ-arithmetic-specific.")
print("If χ4 is comfortably positive (no corner) → the conductor floor-raise removes it; corner is ζ-specific (low conductor).")
print("="*90)

print("\n"+"="*90)
print("§928b FAIR TEST: χ4 at its OWN bottom-edge windows (χ4 first zeros ~6.0,10.2,13.0,16.3) vs interior")
print("="*90)
print(f"{'window':>12} {'4W':>4} {'χ4: #neg':>9} {'min Mσ':>11} {'margin':>11} {'floor':>8}")
for (a,b) in [(4,18),(4,24),(6,20),(2,16),(8,26)]:
    for fourW in [3.0,4.0]:
        c=build(a,b,fourW,'chi4')
        print(f"{f'[{a},{b}]':>12} {fourW:>4.1f} {c[1]:>9} {c[2]:>+11.2e} {c[3]:>+11.2e} {c[4]:>8.3f}")
print("If χ4 shows #neg=1 razor at ITS bottom edge → razor is UNIVERSAL (each L at its own spectrum edge).")
print("="*90)
