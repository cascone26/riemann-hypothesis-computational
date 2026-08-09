#!/usr/bin/env python3
"""§894: arb-CERTIFY the §889 structure at a genuine corner [14,34] 4W=4 N=8. Two rigorous claims:
 (A) the exact-symbol form M_σ = arch−prime has EXACTLY ONE negative eigenvalue (inertia: one negative LDLT pivot,
     radius not crossing 0) — the single-negative-eigenvalue structure, certified.
 (B) the FULL Weil form Q = arch+pole−prime is PSD (all LDLT pivots > 0) — Weil positivity at this corner, certified.
Together: the rank-≤2 pole rigorously lifts the single negative eigenvalue of arch−prime to make Q ⪰ 0."""
import sys,os,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from scipy.special import digamma
from scipy.linalg import eigh
from flint import arb,ctx
import arb_fejer_cert2 as ac, arb_zerofree_cert3 as zc
a,b,fourW,N,rc=14,34,4.0,8,300
W=fourW/4
CEN=np.linspace(a,b,N)
# float sanity
dr=0.004; RG=np.arange(-(b+60),(b+60),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi)
def fej(x):
    y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
Pm=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gf=(Pm@Pm.T)*dr
archf=(Pm*OM[None,:])@Pm.T*dr/(2*np.pi)
from numpy import log,arange,sqrt
def vm(M):
    lam=np.zeros(M+1); isc=np.zeros(M+1,dtype=bool)
    for p in range(2,M+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=log(p); pk=p
            while pk<=M: lam[pk]=lp; pk*=p
    return lam
LAM=vm(2999); NS=arange(2,3000); LN=log(NS); LW=LAM[2:]/sqrt(NS); m=LN<fourW; ln=LN[m]; lw=LW[m]
cu=np.cos(np.outer(ln,RG)); primef=np.zeros((N,N))
for j in range(N):
    gp=((Pm[j][None,:]*Pm)@cu.T)*dr/(2*np.pi); primef[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
primef=0.5*(primef+primef.T)
gw,gv=eigh(Gf); Bf=gv@np.diag(1/np.sqrt(gw))@gv.T
evMs=eigh(Bf@(archf-primef)@Bf,eigvals_only=True)
print(f"§894 float check [{a},{b}] 4W={fourW} N={N}: M_σ eigs={np.array2string(evMs,precision=3)}; #neg={int((evMs<0).sum())}; full margin={eigh(Bf@(archf+2.0*np.real(np.outer([ (np.sin(W*(0.5j-c))/(W*(0.5j-c)))**2+(np.sin(W*(0.5j+c))/(W*(0.5j+c)))**2 for c in CEN],[ (np.sin(W*(0.5j-c))/(W*(0.5j-c)))**2+(np.sin(W*(0.5j+c))/(W*(0.5j+c)))**2 for c in CEN]))-primef)@Bf,eigvals_only=True)[0]:.3e}")

# ---- arb ----
ctx.prec=260; ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
cs=[arb(float(c)) for c in CEN]
# build arch, prime, pole arb matrices
AR=[[None]*N for _ in range(N)]; PR=[[None]*N for _ in range(N)]; PO=[[None]*N for _ in range(N)]; mr=0.0
for i in range(N):
    for j in range(i,N):
        AR[i][j]=AR[j][i]=zc.arch_entry_v3(cs[i],cs[j],arb(rc))
        PR[i][j]=PR[j][i]=ac.prime_entry(cs[i],cs[j])
        PO[i][j]=PO[j][i]=ac.pole_entry(cs[i],cs[j])
        mr=max(mr,float(AR[i][j].rad()),float(PR[i][j].rad()),float(PO[i][j].rad()))
# whitening B (float) applied to arb matrices: Msig_arb = B (AR-PR) B ; Q_arb = B (AR+PO-PR) B
Ba=[[arb(float(Bf[i][j])) for j in range(N)] for i in range(N)]
def whiten(Mrows):
    T=[[arb(0)]*N for _ in range(N)]
    for i in range(N):
        for j in range(N):
            s=arb(0)
            for k in range(N):
                for l in range(N): s+=Ba[i][k]*Mrows[k][l]*Ba[l][j]
            T[i][j]=s
    return T
def ldlt_inertia(M):
    n=len(M); d=[None]*n; L=[[arb(0)]*n for _ in range(n)]; npos=0; nneg=0; indef=False
    for i in range(n):
        s=M[i][i]
        for k in range(i): s-=L[i][k]*L[i][k]*d[k]
        if s>0: npos+=1
        elif s<0: nneg+=1
        else: indef=True; break
        d[i]=s
        for jj in range(i+1,n):
            v=M[jj][i]
            for k in range(i): v-=L[jj][k]*L[i][k]*d[k]
            L[jj][i]=v/d[i]
    return npos,nneg,indef,d
MsA=whiten([[AR[i][j]-PR[i][j] for j in range(N)] for i in range(N)])
QA =whiten([[AR[i][j]+PO[i][j]-PR[i][j] for j in range(N)] for i in range(N)])
p1,n1,ind1,d1=ldlt_inertia(MsA)
p2,n2,ind2,d2=ldlt_inertia(QA)
print(f"arb maxrad(entries)={mr:.1e}, prec={ctx.prec}")
print(f"(A) M_σ=arch−prime inertia: +{p1} / −{n1} {'(INDEF pivot swallowed by radius)' if ind1 else 'CERTIFIED'}  → exactly one neg eig: {n1==1 and not ind1}")
print(f"(B) Q=arch+pole−prime inertia: +{p2} / −{n2} {'(pivot swallowed)' if ind2 else 'CERTIFIED'}  → Q ⪰ 0 (PSD): {n2==0 and not ind2}")
if n2==0 and not ind2 and n1==1 and not ind1:
    print("=> RIGOROUS: arch−prime has exactly ONE negative eigenvalue, and the rank-≤2 pole lifts it to make Q ⪰ 0. §889-891 certified.")
