#!/usr/bin/env python3
"""§883: arb-CERTIFY Q_SS >= 0 (arch|_S dominates prime|_S) at the corner [14,30] 4W=4 N=6. Compute strong-prime
subspace basis V_S (float), build AP-prime in arb, project Q_SS = V_S^T (AP-prime) V_S, certify PSD (LDLT pivots>0).
Certifies Weil positivity on the critical subspace, unconditionally (zero-free)."""
import sys,os,numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from scipy.special import digamma
from scipy.linalg import eigh
from flint import arb,ctx
import arb_fejer_cert2 as ac, arb_zerofree_cert3 as zc
a,b,fourW,N,rc=14,30,4.0,6,100
W=fourW/4
# --- float: get strong-prime subspace basis V_S (in original G-metric space) ---
dr=0.004; RGRID=np.arange(-(b+100),(b+100),dr); OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
CEN=np.linspace(a,b,N)
def fej(x):
    y=W*x; out=np.ones_like(y); nz=np.abs(y)>1e-12; out[nz]=(np.sin(y[nz])/y[nz])**2; return out
Pm=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]); Gf=(Pm@Pm.T)*dr
archf=(Pm*OMEGA[None,:])@Pm.T*dr/(2*np.pi)
mask=np.log(np.arange(2,3000))<fourW; NS=np.arange(2,3000)
def vm(M):
    lam=np.zeros(M+1); isc=np.zeros(M+1,dtype=bool)
    for p in range(2,M+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=M: lam[pk]=lp; pk*=p
    return lam
LAM=vm(2999); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS); m2=LOGN<fourW; ln=LOGN[m2]; lw=LAMW[m2]
cosu=np.cos(np.outer(ln,RGRID)); primef=np.zeros((N,N))
for j in range(N):
    gp=((Pm[j][None,:]*Pm)@cosu.T)*dr/(2*np.pi); primef[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
primef=0.5*(primef+primef.T)
gw,gv=eigh(Gf); Bf=gv@np.diag(1/np.sqrt(gw))@gv.T
APw=Bf@archf@Bf; primew=Bf@primef@Bf
def fc(c):
    y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
fp=np.array([fc(c) for c in CEN]); polef=2.0*np.real(np.outer(fp,fp)); APw=APw+Bf@polef@Bf
floor=eigh(APw,eigvals_only=True)[0]; mu,Wv=eigh(primew); idx=np.argsort(mu)[::-1]; k=int((mu[idx]>floor).sum())
Ws=Wv[:,idx[:k]]; VS=Bf@Ws   # N x k float basis of strong-prime subspace in original space
print(f"§883 [{a},{b}] 4W={fourW} N={N}: k={k} strong-prime; float Q_SS min={eigh(VS.T@ (Bf@(archf+polef-primef)@Bf) @VS,eigvals_only=True)[0]:.3e}")
# --- arb: build AP-prime full NxN, project to Q_SS, certify PSD ---
ctx.prec=180; ac.W=arb(fourW)/4; ac.TW=2*ac.W; ac.c_tri=arb.pi()/ac.W
cs=[arb(float(c)) for c in CEN]
M=[[None]*N for _ in range(N)]; mr=0.0
for i in range(N):
    for j in range(i,N):
        M[i][j]=M[j][i]=ac.pole_entry(cs[i],cs[j])+zc.arch_entry_v3(cs[i],cs[j],arb(rc))-ac.prime_entry(cs[i],cs[j])
        mr=max(mr,float(M[i][j].rad()))
# Q_SS = VS^T M VS  (k x k arb). Note M is in ORIGINAL space (arch+pole-prime), and VS is in original space (with Bf).
# But APw used Bf-whitened; VS=Bf Ws already includes Bf, and M is un-whitened AP-prime. VS^T M VS is correct (the
# Weil form on the subspace spanned by VS columns, which are test-function coefficient vectors).
VSa=[[arb(float(VS[p,c])) for c in range(k)] for p in range(N)]
QSS=[[arb(0)]*k for _ in range(k)]
for i in range(k):
    for j in range(k):
        s=arb(0)
        for p in range(N):
            for q in range(N): s+=VSa[p][i]*M[p][q]*VSa[q][j]
        QSS[i][j]=s
# certify PSD: LDLT, all pivots > 0
d=[None]*k; L=[[arb(0)]*k for _ in range(k)]; psd=True
for i in range(k):
    sd=QSS[i][i]
    for kk in range(i): sd-=L[i][kk]*L[i][kk]*d[kk]
    if not (sd>0): psd=False; break
    d[i]=sd
    for j2 in range(i+1,k):
        v=QSS[j2][i]
        for kk in range(i): v-=L[j2][kk]*L[i][kk]*d[kk]
        L[j2][i]=v/d[i]
print(f"arb: maxrad(M)={mr:.1e}; Q_SS pivots: {[f'{float(x.mid()):.2e}(+-{float(x.rad()):.0e})' for x in d if x is not None]}")
print(f"Q_SS >= 0 CERTIFIED (all pivots > 0): {psd}  => arch|_S dominates prime|_S UNCONDITIONALLY on the critical subspace" if psd else "NOT certified (radii swallowed a pivot)")
