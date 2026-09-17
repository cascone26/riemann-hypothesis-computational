#!/usr/bin/env python3
"""§934: does the zero-free Weil form Q=arch+pole−prime have CLOSED-FORM SOS structure (like the function-field Gram VV*)?
At small comfortable configs, compute the LDL^T / Cholesky SOS in high precision and inspect the L-vectors for recognizable
arithmetic structure. Also test: does the SOS L (or its normalized form) GENERALIZE across configs (same pattern, diff
constants)? Structure that generalizes = a structural certificate. Opaque/generic = no closed-form (expected if RH-hard)."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh, cholesky
np.set_printoptions(precision=4,suppress=True,linewidth=140)
def build(a,b,fourW,N,dr=0.003):
    W=fourW/4; RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(P@P.T)*dr
    arch=(P*OM[None,:])@P.T*dr/(2*np.pi)
    def vm(M):
        lam=np.zeros(M+1); isc=np.zeros(M+1,dtype=bool)
        for p in range(2,M+1):
            if not isc[p]:
                isc[2*p::p]=True; lp=np.log(p); pk=p
                while pk<=M: lam[pk]=lp; pk*=p
        return lam
    NS=np.arange(2,3000); LOGN=np.log(NS); LW=vm(2999)[2:]/np.sqrt(NS); m=LOGN<fourW; ln=LOGN[m]; lw=LW[m]
    cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N): prime[j,:]=2.0*(((P[j][None,:]*P)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    return arch+pole-prime, Gm
print("="*88)
print("§934 closed-form SOS structure hunt: whitened Q=B(arch+pole−prime)B, Cholesky L, inspect for arithmetic structure")
print("="*88)
Ls=[]
for (a,b,fourW,N) in [(30,60,2.0,4),(40,70,2.0,4),(30,60,2.0,5)]:
    Q,Gm=build(a,b,fourW,N); gw,gv=eigh(Gm); Bm=gv@np.diag(1/np.sqrt(gw))@gv.T; Qw=Bm@Q@Bm
    ev=eigh(Qw,eigvals_only=True)
    try:
        Lc=cholesky(Qw,lower=True); Ls.append(Lc)
        print(f"\n[{a},{b}] 4W={fourW} N={N}: min eig={ev[0]:.4f} (PSD:{ev[0]>0}). Cholesky L (lower):")
        print(Lc)
        print(f"  L diagonal (the 'pivots' √): {np.diag(Lc)}")
        print(f"  ratios L[i+1,i]/L[i,i]: {[Lc[i+1,i]/Lc[i,i] for i in range(N-1)]}")
    except Exception as e:
        print(f"[{a},{b}] not PD: {e}")
print("="*88)
print("If L-entries show clean ratios / a recognizable pattern that GENERALIZES across configs → closed-form SOS structure.")
print("If opaque/config-specific → no closed-form zero-free SOS (expected; would-be closed form = RH certificate).")
print("="*88)
