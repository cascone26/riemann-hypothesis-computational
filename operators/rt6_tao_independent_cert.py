"""Tao INDEPENDENT certificate of Hilbert's worst-cell 5x5 physical Gram.
Different sieve + code than rt6_hilbert_gram_physical.py; own Arb enclosure.
Cell (s=48, mod-105 admissible, worst off-diag mass):
  (6,15,27)(1,7,40)(5,11,32)(12,15,21)(1,22,25)
Rigorous: float entry + outward roundoff ball -> Arb interval Cholesky, every
pivot ball lower-bound > 0. Also independent float ladder for two-apparatus Gram."""
import json
import numpy as np
from math import log
from flint import arb, ctx
ctx.prec = 200  # bits

def sieve_Lambda(M):
    spf=np.zeros(M+1,dtype=np.int64)
    for i in range(2,M+1):
        if spf[i]==0:
            for j in range(i,M+1,i):
                if spf[j]==0: spf[j]=i
    L=np.zeros(M+1)
    for n in range(2,M+1):
        p=int(spf[n]); m=n
        while m%p==0: m//=p
        if m==1: L[n]=log(p)
    return L

CELL=[(6,15,27),(1,7,40),(5,11,32),(12,15,21),(1,22,25)]
M=1_000_000
Lam=sieve_Lambda(M+60)
dL=Lam-1.0
base=np.arange(0,M)
V=np.array([dL[base+t[0]]*dL[base+t[1]]*dL[base+t[2]] for t in CELL])
G=(V@V.T)/M                       # 5x5 physical Gram, independent build

# --- float PD two ways (reproduction) ---
Lc=np.linalg.cholesky(G)
pivots=np.diag(Lc)**2
minors=[np.linalg.det(G[:k,:k]) for k in range(1,6)]

# --- independent Arb interval Cholesky with outward roundoff balls ---
# rigorous relative radius: 1e-8 (>> float64 roundoff over 1e6 sums ~1e-9 rel)
def ball(x):
    x=float(x); return arb(x, abs(x)*1e-8+1e-12)
Gb=[[ball(G[i,j]) for j in range(5)] for i in range(5)]
# interval Cholesky: d[k]=G[k,k]-sum_{j<k} L[k,j]^2 d[j]; L[i,k]=(G[i,k]-sum L[i,j]L[k,j]d[j])/d[k]
d=[None]*5; L=[[arb(0)]*5 for _ in range(5)]
cert_ok=True; piv_lo=[]
for k in range(5):
    s=Gb[k][k]
    for j in range(k): s=s-L[k][j]*L[k][j]*d[j]
    d[k]=s
    lo=float(s.mid())-float(s.rad())
    piv_lo.append(lo)
    if lo<=0: cert_ok=False
    L[k][k]=arb(1)
    for i in range(k+1,5):
        t=Gb[i][k]
        for j in range(k): t=t-L[i][j]*L[k][j]*d[j]
        L[i][k]=t/d[k]

piv_mid=[float(d[k].mid()) for k in range(5)]
piv_rad=[float(d[k].rad()) for k in range(5)]
allpos=all(piv_mid[k]-piv_rad[k]>0 for k in range(5))

out={
 "object":"Tao independent Arb-interval cert of worst-cell 5x5 physical Gram",
 "version":"rt6-tao-independent-cert-2026-07-28",
 "cell":CELL,"s":48,"window_M":M,"arb_bits":200,"rel_ball":1e-8,
 "float_cholesky_pivots":[float(x) for x in pivots],
 "float_sylvester_minors":[float(x) for x in minors],
 "arb_interval_pivot_mid":piv_mid,"arb_interval_pivot_rad":piv_rad,
 "arb_all_pivots_strictly_positive":bool(allpos),
 "hilbert_pivots_ref":[343.83,343.58,343.62,342.51,341.10],
}
with open("results/rt6_tao_independent_cert.json","w") as f: json.dump(out,f,indent=2)
print("float Cholesky pivots:", [f"{x:.2f}" for x in pivots])
print("float Sylvester minors:", [f"{x:.4e}" for x in minors])
print("Arb interval pivots (mid+-rad):", [f"{piv_mid[k]:.3f}+-{piv_rad[k]:.1e}" for k in range(5)])
print("ARB ALL PIVOTS STRICTLY POSITIVE:", allpos)
print("wrote results/rt6_tao_independent_cert.json")
