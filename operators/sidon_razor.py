#!/usr/bin/env python3
"""§944 (FIXED causal test of the roundtable idea): does the ARITHMETIC (Sidon/unique-factorization) frequency structure
of the primes CAUSE the razor? REAL Weil form (full Λ(n)/√n, prime powers → true razor ~1e-9) vs SCRAMBLED: same weight
multiset on RANDOM frequencies (breaks the multiplicative/Sidon structure). If REAL razor but SCRAMBLED comfortable ⇒ the
arithmetic IS the razor's cause (a NEW causal handle for the additive-combinatorial route). If both razor ⇒ not the cause."""
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
NMAX=3000; LAM=vm(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)
def margin(freqs,weights,a=14,b=44,fourW=4.0,frac=0.7,dr=0.004):
    W=fourW/4; L=b-a; N=max(4,int(round(frac*L*2*W/np.pi)))
    RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    Pm=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(Pm@Pm.T)*dr
    arch=(Pm*OM[None,:])@Pm.T*dr/(2*np.pi)
    m=(freqs>0)&(freqs<fourW); ln=freqs[m]; lw=weights[m]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N): prime[j,:]=2.0*(((Pm[j][None,:]*Pm)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    return float(eigh(arch+pole-prime,Gm,eigvals_only=True)[0])
# REAL: freqs = log n (prime powers), weights = Λ(n)/√n
mreal=margin(LOGN,LAMW)
# select the in-band ones (log n < 4) to preserve exact count & weight multiset
mb=(LOGN>0)&(LOGN<4.0); band_freq=LOGN[mb]; band_w=LAMW[mb]; nf=len(band_freq)
print("="*82)
print("§944 does the arithmetic frequency structure CAUSE the razor? REAL vs SCRAMBLED (10 seeds)")
print("="*82)
print(f"REAL Weil margin [14,44]4W4 (Λ(n)/√n at log n): {mreal:+.4e}   (in-band terms: {nf})")
print(f"\nSCRAMBLED: same {nf} weights on RANDOM frequencies in (0.4,4.0):")
rng=np.random.default_rng(1); scr=[]
for s in range(10):
    rf=np.sort(rng.uniform(0.4,3.99,nf))
    ms=margin(np.concatenate([rf,[99.0]]),np.concatenate([band_w,[0.0]]))
    scr.append(ms)
scr=np.array(scr)
print(f"  scrambled margins: {np.array2string(scr,precision=3)}")
print(f"  scrambled: mean={scr.mean():+.3e}, min={scr.min():+.3e}, max={scr.max():+.3e}, #razor(|·|<1e-4)={int((np.abs(scr)<1e-4).sum())}")
print("="*82)
print(f"REAL |margin|={abs(mreal):.1e}. If scrambled |margins| are ALL ≫ real (comfortable, not razor) ⇒ the arithmetic")
print(f"(unique-factorization/Sidon) frequency structure CAUSES the razor — a NEW causal handle. If scrambled also razor ⇒ no.")
print("="*82)

print("\n"+"="*82)
print("§944b DISENTANGLE: is it the EXACT arithmetic positions, the density, or the weight-pairing?")
print("="*82)
# D: perturb log-prime positions by δ (keep weights with their perturbed freq) — tests fragility to exact positions
print("D) perturb log-primes log p → log p + δ (keeps density+weights-at-positions, breaks EXACT unique-factorization):")
for delta in [0.001,0.01,0.05,0.15]:
    ms=[]
    for s in range(6):
        pert=band_freq+rng.normal(0,delta,nf)
        ms.append(margin(np.concatenate([pert,[99.0]]),np.concatenate([band_w,[0.0]])))
    ms=np.array(ms); print(f"   δ={delta:.3f}: margins {np.array2string(ms,precision=3)}  mean={ms.mean():+.2e}")
# E: keep EXACT log-prime positions, PERMUTE the weights among them — tests if the specific Λ(n)/√n-at-log n pairing matters
print("\nE) keep exact log-prime positions, PERMUTE weights among them (breaks weight↔frequency arithmetic pairing):")
ms=[]
for s in range(8):
    pw=rng.permutation(band_w)
    ms.append(margin(np.concatenate([band_freq,[99.0]]),np.concatenate([pw,[0.0]])))
ms=np.array(ms); print(f"   permuted-weight margins: {np.array2string(ms,precision=3)}  mean={ms.mean():+.2e}, #razor={int((np.abs(ms)<1e-4).sum())}")
print("="*82)
print("READ: D razor-stable under small δ but breaks at large δ ⇒ robust general arithmetic structure, not knife-edge exact.")
print("     E negative ⇒ the SPECIFIC weight Λ(n)/√n AT log n (the arithmetic pairing) is essential, not just the position set.")
print("="*82)
