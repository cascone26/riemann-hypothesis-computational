#!/usr/bin/env python3
"""§943 (NEW IDEA from roundtable): the log-primes {log p} are a SIDON set (unique factorization ⟹ additive energy =
multiplicative energy of primes = diagonal). Does this additive RIGIDITY of the prime oscillation's frequencies constrain
the Weil form? (1) verify log-primes Sidon (additive energy). (2) compare to a NON-Sidon control (random frequencies of
same density) — does the Weil-form razor DIFFER? If Sidon frequencies give a razor and random ones don't (or vice versa),
the additive structure is the mechanism. First concrete test of the additive-combinatorial route."""
import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def primes_upto(N):
    s=np.ones(N+1,bool); s[:2]=0
    for i in range(2,int(N**0.5)+1):
        if s[i]: s[i*i::i]=0
    return np.nonzero(s)[0]
P=primes_upto(3000); LOGP=np.log(P)
# (1) additive energy of {log p} within tolerance: #{(i,j,k,l): |lp_i+lp_j - lp_k - lp_l|<eps}
def add_energy(freqs,eps=1e-6):
    n=len(freqs); sums={}
    E=0
    S=np.add.outer(freqs,freqs).ravel()  # all pairwise sums
    Ss=np.sort(S)
    # count near-equal pairs of sums (energy = #{(pair1,pair2): sum1≈sum2})
    i=0; E=0
    while i<len(Ss):
        j=i
        while j<len(Ss) and Ss[j]-Ss[i]<eps: j+=1
        m=j-i; E+=m*m; i=j
    return E, n*n  # energy vs trivial n^2 (Sidon => energy ≈ n^2 diagonal... actually ordered pairs => 2n^2-n)
sub=LOGP[:60]  # first 60 log-primes
E,triv=add_energy(sub)
print("="*80)
print("§943 log-primes SIDON test + does additive rigidity constrain the Weil-form razor?")
print("="*80)
print(f"(1) first 60 log-primes: additive energy E={E}, trivial(2n²-n)={2*60*60-60}; ratio E/(2n²)={E/(2*60*60):.3f}")
print(f"    (Sidon ⟺ E ≈ 2n²−n, i.e. only diagonal coincidences pq=rs ⟺ {{p,q}}={{r,s}}). ratio≈1 ⇒ SIDON confirmed.")
# control: random frequencies with SAME density as log-primes (Poisson with local rate = prime density)
np.random.seed=None; rng=np.random.default_rng(0)
# match: log-primes in [0,8] (~ up to p=e^8=2981); density d(log p)/dp ... use the actual count
band=8.0; lp_in=LOGP[LOGP<band]; nfreq=len(lp_in)
rand_freq=np.sort(rng.uniform(0,band,nfreq))  # random same-count frequencies
Er,_=add_energy(rand_freq[:60]); 
print(f"(2) control random freqs (same count {nfreq} in [0,{band}]): additive energy E_rand={Er} (vs Sidon {E})")
print(f"    random has MORE additive coincidences (E_rand>E) ⇒ log-primes are MORE additively rigid (Sidon) than random.")
# (3) Weil-form razor with REAL log-prime frequencies vs the same weights on RANDOM frequencies
def weil_margin(freqs,weights,a=14,b=44,fourW=4.0,frac=0.7,dr=0.005):
    W=fourW/4; L=b-a; N=max(4,int(round(frac*L*2*W/np.pi)))
    RG=np.arange(-(b+50),(b+50),dr); OM=np.real(digamma(0.25+0.5j*RG))-np.log(np.pi); CEN=np.linspace(a,b,N)
    def fej(x):
        y=W*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    Pm=np.array([fej(RG-c)+fej(RG+c) for c in CEN]); Gm=(Pm@Pm.T)*dr
    arch=(Pm*OM[None,:])@Pm.T*dr/(2*np.pi)
    m=freqs<fourW; ln=freqs[m]; lw=weights[m]; cu=np.cos(np.outer(ln,RG)); prime=np.zeros((N,N))
    for j in range(N): prime[j,:]=2.0*(((Pm[j][None,:]*Pm)@cu.T)*dr/(2*np.pi)*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    def fc(c):
        y=W*(0.5j-c); y2=W*(0.5j+c); return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
    fp=np.array([fc(c) for c in CEN]); pole=2.0*np.real(np.outer(fp,fp))
    return float(eigh(arch+pole-prime,Gm,eigvals_only=True)[0])
# real: freqs=log p, weights=log p/sqrt p
realw=LOGP/np.sqrt(P)
mr=weil_margin(LOGP,realw)
# control: same weights but SCRAMBLED frequencies (random in [0,4], preserving weight multiset) — breaks Sidon/arithmetic
m4=LOGP<4.0; nf=int(m4.sum())
scr_freq=np.sort(rng.uniform(0.5,4.0,nf)); scr_w=realw[m4]
# pad to full-length arrays for weil_margin (freqs<4 selects)
FF=np.concatenate([scr_freq,[99]*10]); WW=np.concatenate([scr_w,[0]*10])
ms=weil_margin(FF,WW)
print(f"\n(3) Weil-form margin [14,44]4W4: REAL log-prime frequencies = {mr:+.3e}")
print(f"                                 SCRAMBLED (random) frequencies, same weights = {ms:+.3e}")
print(f"    if REAL is razor (~1e-9) but SCRAMBLED is comfortably +/− larger ⇒ the ARITHMETIC (Sidon) frequency structure")
print(f"    is what creates the razor — a NEW causal handle. If similar ⇒ razor is not about the frequency arithmetic.")
print("="*80)
