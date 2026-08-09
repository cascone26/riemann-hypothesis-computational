#!/usr/bin/env python3
"""§886: the Szegő SYMBOL negative-set. σ(r)=Ω(r)/(2π) − P(r), Ω=Re ψ₀(¼+ir/2)−ln π, P(r)=(1/π)Σ_{log n<4W}
Λ(n)/√n cos(r log n). mean_r σ = floor (prime cos-terms average to 0). If σ≥0 everywhere → Weil positivity is
a trivial nonneg-symbol Toeplitz fact. If σ dips NEGATIVE on a bounded set → positivity is the nontrivial
band-limited/Szegő phenomenon, and the razor should track the negative mass vs floor. Measure min σ, neg-set
width, neg mass ∫|σ_-|, across heights (floor climbs, neg excursions bounded → corner is where they cross)."""
import numpy as np
from scipy.special import digamma
np.set_printoptions(precision=4,suppress=True)

def vonmangoldt(N):
    lam=np.zeros(N+1); isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True; lp=np.log(p); pk=p
            while pk<=N: lam[pk]=lp; pk*=p
    return lam
NMAX=200000; LAM=vonmangoldt(NMAX); NS=np.arange(2,NMAX+1); LOGN=np.log(NS); LAMW=LAM[2:]/np.sqrt(NS)

def symbol(a,b,fourW,dr=0.002):
    r=np.arange(a,b,dr)
    Omega=np.real(digamma(0.25+0.5j*r))-np.log(np.pi)
    floor_r=Omega/(2*np.pi)
    mask=LOGN<fourW; ln=LOGN[mask]; lw=LAMW[mask]
    P=(1/np.pi)*(lw[:,None]*np.cos(np.outer(ln,r))).sum(axis=0)
    sig=floor_r-P
    return r,sig,floor_r,P

print("="*100)
print("§886 Szegő symbol σ(r)=Ω/2π − P(r): is it pointwise negative at the corner? does razor track neg-mass?")
print("="*100)
print(f"{'window':>14} {'4W':>4} {'mean σ (=floor)':>15} {'min σ':>11} {'max σ':>10} {'neg-frac':>9} {'neg-mass∫|σ-|':>13}")
for (a,b,fourW) in [(14,44,4.0),(14,44,3.0),(14,44,5.0),(14,44,2.0),
                    (60,90,4.0),(120,150,4.0),(300,330,4.0),(1000,1030,4.0),(5000,5030,4.0)]:
    r,sig,fl,P=symbol(a,b,fourW)
    neg=sig<0; negfrac=neg.mean(); negmass=(-sig[neg]).sum()*(r[1]-r[0]) if neg.any() else 0.0
    print(f"{f'[{a},{b}]':>14} {fourW:>4.1f} {sig.mean():>15.5f} {sig.min():>+11.4f} {sig.max():>+10.4f} {negfrac:>9.3f} {negmass:>13.4e}")
print("="*100)
print("σ.min<0 → corner positivity is the NONTRIVIAL band-limited/Szegő fact (nonneg-symbol triviality fails).")
print("As height climbs, floor=mean σ grows (∝ ln t/2π); watch where min σ crosses 0 → the floor>pnorm robust onset.")
print("="*100)
