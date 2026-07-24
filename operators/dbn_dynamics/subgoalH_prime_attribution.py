#!/usr/bin/env python3
"""
subgoalH_prime_attribution.py — LEDGER §156 (Lens 2, make-or-break): the prime homotopy (§154) showed arch-only λ_min=
-4.97 lifted to the knife-edge (~0) by the FULL prime sum at t=1. WHICH primes do the lifting? Add primes one at a time
(2,3,5,7,...) at full strength and watch λ_min climb. FEW-BODY (2,3,5 do most, tail small) ⟹ the mechanism is a finite,
analyzable perturbation (real leverage). COLLECTIVE (every prime matters comparably) ⟹ homotopy is a visualization, wall.
Per-prime screw contribution (incl. powers): g_p(d) = Σ_{k≥1, p^k≤e^{|d|}} (log p)·p^{-k/2}·(|d| − k·log p).
Compact G_a/N low spectrum; cumulative arch + {first m primes}.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g_no_primes, vonmangoldt_upto

a=2.0; M=1600
Nmax=int(np.floor(np.exp(2*a)))+1
x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
dif=np.arange(-M,M+1)*h; adif=np.abs(dif)
Lam=vonmangoldt_upto(min(Nmax,30_000_000))
g_arch=np.array([screw_g_no_primes(d,Lam,Nmax) for d in dif])

def primes_upto(n):
    sieve=np.ones(n+1,bool); sieve[:2]=False
    for i in range(2,int(n**.5)+1):
        if sieve[i]: sieve[i*i::i]=False
    return np.where(sieve)[0]

def g_prime_single(p):
    # contribution of prime p (all powers p^k ≤ e^{|d|}) on the grid
    out=np.zeros_like(dif); lp=np.log(p); k=1
    while p**k <= np.exp(adif.max()):
        mask = (p**k) <= np.exp(adif)
        out += mask * (lp * p**(-k/2.0) * (adif - k*lp))
        k+=1
    return out

# compact G_a/N whitening setup (fixed)
idx=np.arange(M+1); X=x[:,None]; Y=x[None,:]
N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
sw=np.sqrt(w); Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
c=sw/np.linalg.norm(sw); Pp=np.eye(M+1)-np.outer(c,c); Np=Pp@Nm@Pp; Np=(Np+Np.T)/2
ev_N,V_N=np.linalg.eigh(Np); keep=ev_N>1e-12*ev_N.max(); Sw=V_N[:,keep]/np.sqrt(ev_N[keep])
def lam_min(gv):
    G=gv[(idx[:,None]-idx[None,:])+M]; Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2
    Gp=Pp@Gm@Pp; Ar=Sw.T@((Gp+Gp.T)/2)@Sw
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[0]

ps=primes_upto(Nmax)
print(f"Prime attribution (a={a}, primes up to {Nmax}: {list(ps)})")
print(f"arch-only λ_min = {lam_min(g_arch):+.4e}\n")
print(f"{'+prime p':>9} {'λ_min':>12} {'lift Δ':>11} {'cum lift %':>10}   {'|gap to 0|':>11}")
gv=g_arch.copy(); prev=lam_min(g_arch); base=prev
lifts=[]
for p in ps:
    gv=gv+g_prime_single(p)
    lm=lam_min(gv); d=lm-prev; lifts.append((p,d))
    cumpct=100*(lm-base)/(0-base)
    print(f"{p:>9} {lm:>12.4e} {d:>+11.4e} {cumpct:>9.1f}% {abs(lm):>11.3e}",flush=True)
    prev=lm
lifts.sort(key=lambda z:-abs(z[1]))
tot=sum(abs(d) for _,d in lifts)
print(f"\ntop lifters (|Δλ_min| share): " + ", ".join(f"p={p}:{100*abs(d)/tot:.0f}%" for p,d in lifts[:6]))
print(f"first 3 primes (2,3,5) share of total lift: {100*sum(abs(d) for p,d in lifts if p in(2,3,5))/tot:.0f}%")
print("\nREAD: if 2,3,5 supply most of the lift (>~60%) & the tail decays fast ⟹ FEW-BODY = analyzable finite mechanism")
print("(real leverage). If lift is spread ~evenly / decays like 1/p ⟹ COLLECTIVE = the homotopy is a visualization, wall.")
