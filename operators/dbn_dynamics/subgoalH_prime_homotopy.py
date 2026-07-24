#!/usr/bin/env python3
"""
subgoalH_prime_homotopy.py — LEDGER §154 (Lens 2, our own object): the PRIME HOMOTOPY. Dial the prime contribution to
the Weil screw function from t=0 (pure archimedean) to t=1 (true zeta) and BEYOND (t>1, over-weighted primes), and watch
the compact-G_a low spectrum morph. The arch-only operator has a LINEAR spectrum (λ_k~k, floor ~0.05, strongly positive,
§143); the primes restructure it into the k² near-null ladder (floor ~1e-6). This tool watches the COLLAPSE MECHANISM:
  g_t(x) = g_arch(x) + t · g_prime(x),   g_prime = screw_g − screw_g_no_primes  (the pure prime sum).
QUESTIONS (each branch = a distinct bankable fact):
  (1) TRAJECTORY: does λ_min(t) glide monotonically to the ~1e-6 floor, or DIVE and RECOVER? (recovery = a restoring
      mechanism = RH content made visible.)
  (2) EXPONENT FLIP: where in t does λ_k~k^p flip from p≈1 (arch) to p≈2 (full)? sharp t* (phase transition) or smooth?
  (3) OVERSHOOT (the sharp one): at t>1 does λ_min go NEGATIVE? If it crosses 0 just past t=1, the TRUE primes sit
      exactly at the positivity boundary (t=1 critical) — connects to the §122 arch↔prime balance & the explicit formula.
      If λ_min stays >0 to t=1.5, there is SLACK (margin) — which is what a proof needs.
Compact G_a/N formulation (physical for arch AND full; avoids the §42 sine-Galerkin arch artifact).
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

a=2.0; M=1600
Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
dif=np.arange(-M,M+1)*h
print(f"precomputing arch & prime screw values on the grid (a={a}, M={M})...", flush=True)
g_arch = np.array([screw_g_no_primes(d,Lam,Nmax) for d in dif])
g_full = np.array([screw_g(d,Lam,Nmax) for d in dif])
g_prime = g_full - g_arch
idx=np.arange(M+1)
X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
sw=np.sqrt(w); Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c)
Np=P@Nm@P; Np=(Np+Np.T)/2
ev_N,V_N=np.linalg.eigh(Np); keep=ev_N>1e-12*ev_N.max()
S=V_N[:,keep]/np.sqrt(ev_N[keep])

def low_spectrum(t, nk=12):
    gv = g_arch + t*g_prime
    G=gv[(idx[:,None]-idx[None,:])+M]
    Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Gp=P@Gm@P
    Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[:nk]

def expo(ev):
    ev=ev[ev>0]; k=np.arange(1,len(ev)+1); kk=k[:8]; ee=ev[:8]
    if len(ee)<3: return np.nan
    return np.polyfit(np.log(kk),np.log(ee),1)[0]

print("PRIME HOMOTOPY  g_t = g_arch + t·g_prime   (t=0 arch → t=1 true ζ → t>1 over-weighted)\n")
print(f"{'t':>5} {'λ_min':>12} {'λ_2':>11} {'λ_3':>11} {'λ_5':>11} {'exponent p':>10}  {'sign':>5}")
ts=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.0,1.05,1.1,1.2,1.3,1.4,1.5]
traj=[]
for t in ts:
    ev=low_spectrum(t); p=expo(ev)
    traj.append((t,ev[0]))
    sign = 'POS' if ev[0]>0 else '**NEG**'
    print(f"{t:>5.2f} {ev[0]:>12.4e} {ev[1]:>11.4e} {ev[2]:>11.4e} {ev[4]:>11.4e} {p:>10.3f}  {sign:>5}",flush=True)
# trajectory analysis
tt=np.array([r[0] for r in traj]); lm=np.array([r[1] for r in traj])
i1=np.argmin(np.abs(tt-1.0))
print(f"\n λ_min at t=1 (true ζ): {lm[i1]:.4e}")
print(f" λ_min MONOTONE decreasing on [0,1]? {'YES' if np.all(np.diff(lm[tt<=1.0])<1e-10) else 'NO — has a dip/recovery'}")
neg=tt[lm<0]
if len(neg): print(f" λ_min goes NEGATIVE first at t={neg[0]:.2f}  ⟹ criticality near there")
else: print(f" λ_min stays POSITIVE through t={tt[-1]:.1f}  ⟹ SLACK/margin beyond the true primes (min λ_min={lm.min():.3e})")
print("\nREAD: monotone-glide vs dip-recovery (restoring force); smooth vs sharp exponent flip (phase transition);")
print("NEGATIVE just past t=1 (true primes critical) vs SLACK to t=1.5 (margin). Each branch = a distinct structural fact.")
