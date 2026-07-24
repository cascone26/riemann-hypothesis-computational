# §38 DECISIVE follow-up (RESEARCH_QUEUE #1): is "Z = unique conical max of Weil margin" REAL or a
# freeze-A / FE-mismatch artifact? Restrict deformations to DENSITY-PRESERVING (sum theta_i = 0, i.e.
# rearrange primes with no net expansion => the archimedean/FE main term stays balanced, freeze-A valid).
# If Z is STILL a strict conical max within this submanifold => real variational characterization.
import numpy as np, scipy.linalg, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import prime_deformation_cartography as PD

N=10; idx=list(range(-N,N+1))
arch=PD.build_arch(idx); B=PD.build_B(idx)
BP=PD.BASE_PRIMES
def margin(primes): return PD.signed_lambda_min(arch,B,PD.comb_from_primes(primes),idx)[0]
def deform(theta): return [p*np.exp(t) for p,t in zip(BP,theta)]
m0=margin(BP)
print(f"baseline lambda_min(Z) = {m0:+.3e}   primes={BP}")

rng=np.random.RandomState(0)
c=np.ones(len(BP))/np.sqrt(len(BP))     # density constraint direction (sum theta_i=0 <=> orthogonal to c)
def proj_out(v): return v - (v@c)*c     # project onto density-preserving subspace
for label, constrained in [("UNCONSTRAINED (reproduces §38)", False), ("DENSITY-PRESERVING (sum theta=0)", True)]:
    for eps in [0.03, 0.08]:
        both_down=0; mixed=0; both_up=0; K=60
        worst=1e9
        for _ in range(K):
            v=rng.randn(len(BP))
            if constrained: v=proj_out(v)
            n=np.linalg.norm(v)
            if n<1e-9: continue
            v=v/n
            dp=margin(deform(+eps*v))-m0
            dm=margin(deform(-eps*v))-m0
            worst=min(worst, dp, dm)
            if dp<0 and dm<0: both_down+=1
            elif dp>0 and dm>0: both_up+=1
            else: mixed+=1
        print(f"  {label:34s} eps={eps}:  both-down(max)={both_down:2d}/{K}  mixed(saddle)={mixed:2d}  both-up(min)={both_up:2d}  | worst dlm={worst:+.2e}")
print("\nIf DENSITY-PRESERVING still shows both-down dominant => Z is a REAL conical max (not freeze-A artifact).")
print("If it flattens to mixed/both-up => §38's max was a density-mismatch artifact.")
