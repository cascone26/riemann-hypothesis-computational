import numpy as np
from scipy.integrate import solve_ivp
Z=np.load('operators/zeros500.npy'); gam=Z[:5]; X=120.0
# Build-constant scan on the beta=0 (Friedrichs) slice, per Grothendieck.
# Frequency scale s: operator freq w=s*2pi, so FP=w^2 in the ODE, match phase w*x,
# true-tail delta = E^2/(2 w x).  S(s) = spread over 5 zeros of phi_k=atan2(B,A) mod pi.
# S=0 <=> all five zeros co-spectral under one extension = Hilbert-Polya condition.
def phi(E,s):
    w=s*2*np.pi; FP=w*w; mu=-E*E
    def rhs(x,y): return [y[1],(2*x*y[1]+(FP*x*x-mu)*y[0])/(1-x*x)]
    xs=np.linspace(X-5,X,600)
    sol=solve_ivp(rhs,[1+1e-4,X],[1.0,(mu-FP)/2.0],rtol=1e-10,atol=1e-12,max_step=0.02,t_eval=xs)
    u=xs*sol.y[0]; ph=w*xs - E*E/(2*w*xs)
    M=np.column_stack([np.sin(ph),np.cos(ph)]); (A,B),*_=np.linalg.lstsq(M,u,rcond=None)
    return np.mod(np.arctan2(B,A),np.pi)
def spread(phis):  # min circular span (max-min after best rigid rotation on the pi-circle)
    p=np.sort(phis%np.pi); gaps=np.diff(np.append(p,p[0]+np.pi)); return np.pi-gaps.max()
print("S(s) build-constant scan, beta=0 slice, X=120.  s=1.0 is the current build (freq=2pi).")
print(" S=0 => Hilbert-Polya (all 5 zeros co-spectral).  Current ~2.2 rad.\n")
best=(9,None)
for s in np.linspace(0.60,1.60,21):
    phis=np.array([phi(g,s) for g in gam]); S=spread(phis)
    if S<best[0]: best=(S,s)
    print(f"  s={s:4.2f} (freq={s:.2f}*2pi): S={S:.4f} rad   phi={np.round(phis,3)}")
print(f"\nMIN over scalar frequency: S={best[0]:.4f} at s={best[1]:.3f}")
print("VERDICT: min_s S >1 rad => (C) scalar cannot cure, wrong-shape; <0.5 => (B) live.")

print("\n--- FINE refine 1.45..1.70, and X-stability check at the min ---")
for s in np.linspace(1.45,1.70,26):
    phis=np.array([phi(g,s) for g in gam]); print(f"  s={s:5.3f}: S={spread(phis):.4f}")
