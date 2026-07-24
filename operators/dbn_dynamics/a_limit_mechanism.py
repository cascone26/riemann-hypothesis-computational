#!/usr/bin/env python3
"""
a_limit_mechanism.py — LEDGER §122 fork (A): does the arch↔prime cancellation have a POSITIVITY LAW at the limit?
Push the compact G_a (§110, best-conditioned tool) to LARGE a and track:
  - λ_a(FULL) = smallest eigenvalue of the full screw Weil form  (stays >0? rate?)
  - λ_a(ARCH) = arch-only (should be O(-e^a), the negativity the primes must cancel)
  - the PRIME OFFSET = how much the primes lift λ_a above the arch value, and whether the margin (λ_a full) has a
    clean law. Mechanism test: is the prime-positivity provably ≥ the arch-negativity with a characterizable margin?
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

def lam_compact(a, gfun, M=1400):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([gfun(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); Gp=P@Gm@P; Np=P@Nm@P
    ev_N,V_N=np.linalg.eigh((Np+Np.T)/2); keep=ev_N>1e-12*ev_N.max()
    S=V_N[:,keep]/np.sqrt(ev_N[keep]); Ar=S.T@((Gp+Gp.T)/2)@S
    return np.linalg.eigvalsh((Ar+Ar.T)/2)[0]

print("Fork (A): arch↔prime positivity at the limit (compact G_a, large a).")
print(f"{'a':>5} {'λ_a(FULL)':>14} {'λ_a(ARCH)':>14} {'prime offset':>14} {'FULL>0?':>8}")
res=[]
for a in [2.0,2.5,3.0,3.5,4.0,4.5,5.0]:
    lf=lam_compact(a, screw_g); la=lam_compact(a, screw_g_no_primes)
    res.append((a,lf,la))
    print(f"{a:>5.1f} {lf:>14.6e} {la:>14.6e} {lf-la:>14.6e} {'YES' if lf>0 else 'NO':>8}",flush=True)
aa=np.array([r[0] for r in res]); lf=np.array([r[1] for r in res])
# rate of λ_a(FULL): log-slope (doubly-exp would show as steepening; simple-exp as constant slope)
sl=np.polyfit(aa[-4:], np.log(np.abs(lf[-4:])), 1)[0]
print(f"\n λ_a(FULL) log-slope (last 4): {sl:+.3f}  (stays >0 = RH-consistent into large a; slope = decay rate of the margin)")
print(" MECHANISM: arch is O(-e^a) negative; primes lift it to a small POSITIVE margin. Is that margin provably >0?")
print(" (= the √x-strength PNT cancellation; §37 meta-theorem: positivity ⟺ reality, likely no non-circular handle.)")
