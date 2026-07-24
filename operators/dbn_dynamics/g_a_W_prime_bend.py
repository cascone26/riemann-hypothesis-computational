#!/usr/bin/env python3
"""
g_a_W_prime_bend.py — LEDGER §112: the METRIC-ARTIFACT-FREE test of Cor 1.6's mechanism. Does the PRIME term bend the
W-comb toward the IRREGULAR zeta pattern? A momentum comb (arch-only) is REGULAR (low gap coefficient-of-variation CV
= std/mean). Zeta zeros are IRREGULAR (high CV). Test: is gap-CV(FULL W) > gap-CV(ARCH-only W), i.e. do primes inject
zeta-like irregularity? This is density-independent (CV normalizes out the comb spacing) — immune to the §111 false
positive. Built on the compact G_a framework (§110, reaches large a).
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

GAMMAS=np.array([14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,48.005151,49.773832])

def build_vhats(a, gfun, M=1500, lam=0.0):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h
    gv=np.array([gfun(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c)
    A=P@(Gm-lam*Nm)@P; A=(A+A.T)/2
    qp=(np.exp(a)-np.exp(x)); qm=-(np.exp(-a)-np.exp(-x))
    up=(np.linalg.pinv(A,rcond=1e-12)@(P@(sw*qp)))/sw
    um=(np.linalg.pinv(A,rcond=1e-12)@(P@(sw*qm)))/sw
    def cumint(u): v=np.zeros(M+1); v[1:]=np.cumsum((u[1:]+u[:-1])/2*h); return v
    vp=cumint(up); vm=cumint(um)
    def vhat(v,z): return np.sum(v*np.exp(1j*z*x)*w)
    return (lambda z: vhat(vp,z)),(lambda z: vhat(vm,z))

def Wzeros(vhp,vhm,theta,zgrid):
    Wz=np.abs(np.array([(z-1j)*vhp(z)+np.exp(1j*theta)*(z+1j)*vhm(z) for z in zgrid]))
    med=np.median(Wz); zs=[]
    for i in range(1,len(zgrid)-1):
        if Wz[i]<Wz[i-1] and Wz[i]<Wz[i+1] and Wz[i]<0.05*med:
            d=(Wz[i-1]-Wz[i+1])/(2*(Wz[i-1]-2*Wz[i]+Wz[i+1])+1e-30)
            zs.append(zgrid[i]+d*(zgrid[1]-zgrid[0]))
    return np.array([z for z in zs if z>1])

def gap_cv(zs):
    if len(zs)<4: return np.nan
    g=np.diff(np.sort(zs)); return np.std(g)/np.mean(g)

# zeta gap-CV in a comparable window (first ~10 zeros)
zeta_cv=gap_cv(GAMMAS)
print(f"Prime-bend test: does the PRIME term inject zeta-like gap-irregularity into the W-comb? (compact G_a)")
print(f"  reference: zeta gap-CV (first 10 γ) = {zeta_cv:.3f}  (irregular);  a pure momentum comb → CV≈0")
print(f"{'a':>5} {'CV(FULL W)':>12} {'CV(ARCH W)':>12} {'FULL/ARCH':>10}  (FULL>ARCH ⇒ primes bend toward zeta)")
zgrid=np.linspace(0.5,50,5000)
for a in [1.5,2.0,2.5,3.0]:
    # average CV over a few θ (the comb pattern shouldn't hinge on one θ)
    vhpF,vhmF=build_vhats(a,screw_g)
    vhpA,vhmA=build_vhats(a,screw_g_no_primes)
    cvF=np.nanmean([gap_cv(Wzeros(vhpF,vhmF,th,zgrid)) for th in np.linspace(0,2*np.pi,8,endpoint=False)])
    cvA=np.nanmean([gap_cv(Wzeros(vhpA,vhmA,th,zgrid)) for th in np.linspace(0,2*np.pi,8,endpoint=False)])
    print(f"{a:>5.1f} {cvF:>12.4f} {cvA:>12.4f} {cvF/cvA if cvA>1e-6 else np.inf:>10.2f}",flush=True)
print("\nREAD: FULL-W gap-CV > ARCH-W gap-CV ⇒ the prime term genuinely injects irregularity (bending the regular")
print("momentum comb toward the irregular zeta pattern) = Cor 1.6's mechanism, density-independent (no §111 artifact).")
