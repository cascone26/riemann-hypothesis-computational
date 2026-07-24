#!/usr/bin/env python3
"""
weyl_mfunction.py — LEDGER §113: the WEYL m-FUNCTION of the Suzuki/Connes construction, the potential PROOF mechanism
(§107d). m_a(z) := v̂_+(z)/v̂_-(z), where v_± are the deficiency solutions (compact G_a build, §110/§111).
Weyl theory: eigenvalues of D_{a,θ} = {z : m_a(z) = -e^{iθ}(z+i)/(z-i)} (m meets the b.c. Möbius curve). m_a is
HERGLOTZ (self-adjoint problem ⟹ Im m_a>0 in C+, real poles/zeros interlacing on R). Cor 1.6 ⟺ m_a converges to a
limiting Herglotz m_∞ whose meeting-set with the b.c. = Γ (zeta zeros). KEY: Herglotz-limit STABILITY ⟹ if m_a→m_∞
locally uniformly with m_a Herglotz, then m_∞ Herglotz ⟹ its spectrum REAL ⟹ RH. So the reality transfers for FREE
once convergence holds. TEST: (1) is m_a actually Herglotz (Im m_a>0 in C+)? (2) do its POLES (zeros of v̂_-) converge
/ stabilize toward the zeta zeros γ as a grows? (3) is the convergence m_a(z) locally-uniform-looking?
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
GAM=np.array([14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073])

def vhat_pm(a, M=1500, lam=0.0):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]; N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w); Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2; Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c); A=P@(Gm-lam*Nm)@P; A=(A+A.T)/2
    qp=(np.exp(a)-np.exp(x)); qm=-(np.exp(-a)-np.exp(-x))
    Ai=np.linalg.pinv(A,rcond=1e-12)
    up=(Ai@(P@(sw*qp)))/sw; um=(Ai@(P@(sw*qm)))/sw
    def cumint(u): v=np.zeros(M+1); v[1:]=np.cumsum((u[1:]+u[:-1])/2*h); return v
    vp=cumint(up); vm=cumint(um)
    def vhat(v,z): return np.sum(v*np.exp(1j*z*x)*w)
    return (lambda z: vhat(vp,z)),(lambda z: vhat(vm,z))

def real_zeros(f, zg):
    F=np.array([f(z) for z in zg]); zs=[]
    for i in range(len(zg)-1):
        a0,b0=F[i],F[i+1]
        if a0.real*b0.real<0 and abs(a0.imag)+abs(b0.imag)<10*(abs(a0.real)+abs(b0.real)):
            zs.append(zg[i]-a0.real*(zg[i+1]-zg[i])/(b0.real-a0.real))
    return zs

print("WEYL m-function m_a(z)=v̂_+/v̂_-  (compact G_a). Herglotz check + pole convergence toward γ.")
zg=np.linspace(1,46,4500)
prev_poles=None
for a in [1.5,2.0,2.5,3.0]:
    vhp,vhm=vhat_pm(a)
    m=lambda z: vhp(z)/vhm(z)
    # (1) Herglotz: sample Im m_a at z = t+0.3i for t in [2,40]; should be >0 (Herglotz maps C+->C+)
    ts=np.linspace(2,40,200); imvals=np.array([m(t+0.3j).imag for t in ts])
    frac_pos=np.mean(imvals>0)
    # (2) poles of m_a = zeros of v̂_- (real); do they approach γ?
    poles=np.array(sorted([z for z in real_zeros(vhm,zg) if z>1]))
    # match poles to γ (nearest), report mean dist for first 5 γ, and pole spacing-CV
    if len(poles)>=4:
        dgam=np.mean([min(abs(g-poles)) for g in GAM[:5]])
        cv=np.std(np.diff(poles))/np.mean(np.diff(poles))
    else: dgam,cv=np.nan,np.nan
    drift=(np.nan if prev_poles is None or len(poles)<3 or len(prev_poles)<3
           else np.mean(np.abs(poles[:min(len(poles),len(prev_poles))]-prev_poles[:min(len(poles),len(prev_poles))])))
    print(f"  a={a}: Herglotz Im(m)>0 frac={frac_pos:.2f} | #poles={len(poles)} | poles[:6]={[round(p,2) for p in poles[:6]]}",flush=True)
    print(f"         nearest-γ mean(5)={dgam:.3f}  pole-gap-CV={cv:.3f}  drift-from-prev-a={drift}",flush=True)
    prev_poles=poles
print("\nREAD: (1) Im(m_a)>0 fraction ~1 ⇒ m_a is HERGLOTZ (the free-reality property). (2) do poles STABILIZE with a")
print("(drift→0) and carry zeta-like gap-CV? Stabilizing Herglotz poles = the Weyl-convergence proof-shape being real.")
