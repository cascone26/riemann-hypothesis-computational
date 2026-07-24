#!/usr/bin/env python3
"""
g_a_W.py — LEDGER §111: W(a,θ;z) rebuilt on the COMPACT G_a framework (§110), to test Cor 1.6 (W→z²ξ/ξ′) into the
LARGE-a regime the sine-Galerkin A_a couldn't reach (§42/§45 stopped ~a=2, match 0.025).
Construction (compact, well-conditioned): the deficiency solutions v_± of (A_a−λ)v=e^{±x} become, in u=v′ (mean-zero):
    (G − λN) u_± = q_±,   q_±(t) = ±(e^{±a} − e^{±t}),   v_±(x)=∫_{−a}^x u_±,
with G (kernel g(x−y)), N (inverse Neumann Laplacian, explicit continuous kernel). λ<λ_a (take λ=0 since λ_a>0, §110).
Then v̂_±(z)=∫v_± e^{izx}dx and W(a,θ;z)=(z−i)v̂_+ + e^{iθ}(z+i)v̂_-. Zeros of W are real (Thm 1.5) = spectrum of D_{a,θ};
Cor 1.6: with θ=θ(a), they → γ (zeta zeros) as a→∞ ⟹ RH. Test: does the match keep tightening past a=2 now?
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

GAMMAS = np.array([14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073])

def build_vhats(a, M=1500, lam=0.0):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    X=x[:,None]; Y=x[None,:]
    N=(X**2+Y**2)/(4*a)-np.abs(X-Y)/2+a/6
    sw=np.sqrt(w)
    Gm=sw[:,None]*G*sw[None,:]; Gm=(Gm+Gm.T)/2
    Nm=sw[:,None]*N*sw[None,:]; Nm=(Nm+Nm.T)/2
    c=sw/np.linalg.norm(sw); P=np.eye(M+1)-np.outer(c,c)
    A=P@(Gm-lam*Nm)@P; A=(A+A.T)/2
    # rhs q_±(t)=±(e^{±a}-e^{±t}); in weighted, mean-zero
    qp=(np.exp(a)-np.exp(x)); qm=-(np.exp(-a)-np.exp(-x))
    qpw=P@(sw*qp); qmw=P@(sw*qm)
    # solve A u = q (pinv on mean-zero range)
    up=np.linalg.pinv(A,rcond=1e-12)@qpw; um=np.linalg.pinv(A,rcond=1e-12)@qmw
    up=up/sw; um=um/sw                                   # unweight -> u_±(x_i)  (guard sw>0 interior)
    # v_±(x)=∫_{-a}^x u  (cumulative trapezoid), enforce v(-a)=0
    def cumint(u):
        v=np.zeros(M+1); v[1:]=np.cumsum((u[1:]+u[:-1])/2*h); return v
    vp=cumint(up); vm=cumint(um)
    def vhat(v,z): return np.sum(v*np.exp(1j*z*x)*w)
    return (lambda z: vhat(vp,z)),(lambda z: vhat(vm,z))

def Wfun(vhp,vhm,theta,z): return (z-1j)*vhp(z)+np.exp(1j*theta)*(z+1j)*vhm(z)

def zeros_of_W(vhp,vhm,theta,zgrid):
    Wz=np.abs(np.array([Wfun(vhp,vhm,theta,z) for z in zgrid]))
    med=np.median(Wz); zs=[]
    for i in range(1,len(zgrid)-1):
        if Wz[i]<Wz[i-1] and Wz[i]<Wz[i+1] and Wz[i]<0.05*med:
            # parabolic refine
            d=(Wz[i-1]-Wz[i+1])/(2*(Wz[i-1]-2*Wz[i]+Wz[i+1])+1e-30)
            zs.append(zgrid[i]+d*(zgrid[1]-zgrid[0]))
    return zs

def match_to_gamma(zs):
    zs=[z for z in zs if z>1]
    if not zs: return None
    return np.mean([min(abs(g-np.array(zs))) for g in GAMMAS[:5]])

print("W(a,θ;z) on the COMPACT G_a framework — Cor 1.6 test into large a (θ-scan, match of low W-zeros to γ_1..γ_5):")
print(f"{'a':>5} {'best θ':>8} {'mean|γ−W| (5 zeros)':>20} {'first 5 W-zeros':>34}")
zgrid=np.linspace(0.5,45,4500)
for a in [1.0,1.5,2.0,2.5,3.0]:
    vhp,vhm=build_vhats(a)
    best=(1e9,None,None)
    for theta in np.linspace(0,2*np.pi,25,endpoint=False):
        zs=zeros_of_W(vhp,vhm,theta,zgrid)
        m=match_to_gamma(zs)
        if m is not None and m<best[0]: best=(m,theta,zs)
    m,th,zs=best
    z5=[round(z,2) for z in sorted([z for z in (zs or []) if z>1])[:5]]
    print(f"{a:>5.1f} {th:>8.3f} {m:>20.4f} {str(z5):>34}",flush=True)
print("\nREAD: does mean|γ−W| keep DROPPING past a=2 (where sine-Galerkin stalled at 0.025)? Tightening → Cor 1.6")
print("convergence W→z²ξ/ξ′ demonstrated into the zeta regime the old operator couldn't reach.")
