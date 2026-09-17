import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.getcwd(),'operators'))
from suzuki_exact_Q import build_Q, E_p, shat
from scipy.linalg import solve

def weyl_pieces(a, K, NG=400, lam_shift=0.2):
    kk=np.arange(1,K+1); Q=build_Q(a,K,NG); lam_a=np.linalg.eigvalsh(Q)[0]/a
    lam=lam_a-lam_shift; Ta=Q-lam*a*np.eye(K)
    fplus=E_p(1.0,a,kk); fminus=E_p(-1.0,a,kk)
    cplus=solve(Ta,fplus,assume_a='sym'); cminus=solve(Ta,fminus,assume_a='sym')
    npl=np.sqrt(abs(cplus@fplus)); nmi=np.sqrt(abs(cminus@fminus)); cminus=cminus*(npl/nmi)
    def vp(z): return np.sum(cplus*shat(z,a,kk))
    def vm(z): return np.sum(cminus*shat(z,a,kk))
    return vp, vm, lam_a

def kappa_kernel(mfun, zpts, tol=1e-9):
    Z=np.array(zpts); mv=np.array([mfun(z) for z in Z]); n=len(Z); K=np.empty((n,n),complex)
    for i in range(n):
        for j in range(n):
            d=Z[i]-np.conj(Z[j]); K[i,j]=np.imag(mv[i])/np.imag(Z[i]) if abs(d)<1e-12 else (mv[i]-np.conj(mv[j]))/d
    K=0.5*(K+K.conj().T); ev=np.linalg.eigvalsh(K); s=max(1.0,abs(ev).max())
    return int(np.sum(ev<-tol*s)), ev

for a in [1.5, 2.0]:
    K=max(80,int(30*a)+40)
    vp,vm,lam_a=weyl_pieces(a,K)
    # Suzuki inner function Theta (|Theta|=1 on R); Weyl m = Cayley transform i(1+Theta)/(1-Theta)
    def Theta(z): return -((z-1j)*vp(z))/((z+1j)*vm(z))
    def m_weyl(z):
        th=Theta(z); return 1j*(1+th)/(1-th)
    def m_naive(z): return vp(z)/vm(z)
    print("="*84); print(f"a={a}  (ω>1 built regime),  λ_a={lam_a:.4e}")
    # CONTROL A: inner function |Theta|->1 as Im z ->0+
    for y in [0.5,0.1,0.01]:
        vals=[abs(Theta(complex(x,y))) for x in np.linspace(6,40,12)]
        print(f"  |Θ| at Im={y:<4}: mean={np.mean(vals):.4f} (→1 on R ⇒ genuine inner fn)")
    # CONTROL B: real-symmetry m(conj z)=conj(m(z))?
    zt=[complex(x,y) for x,y in [(12,0.6),(20,1.1),(31,0.4)]]
    symn=np.mean([abs(m_naive(np.conj(z))-np.conj(m_naive(z))) for z in zt])
    symw=np.mean([abs(m_weyl(np.conj(z))-np.conj(m_weyl(z))) for z in zt])
    print(f"  real-symmetry defect |m(z̄)-conj m(z)|:  naive={symn:.3e}   weyl-Cayley={symw:.3e}")
    # kappa via corrected kernel, UHP test points near axis, in-window
    r=np.random.default_rng(3)
    zp=[complex(x,y) for x,y in zip(r.uniform(6,40,60), r.uniform(0.05,3.0,60))]
    kn,_=kappa_kernel(m_naive,zp); kw,evw=kappa_kernel(m_weyl,zp)
    print(f"  κ(corrected kernel):  naive m=vp/vm → {kn}    Weyl-Cayley m → {kw}")
    print(f"    (Weyl-Cayley kernel eig range: min={evw.min():.2e} max={evw.max():.2e})")

print("="*84)
print("FAITHFULNESS CONTROL: inject off-line conjugate-pair pole(s) into the correct Weyl m — does κ track?")
a=1.5; K=max(80,int(30*a)+40); vp,vm,lam_a=weyl_pieces(a,K)
def Theta(z): return -((z-1j)*vp(z))/((z+1j)*vm(z))
def m_weyl(z):
    th=Theta(z); return 1j*(1+th)/(1-th)
r=np.random.default_rng(3)
zp=[complex(x,y) for x,y in zip(r.uniform(6,40,60), r.uniform(0.05,3.0,60))]
for npairs in [0,1,2,3]:
    centers=[15.0,22.0,31.0][:npairs]; d=0.4
    def m_inj(z, C=centers):
        val=m_weyl(z)
        for c in C: val += 1.0/((c+1j*d)-z) + 1.0/((c-1j*d)-z)
        return val
    k,_=kappa_kernel(m_inj,zp)
    print(f"  {npairs} injected off-line pair(s) → κ = {k}   (expect {npairs})")
# robustness of the base κ=0 under prime-height jitter analog: perturb lam_shift
print("  base κ=0 stability under lam_shift ∈ {0.15,0.2,0.25,0.3}:")
for ls in [0.15,0.2,0.25,0.3]:
    vp2,vm2,_=weyl_pieces(a,K,lam_shift=ls)
    def Th2(z): return -((z-1j)*vp2(z))/((z+1j)*vm2(z))
    def mw2(z):
        t=Th2(z); return 1j*(1+t)/(1-t)
    k,_=kappa_kernel(mw2,zp); print(f"    lam_shift={ls} → κ={k}")
