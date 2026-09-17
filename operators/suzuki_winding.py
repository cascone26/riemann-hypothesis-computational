import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.getcwd(),'operators'))
from suzuki_exact_Q import build_Q, E_p, shat
from scipy.linalg import solve
def weyl_pieces(a,K,NG=400,lam_shift=0.2):
    kk=np.arange(1,K+1);Q=build_Q(a,K,NG);lam_a=np.linalg.eigvalsh(Q)[0]/a
    Ta=Q-(lam_a-lam_shift)*a*np.eye(K)
    cp=solve(Ta,E_p(1.0,a,kk),assume_a='sym');cm=solve(Ta,E_p(-1.0,a,kk),assume_a='sym')
    cm*=np.sqrt(abs(cp@E_p(1.0,a,kk)))/np.sqrt(abs(cm@E_p(-1.0,a,kk)))
    return (lambda z:np.sum(cp*shat(z,a,kk))),(lambda z:np.sum(cm*shat(z,a,kk)))
def winding(f,C):
    v=np.array([f(z) for z in C]);return (np.unwrap(np.angle(v))[-1]-np.unwrap(np.angle(v))[0])/(2*np.pi)
def uhp(R,ymax,n=6000):
    eps=0.02;x=np.linspace(-R,R,n);y=np.linspace(eps,ymax,n//4)
    return np.concatenate([x+1j*eps,R+1j*y,np.linspace(R,-R,n)+1j*ymax,-R+1j*y[::-1]])
print("="*82); print("§987 winding κ handle — match to kernel across a + faithfulness under UHP injection")
print("="*82)
C=uhp(44,8)
print(f"{'a':>4} | {'winding[Θ-1]=κ':>15} | {'winding[Θ]≈zero-count':>21}")
for a in [1.2,1.5,2.0,2.5,3.0]:
    K=max(80,int(30*a)+40);vp,vm=weyl_pieces(a,K)
    def Th(z):return -((z-1j)*vp(z))/((z+1j)*vm(z))
    print(f"{a:4.1f} | {round(winding(lambda z:Th(z)-1,C)):>15} | {winding(Th,C):21.1f}")
print("-"*82)
a=1.5;K=80;vp,vm=weyl_pieces(a,K)
def Th(z):return -((z-1j)*vp(z))/((z+1j)*vm(z))
print("FAITHFULNESS: inject k off-line (UHP) zeros into (Θ-1) via Blaschke (z-uₖ)/(z-ūₖ), uₖ∈UHP:")
for k in range(4):
    us=[15+0.8j,22+1.0j,31+0.9j][:k]
    def g(z,us=us):
        val=Th(z)-1.0
        for u in us: val*= (z-u)/(z-np.conj(u))
        return val
    print(f"  k={k} injected UHP zeros → winding[(Θ-1)·∏] = {round(winding(g,C))}  (expect {k})")
print("="*82)
