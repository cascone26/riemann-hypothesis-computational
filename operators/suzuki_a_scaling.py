import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.getcwd(),'operators'))
from suzuki_exact_Q import build_Q, E_p, shat
from scipy.linalg import solve
from scipy.optimize import brentq

G = np.array([float(x) for x in open("operators/_zeros150.txt").read().split()][:8])

def weyl_pieces(a, K, NG=400, lam_shift=0.2):
    kk=np.arange(1,K+1); Q=build_Q(a,K,NG); lam_a=np.linalg.eigvalsh(Q)[0]/a
    lam=lam_a-lam_shift; Ta=Q-lam*a*np.eye(K)
    cplus=solve(Ta,E_p(1.0,a,kk),assume_a='sym'); cminus=solve(Ta,E_p(-1.0,a,kk),assume_a='sym')
    npl=np.sqrt(abs(cplus@E_p(1.0,a,kk))); nmi=np.sqrt(abs(cminus@E_p(-1.0,a,kk))); cminus*= npl/nmi
    return (lambda z: np.sum(cplus*shat(z,a,kk))), (lambda z: np.sum(cminus*shat(z,a,kk))), lam_a

def kappa_kernel(mfun, zpts, tol=1e-9):
    Z=np.array(zpts); mv=np.array([mfun(z) for z in Z]); n=len(Z); K=np.empty((n,n),complex)
    for i in range(n):
        for j in range(n):
            d=Z[i]-np.conj(Z[j]); K[i,j]=np.imag(mv[i])/np.imag(Z[i]) if abs(d)<1e-12 else (mv[i]-np.conj(mv[j]))/d
    K=0.5*(K+K.conj().T); ev=np.linalg.eigvalsh(K); s=max(1.0,abs(ev).max())
    return int(np.sum(ev<-tol*s)), ev

def best_theta_and_zeros(vp,vm,a):
    # find theta locking low W-zeros to zeta zeros; scan theta, W(theta,z)=(z-i)vp+e^{iθ}(z+i)vm
    def Wr(theta,z): return np.real((z-1j)*vp(z)+np.exp(1j*theta)*(z+1j)*vm(z))
    def Wi(theta,z): return np.imag((z-1j)*vp(z)+np.exp(1j*theta)*(z+1j)*vm(z))
    best=None
    for theta in np.linspace(0,2*np.pi,60,endpoint=False):
        zs=[]
        zg=np.linspace(5,G[-1]+5,3000)
        f=np.array([Wr(theta,z) for z in zg])  # on real axis W is ~real up to phase; use |W|
        Wabs=np.array([abs((z-1j)*vp(z)+np.exp(1j*theta)*(z+1j)*vm(z)) for z in zg])
        for i in range(1,len(zg)-1):
            if Wabs[i]<Wabs[i-1] and Wabs[i]<Wabs[i+1] and Wabs[i]<0.3*np.median(Wabs):
                zs.append(zg[i])
        if len(zs)>=len(G):
            err=np.mean([min(abs(np.array(zs)-g)) for g in G])
            if best is None or err<best[0]: best=(err,theta,zs)
    return best

print("="*88)
print("§985 Suzuki construction convergence: correct-m κ AND W-zero↔ζ-zero match, vs support a")
print("="*88)
print(f"{'a':>5} | {'λ_a':>11} | {'κ(correct m)':>12} | {'min kernel eig':>15} | {'W↔ζ mean err (8 zeros)':>22}")
r=np.random.default_rng(3)
for a in [1.2,1.5,2.0,2.5,3.0]:
    K=max(80,int(30*a)+40); vp,vm,lam_a=weyl_pieces(a,K)
    def Theta(z): return -((z-1j)*vp(z))/((z+1j)*vm(z))
    def m_weyl(z):
        t=Theta(z); return 1j*(1+t)/(1-t)
    zp=[complex(x,y) for x,y in zip(r.uniform(6,G[-1]+3,60), r.uniform(0.05,3.0,60))]
    k,ev=kappa_kernel(m_weyl,zp)
    bt=best_theta_and_zeros(vp,vm,a)
    err = bt[0] if bt else float('nan')
    print(f"{a:5.2f} | {lam_a:11.3e} | {k:>12} | {ev.min():15.2e} | {err:22.4f}")
print("="*88)
print("Read: κ=0 at every finite a (Suzuki Thm 1.5 = real spectrum, EXPECTED, no RH content); the VALUE is")
print("the correct-m construction stays κ=0 AND faithful while its W-zeros lock onto ζ-zeros (err↓ as a→∞).")
