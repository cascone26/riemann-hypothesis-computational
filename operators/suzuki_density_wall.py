import numpy as np, sys, os
sys.path.insert(0, os.path.join(os.getcwd(),'operators'))
from suzuki_exact_Q import build_Q, E_p, shat
from scipy.linalg import solve

G = np.array([float(x) for x in open("operators/_zeros150.txt").read().split()][:40])

def weyl_pieces(a, K, NG=400, lam_shift=0.2):
    kk=np.arange(1,K+1); Q=build_Q(a,K,NG); lam_a=np.linalg.eigvalsh(Q)[0]/a
    lam=lam_a-lam_shift; Ta=Q-lam*a*np.eye(K)
    cplus=solve(Ta,E_p(1.0,a,kk),assume_a='sym'); cminus=solve(Ta,E_p(-1.0,a,kk),assume_a='sym')
    npl=np.sqrt(abs(cplus@E_p(1.0,a,kk))); nmi=np.sqrt(abs(cminus@E_p(-1.0,a,kk))); cminus*= npl/nmi
    return cplus,cminus,kk

def W_zeros(cplus,cminus,kk,a,theta,zmax):
    def Wabs(z): return abs(np.sum(cplus*shat(z,a,kk))*(z-1j)+np.exp(1j*theta)*(z+1j)*np.sum(cminus*shat(z,a,kk)))
    zg=np.linspace(4,zmax,int(zmax*250)); wv=np.array([Wabs(z) for z in zg]); med=np.median(wv)
    zs=[zg[i] for i in range(1,len(zg)-1) if wv[i]<wv[i-1] and wv[i]<wv[i+1] and wv[i]<0.35*med]
    return np.array(zs)

def rvm_count(T):  # Riemann-von Mangoldt N(T)
    return T/(2*np.pi)*np.log(T/(2*np.pi)) - T/(2*np.pi) + 7/8 if T>7 else 0

print("="*94)
print("§986 DENSITY-WALL TEST on the correct construction: does W-zero count track ζ density N(T), or saturate?")
print("="*94)
print(f"{'a':>4} | {'zmax':>5} | {'#ζ-zeros in window':>18} | {'#W-zeros matched(<0.15)':>23} | {'RvM N(zmax)':>11}")
for a in [1.5,2.0,2.5,3.0]:
    K=max(90,int(35*a)+50); cplus,cminus,kk=weyl_pieces(a,K)
    zmax=min(G[-1], 20+15*a)  # window grows modestly with a
    nz=int(np.sum(G<=zmax))
    # best theta by matching the first few
    best=None
    for theta in np.linspace(0,2*np.pi,80,endpoint=False):
        zs=W_zeros(cplus,cminus,kk,a,theta,zmax)
        if len(zs)<3: continue
        matched=sum(1 for g in G[G<=zmax] if np.min(np.abs(zs-g))<0.15)
        if best is None or matched>best[0]: best=(matched,theta,len(zs))
    matched,theta,nW = best if best else (0,0,0)
    print(f"{a:4.1f} | {zmax:5.1f} | {nz:>18} | {matched:>23} | {rvm_count(zmax):11.1f}")
print("="*94)
print("Honest read: if #matched SATURATES (tracks only low zeros, << N(T)) as the window grows, the correct")
print("construction hits the SAME density-mismatch wall Kim 2607.24830 calls fundamental — the false walls")
print("(§982 fragility, §984 counts-zeros) are cleared but the TRUE (field) wall remains. If it tracks N(T), signal.")
