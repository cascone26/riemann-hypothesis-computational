"""Nail the EXACT conditions for 'Toeplitz min-eigenvector zeros on unit circle' before formalizing."""
import numpy as np
np.random.seed(1)
def toep(r):
    n=len(r); return np.array([[r[abs(i-j)] for j in range(n)] for i in range(n)])
def zeros_dev(v):
    roots=np.roots(v[::-1]); return np.max(np.abs(np.abs(roots)-1)) if len(roots) else 0.0

print("A) smooth-density PD Toeplitz, various n — MIN eigenvector:")
for n in [4,6,8,12]:
    m=400; th=np.linspace(-np.pi,np.pi,m); S=0.2+np.abs(np.random.randn(m))
    r=np.array([np.trapezoid(S*np.cos(k*th),th)/(2*np.pi) for k in range(n)])
    T=toep(r); ev,V=np.linalg.eigh(T)
    print(f"   n={n:3d}: min-eig zeros dev={zeros_dev(V[:,0]):.1e}  2nd-eig dev={zeros_dev(V[:,1]):.1e}  top-eig dev={zeros_dev(V[:,-1]):.1e}")

print("\nB) RANDOM real-symmetric Toeplitz (r_k iid, NOT from a density) — is it still PD & on-circle?")
for t in range(4):
    n=6; r=np.random.randn(n); r[0]=abs(r[0])+n  # diagonally dominant → PD
    T=toep(r); ev,V=np.linalg.eigh(T)
    pd = ev[0]>0
    print(f"   trial {t}: PD={pd}  min-eig zeros dev={zeros_dev(V[:,0]):.1e}")

print("\nC) does it need MIN specifically? (compare all eigenvectors, n=6 smooth):")
n=6; m=400; th=np.linspace(-np.pi,np.pi,m); S=0.2+np.abs(np.random.randn(m))
r=np.array([np.trapezoid(S*np.cos(k*th),th)/(2*np.pi) for k in range(n)])
T=toep(r); ev,V=np.linalg.eigh(T)
for i in range(n):
    print(f"   eig#{i} (λ={ev[i]:.3f}): zeros dev={zeros_dev(V[:,i]):.1e}")
