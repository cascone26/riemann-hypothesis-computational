"""
Verify the Connes-van Suijlekom engine mechanism = Pisarenko/Caratheodory-Fejer:
the eigenvector for the SIMPLE smallest eigenvalue of a Hermitian Toeplitz (autocorrelation) matrix
has its z-transform polynomial zeros ON the unit circle. (CS thm: simple+even min eig => real zeros of
Fourier transform = zeros on critical line via the map.)
"""
import numpy as np
np.random.seed(0)

def toeplitz_from_autocorr(r):
    n=len(r); T=np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            T[i,j]=r[abs(i-j)]
    return T

print("Pisarenko/CS engine check: min-eigenvector of Hermitian Toeplitz -> zeros on unit circle?")
print(f"{'case':>22} {'min eig simple?':>15} {'max |zero|-1|':>14} {'all on circle?':>14}")
for trial in range(5):
    # random positive-definite autocorrelation sequence (real, even) of order n
    n=7
    # build from a random spectral density (guarantees PD Toeplitz)
    m=200; freqs=np.linspace(-np.pi,np.pi,m)
    S=0.1+np.abs(np.random.randn(m))  # positive spectral density
    r=np.array([np.trapezoid(S*np.cos(k*freqs),freqs)/(2*np.pi) for k in range(n)])
    T=toeplitz_from_autocorr(r)
    ev,V=np.linalg.eigh(T)
    gap=ev[1]-ev[0]
    simple = gap>1e-9*ev[-1]
    v=V[:,0]                     # min-eigenvalue eigenvector
    # z-transform polynomial coeffs = v; zeros:
    roots=np.roots(v[::-1])      # polynomial v[0]+v[1]z+...+v[n-1]z^{n-1}
    dev=np.max(np.abs(np.abs(roots)-1)) if len(roots) else 0.0
    print(f"{'random PD Toeplitz #%d'%trial:>22} {str(simple):>15} {dev:14.2e} {str(dev<1e-6):>14}")
print("\nIf zeros lie on |z|=1 (dev~0) for simple min eig: CS/Pisarenko engine CONFIRMED in our hands.")
print("This is the mechanism Connes' program rides: finite Weil-Toeplitz min eigenvector -> on-line zeros.")
