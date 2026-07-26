# Claim: the Cauchy corner B_kl = 1/(k+l), k in [25,M] odd, l > M odd, has
#   (a) scale-invariant norm (M-independent),
#   (b) PROVABLE bound ||B|| <= ||B||_HS, and HS^2 = sum 1/(k+l)^2 -> (1/4)log2  (odd spacing),
#       via continuum substitution y=1/s mapping corner to kernel 1/(1+ts) on (0,1)^2, HS^2 = log2.
# Also measure the full u-law coupling corner for comparison with Ein's 0.19-0.24.
import numpy as np, math
b = math.pi/(2*0.549); th = b*math.log(2); beta = 0.624052; c1 = 1.517
def u(m): return -1.0 - beta*np.sin(th*m) + c1/m
for M in (259, 801):
    k = np.arange(25, M+1, 2, dtype=float)
    l = np.arange(M+2, 120001, 2, dtype=float)
    K, L = np.meshgrid(k, l, indexing='ij')
    C = 1.0/(K+L)
    hs2 = (C**2).sum()
    nC = np.linalg.norm(C, 2)
    # full u-law coupling (Loewner form)
    B = (K*u(L) - L*u(K))/((K-L)*(K+L))
    nB = np.linalg.norm(B, 2)
    print(f"M={M:4d}: Cauchy corner norm={nC:.4f}  HS={math.sqrt(hs2):.4f}  (bound sqrt(log2)/2={math.sqrt(math.log(2))/2:.4f})  u-law ||B||={nB:.4f}  ||B||^2={nB*nB:.4f}")
