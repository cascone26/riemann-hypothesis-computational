"""Engage arXiv:2602.20313 (dBN kernel PF-order). Reproduce Theorem 1.1 (5x5 Toeplitz det<0
at u0=0.01,h=0.05, det=-1.847e-9) to validate the kernel, then compute C_5(u0) sign/threshold
(the OPEN, withdrawn-as-uncertified piece; claimed C_5<0 near 0, sign change near u0~0.0311)."""
import mpmath as mp
mp.mp.dps = 50
def Phi(u):
    # eq (2): sum_{n>=1} (2 pi^2 n^4 e^{9u} - 3 pi n^2 e^{5u}) e^{-pi n^2 e^{4u}}
    u = mp.mpf(u)
    tot = mp.mpf(0)
    e4u = mp.e**(4*u); e9u = mp.e**(9*u); e5u = mp.e**(5*u)
    for n in range(1, 60):
        term = (2*mp.pi**2*n**4*e9u - 3*mp.pi*n**2*e5u)*mp.e**(-mp.pi*n**2*e4u)
        tot += term
        if abs(term) < mp.mpf(10)**(-45) and n > 3: break
    return tot
def K(u): return Phi(abs(mp.mpf(u)))
# Theorem 1.1 reproduction
u0, h = mp.mpf('0.01'), mp.mpf('0.05')
M = mp.matrix(5,5)
for i in range(5):
    for j in range(5):
        M[i,j] = K(u0 + (i-j)*h)
detM = mp.det(M)
print("REPRODUCE Thm 1.1: det(M_5) at u0=0.01,h=0.05 =", mp.nstr(detM, 8))
print("  (paper: -1.847846073e-9; match?)")
# lower-order dets D2,D3,D4 (Thm 1.2: positive)
for r in (2,3,4,5):
    Mr = mp.matrix(r,r)
    for i in range(r):
        for j in range(r): Mr[i,j]=K(u0+(i-j)*h)
    print(f"  D_{r}(0.01,0.05) = {mp.nstr(mp.det(Mr),6)}")
print("DONE")
