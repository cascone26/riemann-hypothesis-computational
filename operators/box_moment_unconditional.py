#!/usr/bin/env python3
"""
UN-CIRCULAR box-value moment problem (cycle 330, 2026-08-30).

Goal: compute m_k = sum_rho 1/(rho(1-rho))^k = sum_rho box_rho^{-k},  box_rho = 1/4 + gamma^2,
      WITHOUT any knowledge of zero locations -- straight from the completed zeta xi(s).

Identity (unconditional, from the Hadamard product):
  xi(1/2 + u) is EVEN in u  =>  F(w) := xi(1/2 + sqrt(1/4 - w)) is single-valued analytic in w,
  with zeros exactly at w = box_rho = 1/4 + gamma^2 = rho(1-rho).
  Hadamard:  F(w) = C * prod_rho (1 - w/box_rho)
  =>  -log( F(w)/F(0) ) = sum_{k>=1} (w^k / k) * sum_rho box_rho^{-k} = sum_k (w^k/k) m_k
  =>  m_k = k * [w^k]( -log F(w) ).

The Taylor coefficients of F at w=0 are computed by numerical differentiation of xi (mpmath),
which needs NOTHING about where the zeros are. So the resulting m_k are UNCONDITIONAL.

Then: Hankel PSD of (m_1,m_2,...) <=> {box_rho} is a positive real multiset <=> RH.
(Off-line zero => complex box => complex/failing moments => Hankel not PSD.)

Cross-check: also compute m_k directly from the known zeros (the CIRCULAR cy329 version) and
confirm the unconditional extraction reproduces them -- validates the identity/code, not RH.
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 40

def xi(s):
    # completed zeta: xi(s) = 1/2 s(s-1) pi^{-s/2} Gamma(s/2) zeta(s)
    return mp.mpf('0.5') * s * (s - 1) * mp.power(mp.pi, -s/2) * mp.gamma(s/2) * mp.zeta(s)

K = 6   # number of moments (Hankel n<=3 needs m_1..m_5)
J = 8   # even-series truncation in u (order 2J=16; high-j weight (1/4)^{j-k} tiny)

# --- Taylor coeffs of xi(1/2+u) at u=0 (s=1/2, NO pole); default method, modest order 2J ---
xu = lambda u: xi(mp.mpf('0.5') + u)
tay = mp.taylor(xu, 0, 2 * J)   # default finite-diff; fast & conditioned at order 2J=24
b = [tay[2 * j] for j in range(J + 1)]   # even coefficients b_j = [u^{2j}] xi(1/2+u)
print(f"xi(1/2) = {mp.nstr(b[0], 12)}   b_1=[u^2] = {mp.nstr(b[1], 8)}")

# --- F(w) = xi(1/2+u), u^2 = 1/4 - w  =>  F(w) = sum_j b_j (1/4 - w)^j ; collect [w^k] ---
# c_k = sum_{j>=k} b_j * C(j,k) * (-1)^k * (1/4)^{j-k}
def binom(n, r):
    return mp.factorial(n) / (mp.factorial(r) * mp.factorial(n - r))
c = [mp.mpf(0)] * (K + 1)
for k in range(K + 1):
    s = mp.mpf(0)
    for j in range(k, J + 1):
        s += b[j] * binom(j, k) * mp.power(mp.mpf('0.25'), j - k)
    c[k] = s * mp.power(-1, k)
c0 = c[0]
print(f"F(0)=xi(1) = {mp.nstr(c0, 12)}  (should be 0.5)")

# --- series -log(F/F0):  let g(w)=F/F0 = 1 + sum_{k>=1} d_k w^k,  d_k=c_k/c0 ---
d = [c[k] / c0 for k in range(K + 1)]   # d[0]=1
# L_k = [w^k]( -log g ).  Use recursion for log of a power series.
# log g = sum_{k>=1} L'_k w^k with standard recursion; we want -log g so negate.
# Newton: for g=1+sum d_k w^k, let p_k=[w^k] log g:
#   p_k = d_k - (1/k) sum_{j=1}^{k-1} j p_j d_{k-j}
p = [mp.mpf(0)] * (K + 1)
for k in range(1, K + 1):
    s = d[k]
    for j in range(1, k):
        s -= mp.mpf(j) / k * p[j] * d[k - j]
    p[k] = s
# -log F contributes L_k = -p_k ; m_k = k * L_k
m_unc = [mp.mpf(0)] * (K + 1)
for k in range(1, K + 1):
    m_unc[k] = k * (-p[k])

print("\nUNCONDITIONAL box moments m_k = k*[w^k](-log F)  (from xi only, NO zeros):")
for k in range(1, K + 1):
    print(f"  m_{k} = {mp.nstr(m_unc[k], 15)}")

# --- CROSS-CHECK: m_k directly from known zeros (circular; validates the identity) ---
try:
    G = np.load('/Users/scones/projects/riemann/operators/zeros2000.npy')
except Exception:
    G = np.load('zeros2000.npy')
print(f"\nCross-check from {len(G)} known zeros (box=1/4+gamma^2, both +-gamma => factor 2):")
boxz = 0.25 + G.astype(float)**2
for k in range(1, K + 1):
    tot = 2.0 * np.sum(boxz**(-k))   # factor 2: conjugate pair shares box
    diff = float(m_unc[k]) - tot
    print(f"  k={k}: from_zeros={tot:.15g}  unc-zeros={diff:.4g}")

# --- Hankel PSD test on the UNCONDITIONAL moments ---
# moment sequence for a measure needs m_0 too; m_0 = # zeros = infinite, so test the
# Hankel of the SHIFTED sequence (m_1,m_2,...) which is the Hamburger matrix H_n[i,j]=m_{i+j+1}.
print("\nHankel PSD test on UNCONDITIONAL m_k (H_n[i,j] = m_{i+j+1}):")
for n in range(1, K // 2 + 1):
    H = mp.matrix(n, n)
    for i in range(n):
        for j in range(n):
            H[i, j] = m_unc[i + j + 1]
    ev = mp.eigsy(H, eigvals_only=True)
    mineig = min(ev)
    det = mp.det(H)
    print(f"  n={n}: min_eig={mp.nstr(mineig, 6)}  det={mp.nstr(det, 6)}  {'PSD' if mineig > 0 else 'NOT PSD'}")
