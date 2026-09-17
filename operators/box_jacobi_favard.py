#!/usr/bin/env python3
"""
cycle 331 — Favard/Jacobi construction of the candidate heig from the UN-CIRCULAR box moments.

cy330 gave m_k = sum_{gamma>0} box_gamma^{-k}, box_gamma = 1/4+gamma^2, computed from xi only (no zeros),
verified exact (2 m_1 = 2+euler-log 4pi). The box-moment MEASURE mu = sum_gamma delta_{x_gamma},
x_gamma = 1/box_gamma, has INFINITE mass (m_0 = #zeros = inf). The right FINITE measure is
    dnu = x dmu = sum_gamma x_gamma * delta_{x_gamma},   moments nu_j = int x^j dnu = m_{j+1}.
Its Jacobi matrix J_nu (Favard) is self-adjoint tridiagonal with spectrum {x_gamma} = {1/box_gamma}.
Hence J_nu^{-1} has spectrum {box_gamma} = {1/4+gamma^2} = boxEigen = the forge's heig.

THIS SCRIPT: build (a_n,b_n) of J_nu two ways and compare —
  (A) UN-CIRCULAR: Chebyshev/Stieltjes algorithm on the xi-computed moments nu_j=m_{j+1}. Ill-conditioned
      (moments decay ~box_1^{-k}); report how many (a_n,b_n) survive before precision death.
  (B) CIRCULAR TRUTH: Lanczos on the discrete nodes {x_gamma, weight x_gamma} from zeros2000.npy (stable).
Then: is J_nu POSITIVE (all b_n^2>0 <=> Hankel-PSD <=> RH) and does it show STRUCTURE (a_n,b_n asymptotics
that a xi/prime input could pin) vs pure knife-edge?  Verify-before-ship: (A) vs (B) agreement = the honest
reach of the un-circular route.
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 40
K = 8    # moments m_1..m_K  => nu_0..nu_{K-1}; Chebyshev gives Nrec=K//2=4 recurrence pairs
J = 8    # even-series order 2J=16 for xi Taylor (fast, conditioned — cy330-proven)

def xi(s):
    return mp.mpf('0.5') * s * (s - 1) * mp.power(mp.pi, -s/2) * mp.gamma(s/2) * mp.zeta(s)

# --- unconditional moments m_k from xi Taylor at s=1/2 (cy330 method) ---
xu = lambda u: xi(mp.mpf('0.5') + u)
tay = mp.taylor(xu, 0, 2 * J)
b = [tay[2 * j] for j in range(J + 1)]
def binom(n, r): return mp.factorial(n)/(mp.factorial(r)*mp.factorial(n-r))
c = [mp.mpf(0)]*(K+1)
for k in range(K+1):
    s = mp.mpf(0)
    for j in range(k, J+1):
        s += b[j]*binom(j, k)*mp.power(mp.mpf('0.25'), j-k)
    c[k] = s*mp.power(-1, k)
c0 = c[0]
d = [c[k]/c0 for k in range(K+1)]
p = [mp.mpf(0)]*(K+1)
for k in range(1, K+1):
    s = d[k]
    for j in range(1, k):
        s -= mp.mpf(j)/k*p[j]*d[k-j]
    p[k] = s
m = [mp.mpf(0)]*(K+1)          # m[k] = box moment (per-gamma)
for k in range(1, K+1):
    m[k] = k*(-p[k])
print("unconditional box moments m_k (from xi):")
for k in range(1, K+1):
    print(f"  m_{k} = {mp.nstr(m[k], 14)}")

# --- (A) Chebyshev/Stieltjes: moments nu_j=m_{j+1} -> Jacobi (a_n,b_n), high precision ---
# Wheeler/Chebyshev algorithm on ordinary moments. nu_j, j=0..2N-1 -> a_0..a_{N-1}, b_1..b_{N-1}.
nu = [m[j+1] for j in range(K)]     # nu_0..nu_{K-1}
Nrec = K//2
# sigma[k][l] table (Gautschi's Chebyshev algorithm for monic OPs from moments)
# a_n, b_n via: a_0=nu_1/nu_0, b_0=nu_0; sigma_{-1,l}=0, sigma_{0,l}=nu_l.
sig_prev = [mp.mpf(0)]*(2*Nrec)
sig_cur = [nu[l] for l in range(2*Nrec)] if 2*Nrec<=K else [nu[l] for l in range(K)]
# guard lengths
L = len(nu)
sig_cur = [nu[l] for l in range(L)]
sig_prev = [mp.mpf(0)]*L
aA = []; bA = []
aA.append(nu[1]/nu[0]); bA.append(nu[0])
a_prev = aA[0]
for n in range(1, Nrec):
    sig_new = [mp.mpf(0)]*L
    for l in range(n, L-n):
        sig_new[l] = sig_cur[l+1] - a_prev*sig_cur[l] - (bA[-1] if n>=1 else 0)*sig_prev[l]
    # a_n, b_n
    if sig_cur[n] == 0 or sig_new[n] == 0:
        break
    a_n = sig_new[n+1]/sig_new[n] - sig_cur[n]/sig_cur[n-1] if (n-1)>=0 else sig_new[n+1]/sig_new[n]
    b_n = sig_new[n]/sig_cur[n-1] if (n-1)>=0 else sig_new[n]/sig_cur[n]
    aA.append(a_n); bA.append(b_n)
    sig_prev = sig_cur; sig_cur = sig_new; a_prev = a_n

print(f"\n(A) UN-CIRCULAR Jacobi coeffs from xi-moments (Chebyshev algorithm, dps={mp.mp.dps}):")
print(f"{'n':>3} {'a_n (diag)':>22} {'b_n^2 (offdiag^2)':>24}")
for n in range(len(aA)):
    bn2 = bA[n]
    print(f"{n:>3} {mp.nstr(aA[n],10):>22} {mp.nstr(bn2,10):>24}  {'PSD' if bn2>0 else 'NEG/BREAKDOWN'}")

# --- (B) CIRCULAR TRUTH: Lanczos on nodes x_gamma=1/box, weight x_gamma, from zeros ---
try:
    G = np.load('/Users/scones/projects/riemann/operators/zeros2000.npy')
except Exception:
    G = np.load('zeros2000.npy')
x = 1.0/(0.25 + G.astype(float)**2)   # nodes (each gamma once = per-gamma measure)
w = x.copy()                          # dnu = x dmu  => weight = x_gamma
w = w/ w.sum()                        # normalize (scale irrelevant to a_n; scales b_0 only)
def lanczos(nodes, wt, Kk):
    q = np.sqrt(wt); qp = np.zeros_like(q); a = np.zeros(Kk); bb = np.zeros(Kk); bprev = 0.0
    for k in range(Kk):
        v = nodes*q; aa = np.dot(q, v); a[k] = aa
        v = v - aa*q - bprev*qp
        v = v - q*np.dot(q, v) - qp*np.dot(qp, v)   # reorth
        nb = np.linalg.norm(v)
        if k+1 < Kk: bb[k+1] = nb
        qp = q; q = v/nb if nb > 0 else v; bprev = nb
    return a, bb
aT, bT = lanczos(x, w, Nrec+2)
print(f"\n(B) CIRCULAR TRUTH Jacobi coeffs (Lanczos on {len(G)} zero-nodes, float64):")
print(f"{'n':>3} {'a_n (diag)':>16} {'b_n (offdiag)':>16} {'b_n^2':>16}")
for n in range(min(Nrec+2, 12)):
    print(f"{n:>3} {aT[n]:>16.8e} {bT[n]:>16.8e} {bT[n]**2:>16.8e}")

# --- compare A vs B where both exist ---
print("\nAGREEMENT (A un-circular vs B truth), a_n:")
for n in range(min(len(aA), Nrec)):
    auA = float(aA[n]); rel = abs(auA-aT[n])/max(abs(aT[n]),1e-30)
    print(f"  n={n}: A={auA:.6e}  B={aT[n]:.6e}  rel.err={rel:.2e}")
