#!/usr/bin/env python3
"""
semilocal_jacobi.py — LEDGER §119: compute the SEMILOCAL JACOBI MATRIX that Connes-Moscovici 2310.18423 explicitly
DEFERRED (line 258). The cyclic pair (𝒮, ξ_S) has spectral measure dm_S(s) = |Γ(1/4+is/2)|² ∏_{p∈S}|L_p(1/2-is)|² ds,
L_p(1/2-is)=(1-p^{-(1/2-is)})^{-1}. The orthogonal polynomials for dm_S give the Jacobi matrix (diag α_n=0 by evenness,
off-diag √β_n). This shows HOW the primes shift the archimedean recurrence coefficients (the arch↔prime structure,
§104/§112). VALIDATE: arch-only (S={∞}) must give CM's stated √((n+1/2)(n+1)) (line 227-8), then add primes.
"""
import mpmath as mp
mp.mp.dps = 30

def weight(s, S):
    # arch: |Gamma(1/4+is/2)|^2
    w = abs(mp.gamma(mp.mpf(1)/4 + 1j*s/2))**2
    for p in S:
        if p=='inf': continue
        # |L_p(1/2-is)|^2 = 1/(1 - 2 p^{-1/2} cos(s log p) + 1/p)
        pm=mp.mpf(p)
        w *= 1/(1 - 2*pm**mp.mpf('-0.5')*mp.cos(s*mp.log(pm)) + 1/pm)
    return w

def moments(S, K, L=60):
    # even measure: mu[2k] = 2 int_0^L s^{2k} w ds ; odd=0
    m=[mp.mpf(0)]*(2*K+1)
    for k in range(0, 2*K+1, 2):
        m[k]=2*mp.quad(lambda s: s**k * weight(s,S), [0, 2, 6, 15, 30, L])
    return m

def stieltjes(m, N):
    # recurrence coeffs from moments via inner products <p_i,p_j>=sum p_i[a]p_j[b] m[a+b]
    def ip(p,q):
        return mp.fsum(p[a]*q[b]*m[a+b] for a in range(len(p)) for b in range(len(q)))
    # p_0=1, p_-1=0 ; p_{n+1}(s)=(s-alpha_n)p_n - beta_n p_{n-1}
    p_prev=[mp.mpf(0)]; p_cur=[mp.mpf(1)]
    betas=[]; norms=[ip(p_cur,p_cur)]
    for n in range(N):
        # alpha_n = <s p_n, p_n>/<p_n,p_n>  (=0 for even measure; compute anyway)
        sp=[mp.mpf(0)]+p_cur[:]                      # s*p_n
        alpha=ip(sp,p_cur)/norms[-1]
        beta = norms[-1]/norms[-2] if n>=1 else mp.mpf(0)
        betas.append(beta)
        # p_{n+1} = (s-alpha)p_n - beta p_{n-1}
        sp_a=[c for c in sp]
        for i in range(len(p_cur)): sp_a[i]-=alpha*p_cur[i]
        pn1=sp_a[:]
        for i in range(len(p_prev)): pn1[i]-=beta*p_prev[i]
        p_prev,p_cur=p_cur,pn1
        norms.append(ip(p_cur,p_cur))
    return betas   # beta_n ; off-diagonal Jacobi element = sqrt(beta_n)

print("SEMILOCAL JACOBI (CM 2310.18423 deferred). Off-diagonal √β_n vs the CM arch value √((n+1/2)(n+1)).")
N=6
configs=[('arch {∞}',['inf']), ('{∞,2}',['inf',2]), ('{∞,2,3}',['inf',2,3]), ('{∞,2,3,5}',['inf',2,3,5])]
print(f"{'n':>3} " + " ".join(f"{'√β_'+str(n):>12}" for n in range(1,N)))
cmval=[mp.sqrt((n+mp.mpf(1)/2)*(n+1)) for n in range(N)]
print(f"{'CM arch':>3} " + " ".join(f"{mp.nstr(cmval[n],7):>12}" for n in range(1,N)))
for name,S in configs:
    m=moments(S,N+1)
    b=stieltjes(m,N)
    print(f"{name:>8} " + " ".join(f"{mp.nstr(mp.sqrt(b[n]),7):>12}" for n in range(1,N)), flush=True)
print("\nREAD: (1) does arch-only reproduce CM's √((n+1/2)(n+1))? (validates). (2) how do primes SHIFT √β_n — up? down?")
print("structured? This is the concrete arch↔prime modification of the operator, the piece CM left uncomputed.")
