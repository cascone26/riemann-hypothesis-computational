#!/usr/bin/env python3
"""
cvs_effective.py — LEDGER §167b: the 2×2 effective prime-operator at the §167 crossing, and the
Hellmann–Feynman verification of the §166 cusp slopes.

Q(t) = Q_arch+pole + t·Q_prime, so dQ/dt = Q_prime (matrix built from psi_prime alone). At t=1 the two
lowest even-sector eigenvalues are a near-null pair (§167). The §166 cusp has one-sided slopes; by
Hellmann–Feynman each branch slope = ⟨ψ|Q_prime|ψ⟩. The two crossing states diagonalize Q_prime WITHIN
the 2D near-null subspace, so the 2×2 matrix M = [⟨u_i|Q_prime|u_j⟩] (u_0,u_1 = the near-null pair at t=1)
should have eigenvalues equal to the measured cusp slopes (§166: ≈ +0.70 left, −2.35 right at N=14).

This (a) VERIFIES the cusp is exactly the prime-operator's action in the near-null subspace, and (b) gives
the effective object whose structure governs the crossing height. Usage: cvs_effective.py [N] [dps] [T] [c]
"""
import sys
import mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
import connes_cvs.operator as O

N   = int(sys.argv[1]) if len(sys.argv) > 1 else 18
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
T   = int(sys.argv[3]) if len(sys.argv) > 3 else 200
c   = int(sys.argv[4]) if len(sys.argv) > 4 else 13

mp.mp.dps = dps
L = mp.log(c)
prime_data, _ = O.prime_powers_up_to(int(c))
print(f"§167b CvS effective 2×2 / Hellmann–Feynman — c={c} N={N} dps={dps} T={T}\n", flush=True)

idx = list(range(-N, N+1))
pri={}; prid={}; pol={}; pold={}; arc={}; arcd={}
for n in idx:
    x = mp.mpf(n)
    pri[n]=O.psi_prime(x,L,prime_data); prid[n]=O.psi_prime_deriv(x,L,prime_data)
    pol[n]=O.psi_pole(x,L);            pold[n]=O.psi_pole_deriv(x,L)
    O._hplus_cache_clear()
    arc[n]=O.psi_arch(x,L,T,dps);      arcd[n]=O.psi_arch_deriv(x,L,T,dps)
    O._hplus_cache_clear()

DIM=2*N+1
V=mp.matrix(DIM,N+1); V[N,0]=mp.mpf(1); r2=1/mp.sqrt(2)
for k in range(1,N+1):
    V[N+k,k]=r2; V[N-k,k]=r2

def build(psi, psid):
    Q=mp.matrix(DIM,DIM)
    for i in range(DIM):
        m=i-N
        for j in range(DIM):
            nn=j-N
            Q[i,j]= psid[nn] if m==nn else (psi[m]-psi[nn])/(m-nn)
    for i in range(DIM):
        for j in range(i+1,DIM):
            a=(Q[i,j]+Q[j,i])/2; Q[i,j]=a; Q[j,i]=a
    return V.T*Q*V

# full operator at t=1  and  the pure prime operator (= dQ/dt)
Qe1  = build({n: pri[n]+pol[n]+arc[n] for n in idx}, {n: prid[n]+pold[n]+arcd[n] for n in idx})
Qpe  = build(pri, prid)   # prime piece alone, even sector

ev, vc = mp.eigsy(Qe1)
order = sorted(range(N+1), key=lambda i: ev[i])
print("two smallest even-sector eigenvalues of the TRUE operator at t=1 (the near-null pair):")
print(f"  λ_0 = {mp.nstr(ev[order[0]],8)}   λ_1 = {mp.nstr(ev[order[1]],8)}\n")

# 2D near-null subspace
u = [[vc[r,order[k]] for r in range(N+1)] for k in range(2)]
def mv(A,x):  # A(N+1 x N+1) times vector x
    return [sum(A[r,cc]*x[cc] for cc in range(N+1)) for r in range(N+1)]
def dot(a,b): return sum(a[i]*b[i] for i in range(N+1))

M = mp.matrix(2,2)
for i in range(2):
    Qpu = mv(Qpe, u[i])
    for j in range(2):
        M[i,j] = dot(u[j], Qpu)
M[0,1]=M[1,0]=(M[0,1]+M[1,0])/2
mev = mp.eigsy(M, eigvals_only=True)
print("2×2 effective PRIME operator M = ⟨u_i|Q_prime|u_j⟩ in the near-null subspace:")
print(f"  M = [[{mp.nstr(M[0,0],6)}, {mp.nstr(M[0,1],6)}],")
print(f"       [{mp.nstr(M[1,0],6)}, {mp.nstr(M[1,1],6)}]]")
print(f"  eigenvalues of M (= the two branch slopes dλ/dt at the crossing): {mp.nstr(mev[0],6)}, {mp.nstr(mev[1],6)}")
print("\nVERIFY: these two eigenvalues should match the §166 one-sided cusp slopes (~+0.70 and −2.35 at small N).")
print("If so, the §166 cusp IS exactly the prime operator's action in the near-null subspace (Hellmann–Feynman).")
