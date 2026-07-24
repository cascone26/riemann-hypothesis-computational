#!/usr/bin/env python3
"""
cvs_branches.py — LEDGER §167: resolve the §166 CUSP into its two crossing eigenvalue-branches.

§166 found λ_min(t) has a sharp cusp at t=1 (asymmetric one-sided slopes) ⟹ the ground state is a
CROSSING of two smooth eigenvalue branches λ_a(t), λ_b(t) of the prime-scaled family Q(t). Goal: track
the two SMALLEST even-sector eigenvalues vs t. If they form an X (swap identity) at t=1, the crossing is
confirmed and λ_min(1) is the crossing height — reframing RH-positivity as "crossing height ≥ 0".

Also tracks which branch is the "zeta" branch (ground-state eigenvector at t=1, whose zeros are the γ's)
by following eigenvector overlap across t, to see if the zeta-branch is the one that becomes marginal.

Reuses the §166 component precompute (arch paid once). Usage: cvs_branches.py [N] [dps] [T] [c]
"""
import sys
import mpmath as mp
sys.path.insert(0, ".venv_cvs/lib/python3.14/site-packages")
import connes_cvs.operator as O

N   = int(sys.argv[1]) if len(sys.argv) > 1 else 20
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 50
T   = int(sys.argv[3]) if len(sys.argv) > 3 else 200
c   = int(sys.argv[4]) if len(sys.argv) > 4 else 13

mp.mp.dps = dps
L = mp.log(c)
prime_data, _ = O.prime_powers_up_to(int(c))
print(f"§167 CvS branch-crossing — c={c} N={N} dps={dps} T={T}")
print("track the TWO smallest even-sector eigenvalues vs prime-coupling t; is t=1 a crossing (X)?\n", flush=True)

idx = list(range(-N, N + 1))
pri={}; prid={}; pol={}; pold={}; arc={}; arcd={}
for n in idx:
    x = mp.mpf(n)
    pri[n]=O.psi_prime(x,L,prime_data); prid[n]=O.psi_prime_deriv(x,L,prime_data)
    pol[n]=O.psi_pole(x,L);            pold[n]=O.psi_pole_deriv(x,L)
    O._hplus_cache_clear()
    arc[n]=O.psi_arch(x,L,T,dps);      arcd[n]=O.psi_arch_deriv(x,L,T,dps)
    O._hplus_cache_clear()

DIM = 2*N+1
# even-sector projector (matches package): col0 = e_0, col k = (e_k+e_{-k})/sqrt2
V = mp.matrix(DIM, N+1); V[N,0]=mp.mpf(1); r2=1/mp.sqrt(2)
for k in range(1,N+1):
    V[N+k,k]=r2; V[N-k,k]=r2

def eig_even(t):
    """return sorted even-sector eigenpairs (eigs ascending, matching eigvecs as columns)."""
    psi ={n: t*pri[n]+pol[n]+arc[n]  for n in idx}
    psid={n: t*prid[n]+pold[n]+arcd[n] for n in idx}
    Q = mp.matrix(DIM,DIM)
    for i in range(DIM):
        m=i-N
        for j in range(DIM):
            nn=j-N
            Q[i,j]= psid[nn] if m==nn else (psi[m]-psi[nn])/(m-nn)
    for i in range(DIM):
        for j in range(i+1,DIM):
            a=(Q[i,j]+Q[j,i])/2; Q[i,j]=a; Q[j,i]=a
    Qe = V.T*Q*V
    ev, vc = mp.eigsy(Qe)
    order = sorted(range(N+1), key=lambda i: ev[i])
    eigs = [ev[i] for i in order]
    vecs = [[vc[r,i] for r in range(N+1)] for i in order]  # vecs[j] = j-th eigenvector (list)
    return eigs, vecs

def overlap(u, w):
    return abs(sum(u[i]*w[i] for i in range(len(u))))

# reference zeta-branch = ground eigenvector at t=1
e1, v1 = eig_even(mp.mpf(1))
zeta_vec = v1[0]
print(f"{'t':>10}   {'λ_0 (min)':>16}   {'λ_1 (2nd)':>16}   {'gap λ1-λ0':>14}   {'|<gs|ζ>|':>10}")
ts = [mp.mpf(s) for s in
      ["0.9","0.99","0.999","0.9999","1","1.0001","1.001","1.01","1.1"]]
for t in ts:
    eigs, vecs = eig_even(t)
    ov = overlap(vecs[0], zeta_vec)
    gap = eigs[1]-eigs[0]
    print(f"{float(t):>10.4f}   {mp.nstr(eigs[0],8):>16}   {mp.nstr(eigs[1],8):>16}   {mp.nstr(gap,6):>14}   {mp.nstr(ov,6):>10}", flush=True)

print("\nREAD:")
print(" - λ_0 and λ_1 approach then swap (gap→0 at t=1, X-shape) => genuine crossing; λ_min(1)=crossing height")
print(" - |<gs|ζ>| flips 1→0 across t=1 => the ground eigenvector changes identity => confirms cusp not soft min")
print(" - if gap stays open (avoided crossing) => it's a smooth min after all, reinterpret §166 cusp as resolution")
