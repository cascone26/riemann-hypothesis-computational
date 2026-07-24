#!/usr/bin/env python3
"""
cvs_collectivity.py — LEDGER §168: the §156 COLLECTIVITY test on the exact CvS operator.

§156 (sine-Galerkin Weil operator) found the arch-negative → positive lift is DECISIVELY COLLECTIVE: no single
prime dominates; positivity emerges from all primes together. Re-run on the frontier object. For c=13 the complete
prime-power set is {2,3,4,5,7,8,9,11,13} (9 entries) — small enough to test EXACTLY and completely.

Two probes (prime contributions are additive, so each prime-power's contribution is precomputed once):
  (A) CUMULATIVE — add prime-powers one at a time (ascending), λ_min after each. Smooth collective climb from the
      arch-only negative toward the near-null, or one jump attributable to a single prime?
  (B) LEAVE-ONE-OUT — remove each single prime-power from the full set, measure λ_min degradation. Comparable
      degradations ⟹ collective; one dominant ⟹ attributable.

Reuses the §166 arch/pole precompute. Usage: cvs_collectivity.py [N] [dps] [T] [c]
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
pp_labels = [str(int(pd[0])) for pd in prime_data]
print(f"§168 CvS collectivity — c={c} N={N} dps={dps} T={T}")
print(f"prime-powers ≤ {c}: {pp_labels}  ({len(prime_data)} entries)\n", flush=True)

idx = list(range(-N, N+1))
pol={}; pold={}; arc={}; arcd={}
# per-prime-power contribution vectors (value, deriv), computed once each
ppval  = [dict() for _ in prime_data]
ppderv = [dict() for _ in prime_data]
for n in idx:
    x = mp.mpf(n)
    pol[n]=O.psi_pole(x,L); pold[n]=O.psi_pole_deriv(x,L)
    O._hplus_cache_clear(); arc[n]=O.psi_arch(x,L,T,dps); arcd[n]=O.psi_arch_deriv(x,L,T,dps); O._hplus_cache_clear()
    for i, pd in enumerate(prime_data):
        ppval[i][n]  = O.psi_prime(x, L, [pd])
        ppderv[i][n] = O.psi_prime_deriv(x, L, [pd])

DIM=2*N+1
V=mp.matrix(DIM,N+1); V[N,0]=mp.mpf(1); r2=1/mp.sqrt(2)
for k in range(1,N+1): V[N+k,k]=r2; V[N-k,k]=r2

def lam_min(subset):
    """λ_min with prime piece = sum over prime-power indices in `subset` (arch+pole always included)."""
    psi={}; psid={}
    for n in idx:
        pv = sum((ppval[i][n]  for i in subset), mp.mpf(0))
        pd = sum((ppderv[i][n] for i in subset), mp.mpf(0))
        psi[n]  = pv + pol[n]  + arc[n]
        psid[n] = pd + pold[n] + arcd[n]
    Q=mp.matrix(DIM,DIM)
    for i in range(DIM):
        m=i-N
        for j in range(DIM):
            nn=j-N; Q[i,j]=psid[nn] if m==nn else (psi[m]-psi[nn])/(m-nn)
    for i in range(DIM):
        for j in range(i+1,DIM):
            a=(Q[i,j]+Q[j,i])/2; Q[i,j]=a; Q[j,i]=a
    lam,_=O.compute_ground_state(Q); return lam

full = list(range(len(prime_data)))

print("(A) CUMULATIVE — add prime-powers ascending; λ_min after each:")
print(f"  {'(arch+pole only)':>22}: {mp.nstr(lam_min([]),8)}")
for k in range(1, len(full)+1):
    lm = lam_min(full[:k])
    print(f"  {'+'.join(pp_labels[:k]):>22}: {mp.nstr(lm,8)}", flush=True)

print("\n(B) LEAVE-ONE-OUT — remove one prime-power from the full set; λ_min (full incl. all = near-null):")
lm_full = lam_min(full)
print(f"  {'FULL':>10}: {mp.nstr(lm_full,8)}")
for i in range(len(prime_data)):
    sub = [j for j in full if j != i]
    lm = lam_min(sub)
    print(f"  {'drop '+pp_labels[i]:>10}: {mp.nstr(lm,8)}   (Δ from full: {mp.nstr(lm-lm_full,6)})", flush=True)

print("\nREAD: smooth cumulative climb + comparable leave-one-out drops ⟹ COLLECTIVE (matches §156).")
print("A single dominant jump / one large leave-one-out drop ⟹ attributable to that prime.")
