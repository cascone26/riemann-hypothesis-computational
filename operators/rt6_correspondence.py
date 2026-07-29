#!/usr/bin/env python3
"""§408 (meta A1'/A4): the SECOND SPACE via PRIME CORRESPONDENCES. Zeros define a vector v_gamma;
prime n acts by U_n = multiplication by n^{i*gamma} = e^{i gamma log n} (translation by log n in the dual).
Build the form Q[n,m] = <U_n w, U_m w>_metric over the zero-sampled space, with (A) the PLAIN Hermitian
metric (should be PSD Gram) and (B) the FE-TWISTED metric g -> 1-g i.e. pairing U_n against the FE-partner
(gamma -> -gamma reflection / rho->1-rho). If the FE-twisted correspondence form is INDEFINITE where the
plain one is PSD, THAT indefiniteness is the second-space (Hodge-index-like) structure we've been missing.
Also check the algebra: do the U_n commute (abelian, dead) or is there nonabelian content?"""
import numpy as np
import l2_edge_sweep as L
Z = L.ZEROS
NM = 60
lam = np.zeros(NM+1); comp = np.zeros(NM+1, dtype=bool)
for p in range(2, NM+1):
    if not comp[p]:
        comp[2*p::p]=True; lp=np.log(p); pk=p
        while pk<=NM: lam[pk]=lp; pk*=p
primes = [n for n in range(2, NM+1) if lam[n] > 0]
w = np.array([lam[n]/np.sqrt(n) for n in primes])       # prime weights
logn = np.array([np.log(n) for n in primes])
# zero-sampled "state": weight each zero by a window (finite support), unit for now
gam = Z
# (A) plain Hermitian correspondence Gram: Q_A[a,b] = sum_gamma e^{i gamma (logn_a - logn_b)}  (= zero FT at freq diff)
D = logn[:,None]-logn[None,:]
Q_A = np.real(np.array([[np.sum(np.exp(1j*gam*(logn[a]-logn[b]))) for b in range(len(primes))] for a in range(len(primes))]))
# (B) FE-twisted: pair U_n(gamma) against U_m(-gamma) [rho->1-rho reflects gamma->-gamma on the Re=1/2 line's imag]
Q_B = np.real(np.array([[np.sum(np.exp(1j*gam*(logn[a]+logn[b]))) for b in range(len(primes))] for a in range(len(primes))]))
def sig(M, tol=1e-9):
    ev=np.linalg.eigvalsh((M+M.T)/2); t=tol*max(1,abs(ev).max())
    return int((ev>t).sum()), int((ev<-t).sum()), int((abs(ev)<=t).sum())
print(f"{len(primes)} primes to {NM}. Zero-sampled correspondence forms:")
print(f"  (A) plain U_n vs U_m  [freq diff log n - log m]: signature {sig(Q_A)}  lmin {np.linalg.eigvalsh((Q_A+Q_A.T)/2)[0]:+.3e}")
print(f"  (B) FE-twisted U_n vs U_m(-g) [freq SUM log n+log m]: signature {sig(Q_B)}  lmin {np.linalg.eigvalsh((Q_B+Q_B.T)/2)[0]:+.3e}")
# weighted by prime weights (the actual Weil form): W Q W
WA = np.outer(w,w)*Q_A; WB = np.outer(w,w)*Q_B
print(f"  (A) weighted (Lambda/sqrt n): signature {sig(WA)}")
print(f"  (B) weighted FE-twisted: signature {sig(WB)}")
# nonabelian check: do U_n as operators on zero-space commute?
# represent U_n on a truncated Fourier basis e^{i k a gamma}; commutator norm
K = 12; a = 0.3
basis_freqs = np.arange(-K, K+1)*a
def Uop(ln):
    # U_n shifts frequency by ln in the e^{i f gamma} basis -> permutation-ish; approx via Gram in this basis
    return np.diag(np.exp(1j*basis_freqs*ln))
U2, U3 = Uop(np.log(2)), Uop(np.log(3))
comm = np.linalg.norm(U2@U3 - U3@U2)
print(f"  [U_2,U_3] commutator norm on diagonal rep = {comm:.3e} (0 = abelian = dead per W4; >0 = nonabelian)")
print("\nreading: if (B) FE-twisted is INDEFINITE (neg eigenvalues) while (A) is PSD, the second-space")
print("indefinite structure EXISTS in the prime-correspondence/FE-twist. If both PSD or (B) trivially")
print("degenerate, the naive correspondence collapses too (W1 confirmed deep, needs true arithmetic surface).")
