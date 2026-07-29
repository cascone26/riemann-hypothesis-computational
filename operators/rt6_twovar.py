#!/usr/bin/env python3
"""§410 (ATTACK: the second space via the genuine 2-variable Weil form). NOT a tensor, NOT random-equivalent:
the two variables are COUPLED through the prime pair-correlation (Hardy-Littlewood singular series S(h)),
the arithmetic content of the DOUBLE explicit formula. Build the 2-variable quadratic form on test pairs and
compute its SIGNATURE; discriminate (per §408 control) real zeros vs random/GUE.

Construction: 2-variable test K = sum_{jk} c_{jk} F_j(x)F_k(y). The double explicit formula's PRIME side is
  sum_{n,m} Lambda(n)Lambda(m)/sqrt(nm) * [ g(log n,log m) ]  where the CONNECTED (2-variable-irreducible)
part is weighted by the pair correlation. We build the CONNECTED kernel:
  B2[(j,k),(j',k')] = sum over zero PAIRS (g1,g2) with g1!=g2 of F_j(g1)F_k(g2) F_j'(g1)F_k'(g2)
     MINUS the disconnected (diagonal-product) part  =>  the genuine 'off-diagonal C x C' object.
Equivalently: B2 = (pairGram) - (B tensor B) restricted to connected; its indefiniteness (if RH-special)
is the second-space Hodge-index signature. Discriminate vs random/GUE zero surrogates."""
import numpy as np
import l2_edge_sweep as L
rng = np.random.default_rng(3)
w = 3.0
cs = L.probes(L.BAND, 4.0)   # sub-critical, n probes
n = len(cs)
def pairforms(gam):
    # F evaluated at each zero: E[i, a] = F_{c_i}(gamma_a)
    E = np.array([[np.exp(-((g-c)**2)/(2*w*w))+np.exp(-((g+c)**2)/(2*w*w)) for g in gam] for c in cs])  # n x Nz
    Nz = E.shape[1]
    # full pair Gram over ALL ordered pairs (a,b): B2full[(i,k),(i',k')] = sum_{a,b} E[i,a]E[k,b]E[i',a]E[k',b]
    #   = (sum_a E[i,a]E[i',a]) (sum_b E[k,b]E[k',b]) = B[i,i'] B[k,k']  -> tensor (disconnected)
    B = E @ E.T                                   # n x n  (= zero Gram, our 1-var form)
    # CONNECTED = restrict to a==b diagonal pairs MINUS ... actually the connected 2-pt is the a==b part:
    # sum_{a} E[i,a]E[k,a]E[i',a]E[k',a]  (both variables hit the SAME zero) = the 4-index diagonal contraction
    C = np.einsum('ia,ka,ja,la->ikjl', E, E, E, E).reshape(n*n, n*n)   # connected (coincident) part
    Bt = np.kron(B, B)                            # disconnected tensor
    return B, C, Bt
def sig(M, tol=1e-9):
    ev = np.linalg.eigvalsh((M+M.T)/2); t = tol*max(1, abs(ev).max())
    return int((ev>t).sum()), int((ev<-t).sum()), int((abs(ev)<=t).sum())
print(f"2-variable Weil form, band [14,44] w=3, s=4.0, n={n} probes ({n*n}-dim pair space).")
def report(name, gam):
    B, C, Bt = pairforms(gam)
    # the genuine connected 2-var object: does C - (something) have RH-special indefinite signature?
    # test three natural 2-var forms:
    print(f"  [{name}] connected C: {sig(C)} | tensor Bt: {sig(Bt)} | C - Bt/Nz: {sig(C - Bt/len(gam))}")
    return C, Bt
Creal, Btreal = report("REAL zeros", L.ZEROS)
# controls: random uniform + rigid grid, same count & range
for tr in range(3):
    report(f"random #{tr}", np.sort(rng.uniform(10,102,len(L.ZEROS))))
report("equispaced", np.linspace(10,102,len(L.ZEROS)))
print("\nreading: if the connected 2-var form C (or C-Bt/Nz) has a signature for REAL zeros that DIFFERS")
print("systematically from random/equispaced controls, the second space carries RH-content (the pair")
print("coupling is the surface). If identical to controls, this second-space construction is trivial too.")
