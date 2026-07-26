#!/usr/bin/env python3
"""ATTACK: the Weil form's QUOTIENT BY ITS RADICAL as a candidate arithmetic intersection form with Hodge signature.
Synthesis of P1 (pole = exact rank-2 sig (1,1)) + P2 (S∞). The knife-edge lives INSIDE the radical (near-null);
mod it out and the quotient should have a genuine MARGIN. Pole-split gives the positive class. Hodge-index-for-FF
proves RH via signature (1,-,-,...) on a NON-degenerate intersection form; here we test whether the quotient form
is that object.

High precision (mpmath) so the Schur / quotient are RELIABLE (float64 Schur was unstable).
Questions:
 (1) Quotient by radical: remove the k near-null directions (radical). Is the remaining form PD with a real GAP
     (margin restored once the knife-edge radical is gone)? What's the gap, and is it stable in c?
 (2) Signature with the pole: is the pole direction u the (or a) positive class? Is quotient-mod-pole neg-def w/ margin?
 (3) Does the radical DIMENSION stay bounded (→ finite-dim quotient = finite polarization) or grow with c?
"""
import importlib.util, os
import mpmath as mp
mp.mp.dps=60
GPATH=os.path.expanduser('~/projects/riemann/external/groskin_2607.02828/code/arb_ldlt_certify.py')
spec=importlib.util.spec_from_file_location('gc',GPATH); g=importlib.util.module_from_spec(spec); spec.loader.exec_module(g)
def to_mpf(x):
    s=str(x)
    if '+/-' in s: s=s.split('+/-')[0].strip('[] ')
    try: return mp.mpf(s)
    except Exception: return mp.mpf(float(x))

def build(c,N,prec=1400):
    A,DIM=g.build_arb_tau(c,N,prec)
    M=mp.matrix(DIM,DIM)
    for i in range(DIM):
        for j in range(DIM): M[i,j]=to_mpf(A[i,j])
    return M,DIM

for c in [13,29]:
    print(f"\n===== c={c} =====",flush=True)
    N=22
    tau,DIM=build(c,N)
    E=mp.eigsy(tau,eigvals_only=True); ev=sorted([mp.mpf(e) for e in E])
    mx=max(abs(ev[0]),abs(ev[-1]))
    # radical = eigenvalues below tol
    tol=mp.mpf(10)**(-8)*mx
    rad=[e for e in ev if abs(e)<tol]; nz=[e for e in ev if abs(e)>=tol]
    k=len(rad)
    gap = min(abs(e) for e in nz)
    print(f" dim={DIM}  lambda range [{mp.nstr(ev[0],4)}, {mp.nstr(ev[-1],4)}]")
    print(f" (3) RADICAL dim (|eig|<1e-8*max) = {k}   (grows with c? track 13 vs 29)")
    print(f" (1) QUOTIENT (radical removed): smallest |nonzero eig| = GAP = {mp.nstr(gap,6)}")
    print(f"     smallest 8 |eig| (the near-null tower): {[mp.nstr(abs(e),4) for e in sorted(ev,key=abs)[:8]]}")
    # (2) signature of the FULL tau: n_pos, n_neg, n_zero (Weil form is PSD so n_neg should be ~0; the radical=zero)
    npos=sum(1 for e in ev if e>tol); nneg=sum(1 for e in ev if e<-tol); nzero=k
    print(f" (2) signature of tau: (+{npos}, -{nneg}, 0:{nzero})   [PSD => n_neg=0; radical=the zeros]")
    # gap-to-radical ratio: is there a CLEAN spectral gap separating radical from the bulk?
    if nz:
        sep = gap/ (max(abs(e) for e in rad)+mp.mpf(10)**-50) if rad else mp.inf
        print(f"     spectral separation gap/max|radical| = {mp.nstr(sep,4)}  (large => clean finite-dim quotient)")
print("\nREAD: radical dim BOUNDED across c + a clean GAP separating it from the bulk => the quotient-by-radical is")
print("a finite-dim non-degenerate form = the finite polarization candidate. Then the Hodge/pole signature is the lever.")
