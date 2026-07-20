#!/usr/bin/env python3
"""VECTOR D probe #1: does τ (or A+B) carry the hard-Lefschetz / sl₂ structure behind the Hodge Index
Theorem — the algebraic engine that makes the function-field intersection form definite?

Two crisp tests:
(1) HODGE–RIEMANN SIGN PATTERN. In an sl₂ Lefschetz structure of weight, the form on the weight-k primitive
    part has sign (−1)^k. So A+B's NEGATIVE eigenvalues should sit at ODD weights (odd frequency index n),
    positive at even — IF the structure is present. If instead the negative modes are just "low |n|"
    regardless of parity, there is NO Lefschetz sign structure.
(2) LEFSCHETZ SELF-ADJOINTNESS. Is there a natural raising operator L (frequency shift +1) that is
    SELF-ADJOINT w.r.t. the form (τ L = Lᵀ τ, i.e. L symmetric in the τ-inner-product)? Hodge theory needs
    the Lefschetz operator to be a τ-isometry-compatible map. Measure ||τL − Lᵀτ|| / ||τL||.
"""
import sys, os, math, json
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "external", "groskin_2607.02828", "code"))
import numpy as np
from flint import ctx
from arb_ldlt_certify import build_arb_tau

def mats(c,N,prec=2000):
    ctx.prec=prec; A,DIM=build_arb_tau(c,N,prec)
    T=np.array([[float(A[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
    return T,DIM

def AB_from_tau():  # we need A+B (arch) separately -> rebuild via the split assembly
    pass

def main():
    import importlib.util
    # reuse the split builder from arch_floor_vs_full
    spec=importlib.util.spec_from_file_location("afvf", os.path.join(os.path.dirname(__file__),"arch_floor_vs_full.py"))
    afvf=importlib.util.module_from_spec(spec); spec.loader.exec_module(afvf)
    N=64
    print(f"=== Vector D probe #1: hard-Lefschetz / sl₂ structure in A+B and τ (N={N}) ===")
    for c in [14,50,200]:
        Aa,Af,DIM=afvf.build_split(c,N,2000)  # Aa=A+B (arb_mat), Af=full τ
        AB=np.array([[float(Aa[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
        TAU=np.array([[float(Af[i,j].mid()) for j in range(DIM)] for i in range(DIM)])
        # TEST 1: parity of negative modes of A+B
        w,V=np.linalg.eigh(AB); neg=[k for k in range(len(w)) if w[k]<0]
        parities=[]
        for k in neg:
            n_dom=int(np.argmax(V[:,k]**2))-N
            parities.append(abs(n_dom)%2)
        odd=sum(parities); even=len(parities)-odd
        # TEST 2: Lefschetz self-adjointness of the +1 frequency shift w.r.t. τ
        L=np.zeros((DIM,DIM))
        for i in range(DIM-1): L[i+1,i]=1.0      # raise weight by 1
        comm=TAU@L - L.T@TAU
        selfadj_resid=np.linalg.norm(comm)/ (np.linalg.norm(TAU@L)+1e-30)
        print(f"\nc={c:4d}: A+B has {len(neg)} negative modes | dominant-freq parity: odd={odd} even={even}")
        print(f"   -> {'HODGE-RIEMANN pattern (all one parity)' if (odd==0 or even==0) else 'MIXED parity => NO Lefschetz sign structure'}")
        print(f"   Lefschetz self-adjointness ||τL−Lᵀτ||/||τL|| = {selfadj_resid:.3f}  "
              f"({'L is τ-self-adjoint (Lefschetz-compatible)' if selfadj_resid<0.05 else 'NOT τ-self-adjoint'}")
    print("\n=== VERDICT (D probe #1) ===")
    print("  If negative modes are MIXED parity AND the shift is not τ-self-adjoint, the naive")
    print("  hard-Lefschetz/sl₂ backbone is ABSENT in τ — Vector D needs the full F₁ geometric object,")
    print("  not a structure already hiding in the operator. (Informative negative.)")

if __name__=="__main__":
    main()
