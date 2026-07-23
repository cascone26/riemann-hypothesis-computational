#!/usr/bin/env python3
"""§192: VERIFY the §191 parity decoupling, cleanly. Basis φ_k=sin(kπ(x+a)/2a) (k=1..K): φ_k EVEN in x iff k ODD,
ODD iff k EVEN. Even kernel ⟹ form is parity-block-diagonal. Split Qf,Qa,Mass into the ODD-k block (=EVEN-v
sector) and EVEN-k block (=ODD-v sector); report per-sector arch #negative and full λ_min. Confirms: which sector
holds the 2.5e-5 ground state, and that EACH sector has exactly ONE arch negative mode (⟹ per-sector scalar secular)."""
import sys, math
import mpmath as mp
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from mpmath_weil_hiprec import weil_matrices, gen_eig_sym
mp.mp.dps=40; K=24; M=600
def subblock(Q, idx):
    n=len(idx); B=mp.matrix(n,n)
    for i in range(n):
        for j in range(n): B[i,j]=Q[idx[i],idx[j]]
    return B
for a in [mp.mpf('0.5'), mp.mpf('0.549')]:
    Qf,Qa,Mass,Vv,xs,w = weil_matrices(a,K,M)
    # k index (1-based actual mode) for column c is c+1; EVEN-v sector = odd actual-k = even column-index c (0-based c: k=c+1 odd ⟺ c even)
    evenV = [c for c in range(K) if (c+1)%2==1]   # k odd  -> v EVEN
    oddV  = [c for c in range(K) if (c+1)%2==0]   # k even -> v ODD
    # verify block-diagonality: off-block coupling norm
    offA = mp.fsum(abs(Qa[i,j]) for i in evenV for j in oddV)
    print(f"\n=== a={mp.nstr(a,4)}  (off-block |Qa| coupling = {mp.nstr(offA,3)} — should be ~0 if parity-diagonal) ===")
    for name, idx in [("EVEN-v (k odd)", evenV), ("ODD-v (k even)", oddV)]:
        Ea,_=gen_eig_sym(subblock(Qa,idx), subblock(Mass,idx))
        Ef,_=gen_eig_sym(subblock(Qf,idx), subblock(Mass,idx))
        nneg=sum(1 for e in Ea if e<0)
        print(f"  {name:>16}: arch #neg={nneg}  arch λ1={mp.nstr(Ea[0],5)}  |  full λmin={mp.nstr(Ef[0],5)}")
print("\nREAD: the sector with full λmin≈2.5e-5 is the CRITICAL one; if each sector has arch #neg=1 ⟹ per-sector scalar")
print("secular is exact. Confirms/CORRECTS §191's sector assignment (which parity is critical).")
