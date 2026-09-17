"""The pinned mechanism: structure of the prime coupling B=prime(+pole) in the ARCHIMEDEAN eigenbasis.
Primes act by SHIFTS (log n) => in the smooth-arch eigenbasis, is B Toeplitz/banded (a convolution), and
does its off-diagonal profile reveal WHY the global cancellation holds? Look for legible structure."""
import numpy as np, sys
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q
a=1.5; K=85
Q,P=build_Q(a,K,NG=400,return_pieces=True)
A=(P["core"]+P["tail"]); A=(A+A.T)/2
B=(P["pole"]+P["prime"]); B=(B+B.T)/2
wa,Va=np.linalg.eigh(A)
# order arch eigenmodes by eigenvalue (most-negative first)
idx=np.argsort(wa); wa=wa[idx]; Va=Va[:,idx]
Bt=Va.T@B@Va
print(f"a={a}, K={K}. Arch eigenvalues sorted (most neg first): {np.round(wa[:8],2)} ... {np.round(wa[-3:],2)}")
# off-diagonal decay by band |i-j|
print("\nOff-diagonal RMS of B (arch eigenbasis) by band distance d=|i-j|:")
for d in [0,1,2,3,5,8,12,20,40]:
    vals=[Bt[i,i+d] for i in range(K-d)]
    print(f"  d={d:3d}: RMS={np.sqrt(np.mean(np.array(vals)**2)):.4f}")
# is it Toeplitz? compare diag-band variance
print("\nToeplitz check (is B_{i,i+d} ~ const along each band?): coeff of variation per band:")
for d in [1,2,3,5]:
    vals=np.array([Bt[i,i+d] for i in range(K-d)])
    cv=vals.std()/(abs(vals.mean())+1e-9)
    print(f"  d={d}: mean={vals.mean():+.3f} std={vals.std():.3f} CV={cv:.2f}")
# KEY: the negative arch modes are the low-index ones. Is B's coupling AMONG them structured?
kneg=int(np.sum(wa<-1e-9))
print(f"\nB restricted to arch-NEGATIVE block ({kneg}x{kneg}, the modes needing lift):")
Bneg=Bt[:kneg,:kneg]
print(f"  diag(Bneg) = {np.round(np.diag(Bneg),2)}")
print(f"  |arch_neg| = {np.round(-wa[:kneg],2)}   (diag(Bneg) should ~ match to cancel)")
print(f"  Bneg off-diag RMS = {np.sqrt(np.mean((Bneg-np.diag(np.diag(Bneg)))**2)):.3f}")
# coupling from neg block to POSITIVE modes (how the lift is 'sourced')
Bnp=Bt[:kneg,kneg:]
print(f"  neg->pos coupling ||Bnp|| = {np.linalg.norm(Bnp):.2f}  (vs ||Bneg||={np.linalg.norm(Bneg):.2f})")
print(f"  => is the positivity sourced by neg-block internal structure or neg->pos coupling?")
