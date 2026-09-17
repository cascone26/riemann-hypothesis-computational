"""Attack the WALL directly: κ = #negative eigenvalues of the Weil form Q=arch-prime (Suzuki).
RH-relevant: is κ robustly PINNED (by what mechanism), and can arch's PSD-ness + the prime structure
FORCE κ to stay bounded? Compute κ(Q), the negative spectrum, and how it scales with a,K."""
import numpy as np, sys
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q

print("κ(Q) = #neg eigenvalues of Weil form arch-prime, + smallest eigenvalues, vs (a,K):")
print(f"{'a':>4} {'K':>4} {'κ=#neg':>7} {'λ_min':>12} {'λ_2':>12} {'λ_3':>12}  {'arch alone κ':>12}")
for a in [0.5,1.0,1.5,2.0,2.5,3.0]:
    for K in [ max(60,int(30*a)+40) ]:
        Q,P=build_Q(a,K,NG=400,return_pieces=True)
        arch=P["pole"]+P["core"]+P["tail"]
        ev=np.linalg.eigvalsh(Q); eva=np.linalg.eigvalsh(arch)
        kap=int(np.sum(ev<-1e-9)); kap_arch=int(np.sum(eva<-1e-9))
        print(f"{a:4.1f} {K:4d} {kap:7d} {ev[0]:12.4e} {ev[1]:12.4e} {ev[2]:12.4e}  {kap_arch:12d}")
print()
# is the pole the culprit? compare arch-prime with pole REMOVED
print("Decompose: arch=pole+core+tail. Test κ of (core+tail-prime) [pole removed] vs full:")
for a in [1.5,2.5]:
    K=max(60,int(30*a)+40)
    Q,P=build_Q(a,K,NG=400,return_pieces=True)
    full=P["pole"]+P["core"]+P["tail"]+P["prime"]
    nopole=P["core"]+P["tail"]+P["prime"]
    poleonly=P["pole"]
    for lab,M in [("full arch-prime",full),("pole removed",nopole),("pole only",poleonly),("core+tail (arch smooth)",P["core"]+P["tail"])]:
        ev=np.linalg.eigvalsh((M+M.T)/2)
        print(f"  a={a} {lab:26s} κ={int(np.sum(ev<-1e-9)):3d}  λ_min={ev[0]:+.3e}  rank_neg_sig={int(np.sum(ev<-1e-9))}")
    print()
