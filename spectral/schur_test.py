"""SHARP test: positivity via Schur complement on the arch-negative block. In arch eigenbasis, sort so
M = A+B = [[N, C],[Cᵀ, P]], N=negative-mode block (arch_neg+Bneg ~ tiny diagonal residual), C=neg->pos
coupling, P=positive block (posdef). M⪰0 iff P⪰0 AND N - C P⁻¹ Cᵀ ⪰ 0. Compute the Schur complement
eigenvalues directly: is positivity a CLEAN Schur statement (proof strategy) or still delicate?"""
import numpy as np, sys
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q
for a in [1.0,1.5,2.0,2.5]:
    K=max(70,int(30*a)+40)
    Q,P_=build_Q(a,K,NG=400,return_pieces=True)
    A=(P_["core"]+P_["tail"]); A=(A+A.T)/2
    B=(P_["pole"]+P_["prime"]); B=(B+B.T)/2
    wa,Va=np.linalg.eigh(A); idx=np.argsort(wa); wa=wa[idx]; Va=Va[:,idx]
    M=Va.T@(A+B)@Va; M=(M+M.T)/2       # full form in sorted arch basis (diagonal A part = wa)
    kneg=int(np.sum(wa<-1e-9))
    N=M[:kneg,:kneg]; C=M[:kneg,kneg:]; Pp=M[kneg:,kneg:]
    lamP=np.linalg.eigvalsh(Pp)[0]
    Schur = N - C@np.linalg.solve(Pp, C.T)
    sev=np.linalg.eigvalsh((Schur+Schur.T)/2)
    print(f"a={a}: kneg={kneg}  λ_min(P_block)={lamP:+.4f}  λ_min(N)={np.linalg.eigvalsh(N)[0]:+.4f}")
    print(f"   Schur complement eigenvalues = {np.round(sev,4)}")
    print(f"   min Schur eig = {sev.min():+.4e}   ({'CLEAN POSITIVE (Schur proof works!)' if sev.min()>-1e-6 else 'NEGATIVE at this trunc -> still delicate/global'})")
    print(f"   ||C||={np.linalg.norm(C):.3f}  ||C||²/λ_min(P)={np.linalg.norm(C)**2/max(lamP,1e-9):.3f} (crude corr bound)")
