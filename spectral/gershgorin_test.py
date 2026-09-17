"""Test whether Weil positivity is LOCAL (diagonally dominant in arch's eigenbasis -> Gershgorin-provable)
or GLOBAL (delicate cancellation). In arch's eigenbasis: A diagonal (a_i), B=prime+pole. Row i:
  Gershgorin lower bound  λ_i >= (a_i + B_ii) - Σ_{j≠i}|B_ij|.
If min over i of this is >= 0, positivity PROVEN locally. Else it's a global cancellation."""
import numpy as np, sys
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q
for a in [1.0,1.5,2.0]:
    K=max(70,int(30*a)+40)
    Q,P=build_Q(a,K,NG=400,return_pieces=True)
    A=(P["core"]+P["tail"]); A=(A+A.T)/2
    B=(P["pole"]+P["prime"]); B=(B+B.T)/2
    wa,Va=np.linalg.eigh(A)
    Bt=Va.T@B@Va                       # B in A's eigenbasis
    diag = wa + np.diag(Bt)            # a_i + B_ii  (the mode-by-mode residual)
    offrow = np.sum(np.abs(Bt),axis=1)-np.abs(np.diag(Bt))   # Σ_{j≠i}|B_ij|
    gersh = diag - offrow
    lam_true=np.linalg.eigvalsh((Q+Q.T)/2)[0]
    print(f"a={a}:")
    print(f"   min diagonal residual (a_i+B_ii) = {diag.min():+.4f}   (all>0 => mode-by-mode OK)")
    print(f"   typical off-diagonal row sum      = {np.median(offrow):.3f}  (max {offrow.max():.3f})")
    print(f"   min Gershgorin bound              = {gersh.min():+.3f}   ({'POSITIVITY PROVEN locally' if gersh.min()>=0 else 'FAILS -> global cancellation, not diagonally dominant'})")
    print(f"   true λ_min(Q)                     = {lam_true:+.3e}")
    print(f"   => off-diagonal / diagonal-residual ratio = {np.median(offrow)/max(diag.min(),1e-9):.0f}x  (>>1 means global)")
