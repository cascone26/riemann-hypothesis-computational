"""NEW probe: the cancellation mechanism of Weil positivity. smooth-arch (core+tail) has κ_arch negative
squares; the full form arch-prime(+pole) has κ=0. HOW does prime+pole cancel arch's negativity?
Characterize: project the prime+pole operator onto arch's NEGATIVE eigenspace — is it rank-matched and
sign-correct to lift exactly those modes to >=0? A clean structure = the positivity mechanism."""
import numpy as np, sys
sys.path.insert(0,'operators')
from suzuki_exact_Q import build_Q

for a in [1.0,1.5,2.0]:
    K=max(70,int(30*a)+40)
    Q,P=build_Q(a,K,NG=400,return_pieces=True)
    smooth = (P["core"]+P["tail"]); smooth=(smooth+smooth.T)/2
    rest   = (P["pole"]+P["prime"]); rest=(rest+rest.T)/2
    full   = (Q+Q.T)/2
    wa,Va = np.linalg.eigh(smooth)          # arch smooth eigenmodes
    neg = wa < -1e-9
    kneg = int(neg.sum())
    Vneg = Va[:,neg]                        # negative eigenspace of smooth-arch
    # project 'rest' (prime+pole) onto arch's negative eigenspace
    Rproj = Vneg.T @ rest @ Vneg            # kneg x kneg
    rw = np.linalg.eigvalsh(Rproj)
    # how much does 'rest' add on the negative modes (diagonal)?
    arch_negvals = wa[neg]
    rest_on_neg = np.diag(Rproj)
    net = arch_negvals + rest_on_neg
    wf = np.linalg.eigvalsh(full)
    print(f"a={a}: κ(smooth-arch)={kneg}, κ(full)={int(np.sum(wf<-1e-9))}")
    print(f"   arch neg eigenvalues:      {np.round(arch_negvals,3)}")
    print(f"   rest(prime+pole) on those: {np.round(rest_on_neg,3)}")
    print(f"   net (arch+rest diag):      {np.round(net,3)}  <- all >=0 means rest exactly lifts arch's neg modes")
    print(f"   rank(rest restricted to arch-neg space) eigs: {np.round(rw,3)}")
    print(f"   ||rest|| on arch-neg space / ||rest|| total = {np.linalg.norm(Rproj)/np.linalg.norm(rest):.3f}")
    print()

print("="*70)
print("COMMUTATOR TEST: do arch(smooth) and rest(prime+pole) approximately commute?")
print("If yes -> shared eigenvectors are STRUCTURAL (a common symmetry generator), the mechanism.")
for a in [1.0,1.5,2.0]:
    K=max(70,int(30*a)+40)
    Q,P=build_Q(a,K,NG=400,return_pieces=True)
    smooth=(P["core"]+P["tail"]); smooth=(smooth+smooth.T)/2
    rest=(P["pole"]+P["prime"]); rest=(rest+rest.T)/2
    C=smooth@rest-rest@smooth
    rel=np.linalg.norm(C)/(np.linalg.norm(smooth)*np.linalg.norm(rest))
    # also on the arch-negative subspace only
    wa,Va=np.linalg.eigh(smooth); Vn=Va[:,wa<-1e-9]
    Cn=Vn.T@C@Vn
    reln=np.linalg.norm(Cn)/(np.linalg.norm(Vn.T@smooth@Vn)*np.linalg.norm(Vn.T@rest@Vn)+1e-30)
    print(f"  a={a}: ||[smooth,rest]||/(||s||||r||) = {rel:.4f} (global);  on arch-neg subspace = {reln:.4f}")
    # what does smooth commute with? test vs core alone, tail alone, prime alone
    for lab,M in [("core",P['core']),("tail",(P['tail']+P['tail'].T)/2),("prime",(P['prime']+P['prime'].T)/2),("pole",(P['pole']+P['pole'].T)/2)]:
        Cx=smooth@M-M@smooth
        print(f"      [smooth,{lab}] rel = {np.linalg.norm(Cx)/(np.linalg.norm(smooth)*np.linalg.norm(M)+1e-30):.4f}")
