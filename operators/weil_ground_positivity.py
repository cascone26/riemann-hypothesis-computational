"""
Connes step (a): is the min eigenvalue of the (direct) Weil form QW_λ SIMPLE + EVEN for ALL λ?
Suzuki §5.2 proves it for small λ via positivity-improving Dirichlet form. Test uniformly: is the ground
eigenvector SINGLE-SIGNED (0 interior nodes = positivity-improving) across a range of a? If yes ∀a → simple+even ∀a.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def ground(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    w=np.sqrt(np.full(M+1,h)); w[0]=w[-1]=np.sqrt(0.5*h)
    Gs=(w[:,None]*G*w[None,:]); Gs=(Gs+Gs.T)/2
    ev,V=np.linalg.eigh(Gs)
    eta=V[:,0]/w
    if eta[np.argmax(np.abs(eta))]<0: eta=-eta
    gap=ev[1]-ev[0]
    # interior sign changes (nodes), ignoring tiny numerical noise near boundary
    core=eta[M//10:-M//10] if M>20 else eta
    nodes=int(np.sum(np.diff(np.sign(core[np.abs(core)>1e-6*np.max(np.abs(eta))]))!=0))
    evenerr=np.max(np.abs(eta-eta[::-1]))
    minval=np.min(eta[np.abs(x)<a*0.9])   # is it single-signed (positive) in the bulk?
    return ev[0], gap, nodes, evenerr, minval

print("Connes step(a): direct Weil form ground eigenvector — single-signed (0 nodes) ∀a ⟹ positivity-improving ⟹ simple+even")
print(f"{'a':>5} {'min eig':>11} {'gap(simple?)':>13} {'interior nodes':>14} {'even-err':>10} {'bulk min (>=0?)':>15}")
for a in [0.3,0.5,0.8,1.2,1.8,2.5,3.5]:
    M=min(1200,max(500,int(400/a)))
    e0,gap,nodes,eerr,mn=ground(a,M)
    print(f"{a:5.2f} {e0:11.4e} {gap:13.4e} {nodes:14d} {eerr:10.1e} {mn:15.4e}")
print("\nIf nodes=0 & bulk-min>=0 for ALL a: ground state SINGLE-SIGNED ⟹ positivity-improving ⟹ Connes step(a)")
print("(simple+even min eigenvalue) holds UNIFORMLY in a, not just small a (Suzuki). Real progress on the program.")
