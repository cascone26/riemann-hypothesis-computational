"""Step(a) 'even' part: does the EVEN sector always contain the min eigenvalue of the Weil form? (persymmetric
decouples even/odd; Connes needs the ground state EVEN). Compute even-sector vs odd-sector min across a."""
import numpy as np
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def sectors(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gv[(idx[:,None]-idx[None,:])+M]
    w=np.sqrt(np.full(M+1,h)); w[0]=w[-1]=np.sqrt(0.5*h)
    Gs=w[:,None]*G*w[None,:]; Gs=(Gs+Gs.T)/2
    n=M+1
    # even/odd projectors via the flip J
    J=np.eye(n)[::-1]
    Pe=(np.eye(n)+J)/2; Po=(np.eye(n)-J)/2
    # restrict to even/odd subspaces: eigenvalues of Pe Gs Pe on even range
    # build even/odd basis
    ev_min=np.min(np.linalg.eigvalsh(Pe@Gs@Pe + 1e6*Po))   # penalize odd
    od_min=np.min(np.linalg.eigvalsh(Po@Gs@Po + 1e6*Pe))
    return ev_min, od_min

print("Step(a) even part: even-sector min vs odd-sector min (Weil form). Ground state even iff even<odd.")
print(f"{'a':>5} {'even min':>12} {'odd min':>12} {'even<odd?':>10}")
for a in [0.5,0.8,1.2,1.8,2.5,3.5]:
    M=min(900,max(400,int(350/a)))
    em,om=sectors(a,M)
    print(f"{a:5.2f} {em:12.5e} {om:12.5e} {str(em<om):>10}")
print("\nIf even<odd for ALL a: the Weil ground state is ALWAYS EVEN → step(a)'s 'even' part holds uniformly.")
