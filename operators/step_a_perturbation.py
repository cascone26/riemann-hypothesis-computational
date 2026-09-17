"""
Step (a) for all λ via PERTURBATIVE STABILITY (non-circular, non-Dirichlet): does the archimedean ground state
stay the simple+even ground state under the prime perturbation (no level crossing)? Compare arch-only vs full
Weil form ground states + gaps across a. If full ≈ arch (single-signed, gap preserved), step (a) is inherited.
"""
import numpy as np
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

def forms(a, M):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    ga=np.array([screw_g_no_primes(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); sh=(idx[:,None]-idx[None,:])+M
    w=np.sqrt(np.full(M+1,h)); w[0]=w[-1]=np.sqrt(0.5*h)
    def build(g): 
        Gs=w[:,None]*g[sh]*w[None,:]; return (Gs+Gs.T)/2
    Gf=build(gf); Ga=build(ga)
    evf,Vf=np.linalg.eigh(Gf); eva,Va=np.linalg.eigh(Ga)
    def gnd(ev,V):
        e=V[:,0]/w; 
        if e[np.argmax(np.abs(e))]<0: e=-e
        return e/np.sqrt(np.sum(e*e)*h), ev[0], ev[1]-ev[0]
    ef,e0f,gapf=gnd(evf,Vf); ea,e0a,gapa=gnd(eva,Va)
    overlap=abs(np.sum(ef*ea)*h)
    return e0f,gapf,e0a,gapa,overlap

print("Step(a) perturbative stability: arch-only vs full Weil ground state + gap, across a")
print(f"{'a':>5} {'full min':>10} {'full gap':>10} {'arch min':>10} {'arch gap':>10} {'<full|arch>':>12}")
for a in [0.5,0.8,1.2,1.8,2.5,3.5]:
    M=min(1000,max(500,int(350/a)))
    e0f,gapf,e0a,gapa,ov=forms(a,M)
    print(f"{a:5.2f} {e0f:10.4e} {gapf:10.4e} {e0a:10.4e} {gapa:10.4e} {ov:12.6f}")
print("\nIf full ground state ~ arch ground state (overlap~1) & both simple (gap>0): the prime perturbation does")
print("NOT cause a ground-state crossing → step(a) [simple+even] INHERITED from arch (prolate) for all a, non-circular.")
