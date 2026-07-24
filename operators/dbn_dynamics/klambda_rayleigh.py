# Connes Step 2 core (LEDGER §99 NEXT): is k_lambda a NEAR-NULL vector of the Weil form? Compute the Rayleigh
# quotient QW_lambda(k_lambda)/||k_lambda||^2 DIRECTLY on the known k_lambda (sidesteps the ill-conditioned eigenvector
# of S97). Connes line 1161-1165: k_lambda is the prolate near-radical vector, so QW_lambda should be doubly-exp small
# and track eps(lambda) ~ 1-chi_2. If so, k_lambda ACHIEVES the near-null value = the core of Step 2.
# Convention (S97): the match is Suzuki's v itself (even), and the Weil form on v is Q(v)=int int g(x-y) v'(x)v'(y).
# So we evaluate Q(k_lambda) with k_lambda in the ROLE of v, i.e. using k_lambda'(x).
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from stage3_theta_vs_k import k_lambda_on

def weil_rayleigh_of_k(a, lmbda, M=8000):
    """Q(k)/||k||^2 with k=k_lambda(e^s) on s in (-a,a), Q(v)=int int g(s-s') v'(s)v'(s') ds ds'."""
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    s=np.linspace(-a,a,M+1); h=s[1]-s[0]
    diffs=np.arange(-M,M+1)*h
    gvals=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1); G=gvals[(idx[:,None]-idx[None,:])+M]
    k=k_lambda_on(s,lmbda)
    # enforce Dirichlet-ish: k should ~vanish at +/-a (support edge); subtract endpoints linearly to be safe for v'
    kp=np.gradient(k,s)                          # k'(s)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wk=kp*w
    Q=float(Wk@G@Wk)                             # int int g(s-s') k'(s) k'(s')
    norm2=float(np.sum(k*k*w))                   # ||k||^2 = int |k|^2
    return Q/norm2, Q, norm2

# reference: the smallest Weil eigenvalue lam_a (floored ~1e-7) and Connes 1-chi_2 (the true scale)
def connes_1mchi2(lmbda):
    L=2*np.log(lmbda); return (2**14/3)*np.sqrt(2)*np.pi**5*np.exp(-4*np.pi*np.exp(L)+4.5*L)

print("Rayleigh quotient of the prolate guess k_lambda in the Weil form  (Connes Step-2 core):")
print(f"{'lambda':>7} {'a':>6} {'QW(k)/||k||^2':>15} {'Connes 1-chi_2':>15} {'ratio':>9}")
for lmbda in [1.4,1.6,1.8,2.0,2.2]:
    a=np.log(lmbda)
    R,Q,n2=weil_rayleigh_of_k(a,lmbda)
    c=connes_1mchi2(lmbda)
    print(f"{lmbda:>7.2f} {a:>6.3f} {R:>15.4e} {c:>15.4e} {R/c if c>0 else float('nan'):>9.3f}",flush=True)
print("\nREAD: is QW(k)/||k||^2 doubly-exp small (following 1-chi_2 down)? If yes, k_lambda ACHIEVES the near-null")
print("value = k_lambda is a genuine near-minimizer of the Weil form = the core of Connes Step 2 (k~theta_x).")
print("(A ~constant ratio to Connes' formula across lambda = the doubly-exp RATE captured; constant = normalization.)")
