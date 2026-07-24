# LENS 2 (brainstorm, 2026-07-21): does the PRIME term break the archimedean near-null cluster degeneracy and SELECT
# an EVEN ground state ≈ k_λ?  Degenerate perturbation theory: A_λ = Q_arch + Q_prime. The near-null cluster of
# Q_arch (the archimedean prolate/Sonin skeleton, Connes-Moscovici) is near-degenerate; Q_prime lifts it. KEY: the
# cluster-restricted Q_prime has O(1) entries (float64-IMMUNE) even though the absolute eigenvalues are ~1e-14 — so
# the SELECTION/ordering is well-conditioned even where the level is not (the §84 curse is on the level, not the order).
# Test: build Q_full (screw_g) & Q_arch (screw_g_no_primes) in the sine basis, take Q_arch's lowest-d cluster, project
# Q_prime=Q_full-Q_arch onto it (O(1)?), diagonalize Q_arch+Q_prime IN the cluster, check the ground state: EVEN? ≈k_λ?
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto
sys.path.insert(0,'/Users/scones/projects/riemann/operators/dbn_dynamics')
from stage3_theta_vs_k import k_lambda_on

def sine_weil(a, gfun, Lam, Nmax, M=8000, K=120):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gv=np.array([gfun(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg); Vv=np.sin(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Q=(Wp.T@G@Wp); Q=(Q+Q.T)/2
    return Q/a, x, Vv, w

for lmbda in [1.8, 2.2]:
    a=np.log(lmbda); Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    Hf,x,Vv,w=sine_weil(a,screw_g,Lam,Nmax)
    Ha,_,_,_ =sine_weil(a,screw_g_no_primes,Lam,Nmax)
    Hp=Hf-Ha                                     # prime part
    # archimedean near-null cluster: lowest d eigenvectors of Ha
    ea,Ua=np.linalg.eigh(Ha)
    d=6
    U=Ua[:,:d]                                   # cluster basis (arch near-null)
    Lam_arch=U.T@Ha@U                            # arch levels in cluster (~diag, tiny)
    P=U.T@Hp@U                                   # PRIME perturbation in cluster — O(1)?
    P=(P+P.T)/2
    evc,ecc=np.linalg.eigh(Lam_arch+P)           # full A_λ ground state within the cluster
    cmin=U@ecc[:,0]                              # back to sine coeffs
    vmin=Vv@cmin                                 # ground state on grid
    # parity
    ov_sym=np.trapezoid(vmin*vmin[::-1]*w)/np.trapezoid(vmin*vmin*w)
    parity='EVEN' if ov_sym>0.9 else ('ODD' if ov_sym<-0.9 else 'mixed')
    # compare to k_λ (project to same grid, unit-align)
    k=k_lambda_on(x,lmbda)
    def unit(z): return z/np.sqrt(np.trapezoid(z*z*w))
    vk,kk_=unit(vmin),unit(k)
    ov_k=abs(np.trapezoid(vk*kk_*w))
    print(f"===== lambda={lmbda} (a={a:.3f}) =====")
    print(f" arch cluster levels (lowest {d}): {np.array2string(ea[:d],formatter={'float':lambda z:f'{z:.2e}'})}")
    print(f" PRIME perturbation P in cluster: entry scale |P|max={np.max(np.abs(P)):.3e}, diag={np.array2string(np.diag(P),formatter={'float':lambda z:f'{z:+.2e}'})}")
    print(f" ground state of (arch+prime) in cluster: parity={parity} (sym={ov_sym:+.3f}), overlap with k_λ = {ov_k:.4f}")
    # does prime alone reorder? show eigenvector composition in arch-cluster basis
    comp=ecc[:,0]/np.max(np.abs(ecc[:,0]))
    print(f" ground-state composition in arch-cluster modes: {np.array2string(comp,formatter={'float':lambda z:f'{z:+.2f}'})}")
print("\nREAD: (1) is |P| O(1) while arch levels ~1e-7? (float64-immune selection). (2) does the prime perturbation")
print("give an EVEN ground state? (3) overlap with k_λ -> is k_λ the prime-selected direction? YES on all = Lens 2.")
