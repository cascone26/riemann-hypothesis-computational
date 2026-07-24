#!/usr/bin/env python3
"""
subgoalH_operator_structure.py — LEDGER §144 (Lens 1→5): is the bottom of the compact Weil form GENUINELY an effective
local Schrödinger operator −d²/dξ²+V(ξ), or only heuristically k²-like? §143 found λ_k≈C·k² with prime-distorted sine
modes. A local 2nd-order operator −d²/dξ²+V, in the sine basis {v_k=sin(kπ(x+a)/2a)}, looks like: DIAGONAL ∝ k² (the
−d²) + a TOEPLITZ/HANKEL band from V (multiplication operator ⟹ matrix elements = sine-transform of V, decaying
off-diagonally iff V is smooth). So build the Weil form in the sine basis and TEST:
  (1) diagonal M_kk vs k — does it grow ∝ k² (with what coefficient)?
  (2) off-diagonal decay — is M banded (LOCAL ⟹ real differential operator) or dense (NONLOCAL ⟹ Schrödinger picture
      is only heuristic)? Quantify band-energy fraction.
  (3) if banded, reconstruct the effective potential V from the off-diagonal (Toeplitz) part and print its shape —
      does V carry visible prime structure?
Uses the sine-Galerkin Weil matrix Q_kl=∫∫ g(x−y) v_k'(x) v_l'(y); operator M = Q/a (mass ‖v‖²=aΣc²). Eigenvalues of M
are the λ_k (cross-check vs §143). FULL vs ARCH-only to isolate the prime part of the potential.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

def weil_sine_matrix(a, gfun, K=30, M=3000):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h; gv=np.array([gfun(d,Lam,Nmax) for d in diffs]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)                         # v_k'(x)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h; Wp=Vp*w[:,None]
    Q=Wp.T@G@Wp; Q=(Q+Q.T)/2
    return Q/a                                             # operator M in the sine basis (eigenvalues = λ_k)

def analyze(M, tag):
    K=M.shape[0]; k=np.arange(1,K+1)
    diag=np.diag(M)
    # fit diagonal ∝ k^p
    pos=diag>0
    p,logc=np.polyfit(np.log(k[pos][:20]), np.log(diag[pos][:20]),1)
    # band-energy: fraction of Frobenius energy within |k-l|<=b
    fro=np.sqrt(np.sum(M**2))
    fracs={}
    for b in [0,1,2,4,8]:
        mask=np.abs(np.subtract.outer(k,k))<=b
        fracs[b]=np.sqrt(np.sum((M*mask)**2))/fro
    ev=np.linalg.eigvalsh(M)
    print(f"  [{tag}] diag∝k^{p:.2f} (coef {np.exp(logc):.2e}); λ_min={ev[0]:.3e}; "
          f"band-frac |k-l|≤{{0,1,2,4,8}} = {[round(fracs[b],3) for b in [0,1,2,4,8]]}")
    return diag, ev

print("Lens 1→5 §144: is the compact-form bottom a LOCAL Schrödinger operator (banded in sine basis) or nonlocal?\n")
for a in [2.0, 2.5]:
    print(f"=== a={a} ===")
    MF=weil_sine_matrix(a, screw_g); MA=weil_sine_matrix(a, screw_g_no_primes)
    dF,evF=analyze(MF,"FULL"); dA,evA=analyze(MA,"ARCH")
    # the PRIME part of the operator = FULL - ARCH
    MP=MF-MA; dP,_=analyze(MP,"PRIME=FULL-ARCH")
    # effective potential from the Toeplitz (k-l) part of the PRIME operator: V̂(d)=avg of MP over |k-l|=d
    K=MP.shape[0]; Td=[]
    for d in range(0,12):
        vals=[MP[i,i+d] for i in range(K-d)]; Td.append(np.mean(vals))
    print(f"    PRIME Toeplitz band V̂(d), d=0..11: {np.round(Td,4)}")
    # reconstruct V(ξ) on [0,π] from the cosine coefficients Td (real, even) — look for structure
    xi=np.linspace(0,np.pi,200); V=Td[0]+2*sum(Td[d]*np.cos(d*xi) for d in range(1,len(Td)))
    print(f"    reconstructed prime-potential V(ξ): min={V.min():.3f} max={V.max():.3f} mean={np.mean(V):.3f} "
          f"@argmin ξ/π={xi[np.argmin(V)]/np.pi:.2f}")
    print()
print("READ: band-frac →1 quickly (e.g. ≥0.9 by |k-l|≤4) ⟹ the operator is LOCAL = a genuine −d²/dξ²+V Schrödinger op,")
print("and V=FULL−ARCH Toeplitz band IS the prime potential. band-frac staying low (dense matrix) ⟹ nonlocal, the")
print("Schrödinger picture is only heuristic and the real object is an integral operator (still k²-laddered, but not local).")
