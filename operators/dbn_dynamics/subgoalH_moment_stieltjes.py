#!/usr/bin/env python3
"""
subgoalH_moment_stieltjes.py — LEDGER §148 (Lens 1/A, orthogonal-polynomial route): the compact Weil form's bottom is a
JACOBI operator J (§147). Classical OP theory: J ≥ 0 (⟺ RH here) ⟺ its orthogonality measure is supported on [0,∞) ⟺
the MOMENT sequence m_k = ⟨e_1, J^k e_1⟩ = ∫ λ^k dμ(λ) is a STIELTJES sequence (both Hankel matrices [m_{i+j}] and the
shifted [m_{i+j+1}] positive-semidefinite). The moments are explicit — computable from J directly. NEW HANDLE: recast
"prove λ_min>0" as "prove the Weil moment sequence is Stieltjes", and look at HOW the primes shift the moments (FULL vs
ARCH): a sign-DEFINITE prime shift would be a structural positivity mechanism.
This tool: build J (derivative-Galerkin Weil operator, §147), compute moments, test Stieltjes (Hankel margins), and
print the arch→full moment shift to see if the primes move the sequence in one direction.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

def weil_operator(a, gfun, K=24, M=2400):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([gfun(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    kk=np.arange(1,K+1); Vp=(kk*np.pi/(2*a))*np.cos(np.pi*np.outer((x+a),kk)/(2*a))
    Wp=Vp*w[:,None]; Q=Wp.T@G@Wp; return (Q+Q.T)/2/a       # operator, eigenvalues = λ_k>0

def moments_and_stieltjes(C, nmom=10):
    # spectral measure of C w.r.t. its GROUND eigenvector direction is not e_1; use the actual measure of C: moments
    # m_k = mean over the operator = tr via eigen: m_k = Σ_j w_j λ_j^k with weights from a cyclic vector. Use the
    # normalized all-ones-ish cyclic vector b (must be cyclic); take b = lowest few modes mix. Cleanest: use e_1 in the
    # eigenbasis-independent way — moments of the measure dμ(λ)=Σ|<b,ψ_j>|² δ(λ_j). Pick b = first canonical vector.
    ev,U=np.linalg.eigh(C); b=np.zeros(C.shape[0]); b[0]=1.0
    wj=(U.T@b)**2                                            # weights of the spectral measure w.r.t e_1
    m=np.array([np.sum(wj*ev**k) for k in range(2*nmom+1)])
    H0=np.array([[m[i+j] for j in range(nmom)] for i in range(nmom)])
    H1=np.array([[m[i+j+1] for j in range(nmom)] for i in range(nmom)])
    e0=np.linalg.eigvalsh(H0); e1=np.linalg.eigvalsh(H1)
    return ev, m, e0[0], e1[0]

print("Sub-goal H §148: Weil moment sequence & Stieltjes (moment-problem) positivity — how do primes shift the moments?\n")
for a in [1.5, 2.0, 2.5]:
    CF=weil_operator(a, screw_g); CA=weil_operator(a, screw_g_no_primes)
    evF,mF,h0F,h1F = moments_and_stieltjes(CF)
    evA,mA,h0A,h1A = moments_and_stieltjes(CA)
    print(f"=== a={a} ===")
    print(f"  FULL λ_min={evF[0]:+.3e}   Hankel-min[m_ij]={h0F:+.3e}  shifted Hankel-min[m_i+j+1]={h1F:+.3e}  "
          f"{'STIELTJES✓' if h0F>-1e-14 and h1F>-1e-14 else 'not-Stieltjes'}")
    print(f"  ARCH λ_min={evA[0]:+.3e}   Hankel-min={h0A:+.3e}  shifted={h1A:+.3e}")
    # moment shift primes cause: full - arch, first few moments
    shift=mF[:6]-mA[:6]
    print(f"  moments FULL m_0..5 = {np.array2string(mF[:6], formatter={'float':lambda z:f'{z:.3e}'})}")
    print(f"  prime SHIFT (FULL-ARCH) m_0..5 = {np.array2string(shift, formatter={'float':lambda z:f'{z:+.3e}'})}  "
          f"signs={''.join('+' if s>0 else '-' for s in shift)}")
    print()
print("READ: both Hankel mins ≥0 ⟹ Stieltjes ⟹ measure on [0,∞) ⟹ positivity (consistency check). KEY: if the prime")
print("moment SHIFT has a consistent SIGN pattern across a & k, the primes move the moment sequence in a definite")
print("direction — a candidate structural mechanism (does the prime shift PRESERVE the Stieltjes property?).")
