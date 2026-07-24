#!/usr/bin/env python3
"""
prolate_weil_v1.py — LEDGER §163 (tool rebuild, increment 1): the sine-Galerkin Weil operator SATURATES (~1e-5 at a≈1.5,
K~40) while the CvS/Groskin prolate-adapted operator reaches ~1e-59 (c=13). Hypothesis: the truncated Weil form is a
BAND-LIMITING operator, so PROLATE spheroidal wave functions (PSWF) — its natural eigenbasis — concentrate the near-null
into the lowest modes and reach FAR deeper than sines. TEST (float64 first, before hi-prec): build Q_{kl}=∫∫ g(x-y)
φ_k'(x)φ_l'(y) in the PSWF basis (bandwidth-scanned) vs the sine basis, same K. Does prolate λ_min drop BELOW the sine
saturation? If yes → prolate is the right basis (proceed to hi-prec/CvS-grade). If no → reconsider.
PSWF from prolate_operator.py (Legendre-Galerkin eigenfns of −d((1−y²)d)+c_band² y² on [−1,1]); map y→x=a y.
"""
import numpy as np, sys
sys.path.insert(0,'/Users/scones/projects/riemann/operators')
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto
from numpy.polynomial import legendre as LEG

def pswf_basis(cband, K, Nleg, ygrid):
    """lowest K PSWF (as functions on ygrid∈[-1,1]) and their y-derivatives, via Legendre-Galerkin."""
    n=np.arange(Nleg)
    # y^2 in orthonormal Legendre basis
    Y=np.zeros((Nleg+2,Nleg+2))
    for m in range(Nleg+1):
        if m+1<=Nleg+1: Y[m+1,m]=(m+1)/(2*m+1)*np.sqrt((2*m+1)/(2*m+3))
        if m-1>=0:      Y[m-1,m]=m/(2*m+1)*np.sqrt((2*m+1)/(2*m-1))
    Y=Y[:Nleg,:Nleg]; Y2=Y@Y
    Aop=np.diag((n*(n+1)).astype(float))+cband**2*Y2; Aop=(Aop+Aop.T)/2
    w,V=np.linalg.eigh(Aop)                      # columns = PSWF in orthonormal-Legendre coeffs
    # evaluate lowest K PSWF + derivatives on ygrid
    Pn=np.array([np.sqrt((2*k+1)/2)*LEG.legval(ygrid, [0]*k+[1]) for k in range(Nleg)])   # (Nleg, ngrid) orthonormal P~_k
    dPn=np.array([np.sqrt((2*k+1)/2)*LEG.legval(ygrid, LEG.legder([0]*k+[1])) for k in range(Nleg)])
    phis=(V[:,:K].T@Pn)      # (K, ngrid)
    dphis=(V[:,:K].T@dPn)    # (K, ngrid)  d/dy
    return phis, dphis

def weil_lammin(a, basis, K=40, M=1600, cband=None):
    Nmax=int(np.floor(np.exp(2*a)))+1; Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    dif=np.arange(-M,M+1)*h; gv=np.array([screw_g(d,Lam,Nmax) for d in dif]); idx=np.arange(M+1)
    G=gv[(idx[:,None]-idx[None,:])+M]
    if basis=='sine':
        kk=np.arange(1,K+1); Vp=(kk*np.pi/(2*a))*np.cos(np.pi*np.outer((x+a),kk)/(2*a))  # (M+1,K)
        Vv=np.sin(np.pi*np.outer((x+a),kk)/(2*a))
        Vp=Vp.T; Vv=Vv.T
    else:  # prolate
        yg=x/a; phis,dphis=pswf_basis(cband,K,max(K+40,120),yg)
        Vv=phis; Vp=dphis/a                                        # φ_k(x), φ_k'(x)=(1/a)φ_k'(y)
    Wp=(Vp*w).T                                                    # (M+1,K)
    Q=Wp.T@G@Wp; Q=(Q+Q.T)/2
    Mass=(Vv*w)@Vv.T; Mass=(Mass+Mass.T)/2
    ev=np.linalg.eigvalsh(np.linalg.solve(Mass,Q)) if False else None
    # symmetric generalized: Mass^{-1/2} Q Mass^{-1/2}
    eM,UM=np.linalg.eigh(Mass); eM=np.clip(eM,1e-14*eM.max(),None); Mih=UM@np.diag(eM**-0.5)@UM.T
    C=Mih@Q@Mih; C=(C+C.T)/2
    return np.linalg.eigvalsh(C)[0]

a=1.5; K=40
print(f"§163 prolate-vs-sine Weil λ_min (a={a}, K={K}, float64):\n")
ls=weil_lammin(a,'sine',K=K)
print(f"  SINE basis:        λ_min = {ls:.4e}   (the saturation baseline)")
print(f"  PROLATE basis (bandwidth scan):")
best=(1e9,None)
for cb in [2.0,4.0,7.0,10.0,15.0,20.0,30.0,45.0]:
    lp=weil_lammin(a,'prolate',K=K,cband=cb)
    tag='  <-- below sine' if lp<ls*0.5 else ''
    if lp<best[0]: best=(lp,cb)
    print(f"     c_band={cb:>5.1f}:  λ_min = {lp:>.4e}{tag}",flush=True)
print(f"\n  best prolate: λ_min={best[0]:.4e} at c_band={best[1]}   vs sine {ls:.4e}   ratio={best[0]/ls:.2e}")
print("\nREAD: if best prolate λ_min ≪ sine λ_min ⟹ PSWF is the better-conditioned basis (reaches deeper toward the")
print("doubly-exp near-null) — the tool rebuild is on the right track; next = high precision (mpmath/Olympus) + tune")
print("bandwidth + verify vs Groskin's λ_min(c). If prolate ≈ sine ⟹ conditioning isn't the basis; reconsider.")
