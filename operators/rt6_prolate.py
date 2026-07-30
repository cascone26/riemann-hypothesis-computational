#!/usr/bin/env python3
"""§462 (BUILD the CCM prolate operator — the concrete Θ). Connes-Moscovici (arXiv:2112.05500) eq(1):
  W = -d/dx[(λ²-x²) d/dx] + (2πλ)² x²   on L²(ℝ), λ=1.
Singular points x=±1 (where λ²-x²=0, coefficient flips sign → indefinite → negative eigenvalues). The paper
(§4, Thm 2.6): the self-adjoint W_sa has NEGATIVE eigenvalues μ_k whose |μ_k|<E² count as N(E)~(E/2π)(log(E/2π)-1)
= RvM ⟹ μ_k ≈ −γ_k² (minus the SQUARES of the zeta zeros). Build W by finite differences, compute eigenvalues,
check if the negative ones ≈ −γ_k². This is the concrete Θ my whole investigation pointed to."""
import numpy as np
Z=np.load('zeros500.npy')
lam=1.0
def build_W(L,N):
    x=np.linspace(-L,L,N); h=x[1]-x[0]
    xh=x[:-1]+h/2                                  # half-grid points
    p=lam**2-xh**2                                 # (λ²-x²) at half-points (negative for |x|>1)
    # -d/dx[p dξ/dx]: symmetric tridiagonal. diag_i = (p_{i-1/2}+p_{i+1/2})/h² ; off_i = -p_{i+1/2}/h²
    main=np.zeros(N); off=np.zeros(N-1)
    for i in range(N):
        pl=p[i-1] if i>0 else 0.0
        pr=p[i] if i<N-1 else 0.0
        main[i]=(pl+pr)/h**2
    off=-p/h**2
    W=np.diag(main)+np.diag(off,1)+np.diag(off,-1)
    W+=np.diag((2*np.pi*lam)**2 * x**2)            # + q(x)
    return W,x
for L,N in [(30,3000),(40,5000)]:
    W,x=build_W(L,N)
    ev=np.linalg.eigvalsh(W)
    neg=ev[ev<0]                                   # negative eigenvalues
    neg_sorted=np.sort(neg)[::-1]                  # closest to 0 first (least negative)
    print(f"\nL={L} N={N}: {len(neg)} negative eigenvalues. Compare -μ_k^{{1/2}} (=√|μ|) to zeta zeros γ_k:")
    print(f"  {'k':>3} {'μ_k (neg eig)':>14} {'√|μ_k|':>9} {'γ_k (zeta)':>10} {'ratio':>7}")
    sq=np.sqrt(-neg_sorted)                         # √|μ_k| should ≈ γ_k
    for k in range(min(12,len(sq))):
        gz=Z[k] if k<len(Z) else np.nan
        print(f"  {k+1:>3} {neg_sorted[k]:14.3f} {sq[k]:9.3f} {gz:10.3f} {sq[k]/gz if gz==gz else np.nan:7.4f}")
print("""
reading: if √|μ_k| (from the prolate operator's NEGATIVE eigenvalues) ≈ γ_k (the zeta zeros), then I have BUILT the
CCM operator correctly and it spectrally realizes the zeros (UV). The match should improve for larger k (UV) and
larger L,N. This is the concrete Θ = (scaling)²+grading realized as the prolate operator. Boundary/extension
subtleties (the exact W_sa domain, eq 16-18) affect the LOW eigenvalues; the UV (large k) is robust.
""")
