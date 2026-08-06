"""Reproduce Connes-Moscovici prolate operator W_λ = -d/dx((λ²-x²)d/dx) + (2π)²x² (arXiv:2112.05500).
Its NEGATIVE eigenvalues should reproduce -(γ²) (squares of zeta zeros) in the UV. Discretize by
finite differences (symmetric Sturm-Liouville) and compare the negative spectrum to -γ² for the
actual first zeros. Probe: is the match EXACT or only UV-asymptotic (the paper claims asymptotic)?"""
import numpy as np
zeros = np.array([14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,
                  43.327073,48.005151,49.773832,52.970321,56.446248,59.347044,60.831779])
def prolate_eigs(lam, X, Ngrid):
    x = np.linspace(-X, X, Ngrid); h = x[1]-x[0]
    xm = (x[:-1]+x[1:])/2                     # midpoints
    p_mid = lam**2 - xm**2                     # p(x)=λ²-x² at midpoints
    q = (2*np.pi)**2 * x**2                    # q(x)=(2π)²x²
    # -(p ξ')' via symmetric FD: tri-diagonal
    main = np.zeros(Ngrid); off = np.zeros(Ngrid-1)
    for i in range(Ngrid):
        pl = p_mid[i-1] if i>0 else 0.0
        pr = p_mid[i]   if i<Ngrid-1 else 0.0
        main[i] = (pl+pr)/h**2 + q[i]
    for i in range(Ngrid-1):
        off[i] = -p_mid[i]/h**2
    from scipy.linalg import eigh_tridiagonal
    ev = eigh_tridiagonal(main, off, eigvals_only=True)
    return np.sort(ev)
for lam in [1.0, np.sqrt(2)]:
    ev = prolate_eigs(lam, 6.0, 6000)
    neg = ev[ev<-1][:14]   # most-negative eigenvalues
    print(f"\nλ={lam:.4f}: most-negative eigenvalues vs -γ² (first zeros):")
    print(f"  {'n':>2} {'E_n(neg)':>12} {'-γ_n²':>12} {'sqrt(-E_n)':>11} {'γ_n':>10} {'ratio':>8}")
    negsorted = np.sort(neg)  # ascending (most negative first)
    for n in range(min(10,len(zeros))):
        En = negsorted[n] if n < len(negsorted) else float('nan')
        g2 = -zeros[n]**2
        sq = np.sqrt(-En) if En<0 else float('nan')
        print(f"  {n+1:>2} {En:>12.2f} {g2:>12.2f} {sq:>11.3f} {zeros[n]:>10.3f} {sq/zeros[n] if En<0 else 0:>8.4f}")
