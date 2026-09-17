import numpy as np
from numpy.polynomial.hermite import hermval
# CCM Lemma 7.2: prolate eigenfunctions of PW_λ = -d/dx((λ²-x²)d/dx) + (2πλx)² on [-λ,λ]
# converge to Hermite (harmonic-oscillator) eigenfunctions with O(λ^{-2}) error as λ→∞.
# Near x=0: PW_λ ≈ λ²[-∂² + (2π)²x²], eigenfns → φ_n(x)=H_n(√(2π)x)e^{-πx²}.

def hermite_fn(n, x):
    # normalized harmonic-oscillator eigenfunction for -∂²+(2π)²x² (freq ω=2π): φ_n = H_n(√ω x) e^{-ω x²/2}
    w = 2*np.pi
    c = np.zeros(n+1); c[n] = 1.0
    Hn = hermval(np.sqrt(w)*x, c)
    phi = Hn*np.exp(-0.5*w*x*x)
    nrm = np.sqrt(np.trapezoid(phi*phi, x))
    return phi/ (nrm if nrm>0 else 1)

def prolate_eigs(lam, M=2000, nmodes=6):
    x = np.linspace(-lam, lam, M+1); h = x[1]-x[0]
    # FD for -d/dx((λ²-x²) d/dx): use midpoint coefficients p(x)=λ²-x²
    xm = 0.5*(x[:-1]+x[1:]); p = lam*lam - xm*xm    # p at midpoints, length M
    N = M+1
    A = np.zeros((N,N))
    for i in range(1,N-1):
        A[i,i-1] = -p[i-1]/h**2
        A[i,i+1] = -p[i]/h**2
        A[i,i]   = (p[i-1]+p[i])/h**2 + (2*np.pi*lam*x[i])**2
    # Dirichlet ends
    A[0,0]=A[-1,-1]=1e12
    ev, evec = np.linalg.eigh(A)
    # take lowest nmodes interior eigenvectors, normalize in L²
    out=[]
    for k in range(nmodes):
        v = evec[:,k].copy()
        nrm = np.sqrt(np.trapezoid(v*v, x))
        out.append((ev[k], v/(nrm if nrm>0 else 1), x))
    return out

print("CCM Lemma 7.2 check: L2 error ‖ψ_n(λ) − φ_n(Hermite)‖ vs λ  (expect O(λ^{-2}))")
print(f"{'lam':>5} " + " ".join(f"err_n{n}".rjust(9) for n in range(4)))
rows=[]
for lam in [3,4,6,8,12]:
    eigs = prolate_eigs(lam, M=2400, nmodes=4)
    errs=[]
    for n,(mu,v,x) in enumerate(eigs):
        phi = hermite_fn(n, x)
        # align sign
        if np.trapezoid(v*phi,x) < 0: v = -v
        err = np.sqrt(np.trapezoid((v-phi)**2, x))
        errs.append(err)
    rows.append((lam,errs))
    print(f"{lam:>5} " + " ".join(f"{e:9.3e}" for e in errs), flush=True)

# fit rate for n=0: err ~ lam^p
import numpy as np
lams=np.array([r[0] for r in rows]); e0=np.array([r[1][0] for r in rows])
good=e0>1e-9
if good.sum()>=3:
    p=np.polyfit(np.log(lams[good]), np.log(e0[good]),1)[0]
    print(f"\nn=0 error scaling: err ~ lam^{p:.2f}   (CCM Lemma 7.2 predicts -2)")
