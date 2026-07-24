# Connes Step 2 analytic skeleton (LEDGER §98): the near-null Weil cluster = prolate ANGLE spectrum. Compute the
# Slepian concentration eigenvalues nu_m = chi_m^2 (fraction of Fourier energy of h_{2m,lambda} inside [-lambda,lambda]),
# verify 1-chi_2 vs Connes' EXPLICIT formula (line 1149), and get the chi-gap that must replace the vacuous Weil gap.
# Connes: F(h_{2m,lambda})|_{[-l,l]} = chi_m h_{2m,lambda}, chi_m^2 = nu_m; h_0<->chi_0, h_4<->chi_2. 1-chi_2 is the
# doubly-exponentially small quantity that controls eps(lambda) (Weil ground eigenvalue). Verify: 1-chi_2 ~ (2^14/3)
# sqrt2 pi^5 e^{-4 pi e^L + 9/2 L}, L=2 log lambda.
import numpy as np
from numpy.polynomial.legendre import legval

def y2_matrix(Nleg):
    Y=np.zeros((Nleg+2,Nleg+2))
    for n in range(Nleg+1):
        if n+1<=Nleg+1: Y[n+1,n]=(n+1)/(2*n+1)*np.sqrt((2*n+1)/(2*n+3))
        if n-1>=0:      Y[n-1,n]=n/(2*n+1)*np.sqrt((2*n+1)/(2*n-1))
    return (Y[:Nleg,:Nleg])@(Y[:Nleg,:Nleg])

def prolate_even(lmbda,Nleg):
    c=2*np.pi*lmbda**2; n=np.arange(Nleg)
    A=np.diag((n*(n+1)).astype(float))+c**2*y2_matrix(Nleg); A=(A+A.T)/2
    w,V=np.linalg.eigh(A)
    par=np.array([(-1)**int(np.argmax(np.abs(V[:,k]))) for k in range(V.shape[1])])
    ei=np.where(par>0)[0]
    return V,ei,Nleg

def concentration(V,col,lmbda,Nleg,Nx=6000):
    """nu = fraction of Fourier energy of h (supported [-l,l]) landing back inside [-l,l]."""
    x=np.linspace(-lmbda,lmbda,Nx); dx=x[1]-x[0]
    coeff=np.array([V[m,col]*np.sqrt((2*m+1)/2) for m in range(Nleg)])
    h=legval(x/lmbda,coeff)
    h=h/np.sqrt(np.sum(h*h)*dx)                       # normalize ∫|h|²=1
    # FT hhat(y)=∫ h(x) e^{2πixy} dx, evaluate on y∈[-l,l] (need concentration there)
    y=np.linspace(-lmbda,lmbda,Nx); dy=y[1]-y[0]
    # hhat(y_k) = sum_x h(x) e^{2πi x y_k} dx  -> matrix (use real/imag)
    ph=2*np.pi*np.outer(y,x)
    hhat=(h[None,:]*np.exp(1j*ph)).sum(1)*dx
    nu=np.sum(np.abs(hhat)**2)*dy                     # energy inside [-l,l] (total energy =1 by Parseval)
    return nu

print("PROLATE ANGLE SPECTRum chi_m (Slepian concentration) — the fine structure Connes uses (LEDGER §98).")
print(f"{'lambda':>7} {'L=2ln l':>8} {'chi_0':>10} {'1-chi_0':>12} {'chi_2':>10} {'1-chi_2':>12} {'Connes 1-chi_2':>15}")
for lmbda in [1.4,1.6,1.8,2.0,2.3]:
    Nleg=max(120,int(5*2*np.pi*lmbda**2))
    Vm,ei,Nl=prolate_even(lmbda,Nleg)
    i0,i4=ei[0],ei[2]
    nu0=concentration(Vm,i0,lmbda,Nl); nu2=concentration(Vm,i4,lmbda,Nl)
    chi0,chi2=np.sqrt(min(nu0,1)),np.sqrt(min(nu2,1))
    L=2*np.log(lmbda)
    connes=(2**14/3)*np.sqrt(2)*np.pi**5*np.exp(-4*np.pi*np.exp(L)+4.5*L)
    print(f"{lmbda:>7.2f} {L:>8.3f} {chi0:>10.6f} {1-chi0:>12.3e} {chi2:>10.6f} {1-chi2:>12.3e} {connes:>15.3e}",flush=True)
print("\nREAD: does numeric 1-chi_2 track Connes' explicit doubly-exp formula? And chi-gap (1-chi_0 vs 1-chi_2)")
print("is the fine structure separating h_0 from h_4 in the near-null cluster (the object for the Step-2 bound).")
