# Connes program (LEDGER §94), Step-1 tool: build the PROLATE WAVE OPERATOR pw_lambda and verify Connes' stated
# facts, to get the spectral GAP that transfers even-simplicity to QW_lambda.
# Connes eq (wop): pw_lambda = -d_x((lambda^2 - x^2) d_x) + (2 pi lambda x)^2 on [-lambda,lambda].
# Sub x = lambda*y (y in [-1,1]):  pw = -d_y((1-y^2) d_y) + c^2 y^2,  c = 2*pi*lambda^2  (= a^2, a=sqrt(2pi)lambda).
# In the NORMALIZED Legendre basis {P~_n(y)}: Legendre op has eigenvalues n(n+1); y^2 is a banded matrix.
# PSWF = eigenvectors; eigenvalues chi come from the COMMUTING concentration op P_l F P_l (computed separately).
# This is the textbook Slepian computation -> low convention risk. Verify: (a) eigenvalues SIMPLE (Connes Fact),
# (b) ~4 lambda^2 = 2c/pi concentrated modes, (c) even ground state h_{0,l} and 3rd even h_{4,l} + the GAP.
import numpy as np
from numpy.polynomial import legendre as L

def y2_matrix(Nleg):
    """Matrix of multiplication by y^2 in the ORTHONORMAL Legendre basis P~_n = sqrt((2n+1)/2) P_n, n=0..Nleg-1."""
    # recurrence: y P_n = ((n+1)P_{n+1} + n P_{n-1})/(2n+1). Build Y (mult by y) in orthonormal basis, then Y2=Y@Y.
    Y = np.zeros((Nleg+2, Nleg+2))
    for n in range(Nleg+1):
        # in orthonormal basis, <P~_m| y |P~_n>: y P_n coefficients, rescaled
        if n+1 <= Nleg+1:
            Y[n+1,n] = (n+1)/(2*n+1) * np.sqrt((2*n+1)/(2*n+3))   # coeff of P_{n+1}, orthonormal rescale
        if n-1 >= 0:
            Y[n-1,n] = n/(2*n+1) * np.sqrt((2*n+1)/(2*n-1))
    Y = Y[:Nleg,:Nleg]
    return Y @ Y

def prolate(lmbda, Nleg=400):
    c = 2*np.pi*lmbda**2
    n = np.arange(Nleg)
    Aop = np.diag(n*(n+1).astype(float)) + c**2 * y2_matrix(Nleg)
    Aop = (Aop+Aop.T)/2
    w, V = np.linalg.eigh(Aop)                 # ascending; PSWF eigenvalues (operator, not concentration)
    return c, w, V, n

print("PROLATE WAVE OPERATOR pw_lambda (Legendre-Galerkin). Verifying Connes' stated facts.\n")
print(f"{'lambda':>7} {'c=2pi l^2':>10} {'#conc~2c/pi':>11} {'ev0':>12} {'ev1':>12} {'gap ev1-ev0':>13} {'evens simple?':>13}")
for lmbda in [1.5, 2.0, 3.0, 4.0]:
    c, w, V, n = prolate(lmbda, Nleg=max(400, int(6*2*np.pi*lmbda**2)))
    # even eigenfunctions: dominated by even-n Legendre. Identify parity by dominant component.
    parity = np.array([ (-1)**np.argmax(np.abs(V[:,k])) for k in range(len(w))])
    ev_even = w[parity>0]
    gap = ev_even[1]-ev_even[0]
    conc = 2*c/np.pi
    # simplicity: min consecutive spacing among lowest 8
    minspac = np.min(np.diff(w[:8]))
    print(f"{lmbda:>7.2f} {c:>10.2f} {conc:>11.1f} {ev_even[0]:>12.4e} {ev_even[1]:>12.4e} {gap:>13.4e} "
          f"{'YES' if minspac>1e-6 else 'NO':>13}", flush=True)

print("\nRefs: Connes Fact (line 1123) eigenvalues SIMPLE; ~4 lambda^2 = 2c/pi small nonzero angle eigenvalues.")
print("The GAP ev1-ev0 in the EVEN sector is what transfers even-simplicity from pw_lambda to QW_lambda (A_lambda)")
print("IF k_lambda approximates theta_x within that gap (Connes Step 2, the single remaining crux).")
