import numpy as np
# Li coefficients from zeros: λ_n = Σ_ρ [1-(1-1/ρ)^n], RH ⟺ all λ_n≥0 (Li's criterion).
# Then the MOMENT/JACOBI structure: build Hankel from power sums, get Jacobi b_n^2 (Hankel det ratios),
# check b_n^2>0 (RH-consistent) + Remling convergence (do a_n,b_n settle?).
G = np.load('operators/zeros10k.npy')          # γ_k > 0 (real, RH-assumed positions)
rho = 0.5 + 1j*G                               # ρ = 1/2 + iγ (and conjugate)

def li_coeffs(nmax):
    lam=[]
    for n in range(1,nmax+1):
        term = 1 - (1 - 1/rho)**n
        lam.append(2*np.sum(term.real))        # +conjugate = 2*Re
    return np.array(lam)

NL=60
lam = li_coeffs(NL)
neg = np.where(lam<0)[0]
print(f"Li coefficients λ_1..{NL} from 10k zeros:")
print(f"  λ_1..8: {np.round(lam[:8],5).tolist()}")
print(f"  all positive? {len(neg)==0}  (Li's criterion: RH ⟺ all λ_n>0);  #neg={len(neg)}")
print(f"  λ_n growth: λ_30={lam[29]:.3f} λ_45={lam[44]:.3f} λ_60={lam[59]:.3f} (grow ~ (n/2)(log n - ...) under RH)")

# power sums σ_k = Σ_ρ ρ^{-k} (moments), build Hankel & Jacobi
def sigma(k): return 2*np.sum((rho**(-k)).real)
NM=40
sig = np.array([sigma(k) for k in range(1,2*NM+2)])
# Hankel-based Jacobi via Chebyshev/Stieltjes on modified moments is unstable; use LDL of Hankel for b_n^2 sign
H = np.array([[sig[i+j] for j in range(NM)] for i in range(NM)])   # Hankel of shifted moments
# leading principal minors D_k; Jacobi b_k^2 = D_{k+1}D_{k-1}/D_k^2 ; sign(D_k) tells inertia
dets=[1.0]
for k in range(1,NM+1):
    dets.append(np.linalg.det(H[:k,:k]))
dets=np.array(dets)
signchanges = sum(1 for i in range(1,len(dets)) if dets[i]*dets[i-1]<0)
print(f"\nHankel(σ) principal-minor sign-changes = {signchanges}  (=#neg squares of the moment form; RH-consistent = controlled/small)")
print(f"  first minors signs: {[int(np.sign(d)) for d in dets[:12]]}")
