"""CORRECT the class-C control: Srednicki's class = Altland-Zirnbauer BdG class C (symplectic PHS),
NOT Wigner-Dyson GSE. AZ class C: H = -σ_y Hᵀ σ_y (particle-hole), realized as
  H = [[a, b],[b†, -aᵀ]],  a hermitian, b symmetric.
Spectrum ±symmetric about 0. KEY prediction: BULK level statistics are GUE-like (β=2, <r>≈0.60) —
so bulk <r> CANNOT distinguish class A from class C — while the NEAR-E=0 density has a soft gap ρ(E)~E².
Verify both, then state what this means for testing 'ζ operator is class C'."""
import numpy as np
rng = np.random.default_rng(7)

def classC(n):  # 2n x 2n AZ class-C Hamiltonian
    a = (rng.standard_normal((n,n))+1j*rng.standard_normal((n,n))); a=(a+a.conj().T)/2
    b = (rng.standard_normal((n,n))+1j*rng.standard_normal((n,n))); b=(b+b.T)/2   # symmetric
    top=np.hstack([a,b]); bot=np.hstack([b.conj().T,-a.T])
    return np.vstack([top,bot])
def gue(n):
    A=(rng.standard_normal((n,n))+1j*rng.standard_normal((n,n)))/np.sqrt(2); return (A+A.conj().T)/2
def r_bulk(ev):
    ev=np.sort(ev); m=len(ev); ev=ev[m//4:3*m//4]; s=np.diff(ev)
    rr=np.minimum(s[1:],s[:-1])/np.maximum(s[1:],s[:-1]); return rr.mean()

n=200; T=80
rc=[]; rg=[]
allC=[]
for _ in range(T):
    ev=np.linalg.eigvalsh(classC(n)); rc.append(r_bulk(ev)); allC.append(ev)
    rg.append(r_bulk(gue(2*n)))
rc=np.array(rc); rg=np.array(rg)
print("BULK spacing-ratio <r> (unfolding-free):")
print(f"  AZ class C (correct): <r> = {rc.mean():.4f} ± {rc.std():.4f}")
print(f"  GUE      (class A):   <r> = {rg.mean():.4f} ± {rg.std():.4f}")
print(f"  => bulk separation between A and C: {abs(rc.mean()-rg.mean()):.4f}  "
      f"({abs(rc.mean()-rg.mean())/np.sqrt(rc.std()**2+rg.std()**2):.1f}σ)")

# NEAR-E=0 density: pool eigenvalues near 0 across trials, unfold by bulk mean-spacing, fit ρ(E)~|E|^k
allC=np.concatenate(allC)
# local mean spacing at the origin from the class-C spectra
mid = np.sort(np.abs(allC)); 
# density in small-E bins (symmetric), normalized by the bulk density
bulk_density = len(allC)/(allC.max()-allC.min())
edges=np.linspace(0, 3.0, 16)
cnt,_=np.histogram(np.abs(allC), bins=edges)
width=np.diff(edges); centers=(edges[:-1]+edges[1:])/2
rho = cnt/(width*len(allC))    # per-trial-ish density of |E|
# fit small-E exponent
lo=(centers>0.15)&(centers<1.2)
k=np.polyfit(np.log(centers[lo]), np.log(rho[lo]),1)[0]
print(f"\nNEAR-E=0 density of states, AZ class C:  ρ(E) ~ |E|^{k:.2f}   (class C predicts k≈2; class A gives k≈0)")
print("Interpretation: BULK is GUE-like for BOTH classes (can't distinguish); the class-C fingerprint")
print("is ONLY the near-origin soft gap ρ~E². A single ζ has ONE spectrum — you cannot measure a")
print("near-origin density EXPONENT from one realization; that is a FAMILY (Katz-Sarnak) or ALGEBRAIC")
print("(operator-structure) statement, not a ζ-bulk observable.")
