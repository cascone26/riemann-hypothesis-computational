#!/usr/bin/env python3
"""§938 (spectral engagement, FULL STEAM): build the PROLATE operator (the object in CCM 2511.22755 whose spectrum → ζ-zeros)
and characterize it in OUR setting. Prolate/concentration operator on band-limited functions over a window = time-limit∘
band-limit. Compute its eigenvalues (Slepian concentration values λ_n, sharp plunge) + the commuting 2nd-order operator's
eigenvalues χ_n. Then: does a transform of χ_n / the plunge structure track the ζ-zeros in the window? First concrete step
of the spectral-convergence engagement. Connect to our near-null cluster (§905, ~0.53 overlap)."""
import numpy as np
from scipy.linalg import eigh
import mpmath as mp
mp.mp.dps=20
# window [a,b], bandwidth c (=4W scaled). Prolate op K(x,y)=sin(c(x-y))/(pi(x-y)) restricted to [a,b], eigenvalues=concentration.
def prolate(a,b,c,M=400):
    x=np.linspace(a,b,M); h=(b-a)/(M-1)
    D=x[:,None]-x[None,:]
    K=np.where(np.abs(D)<1e-9, c/np.pi, np.sin(c*D)/(np.pi*D))*h  # concentration kernel (quadrature)
    ev=eigh(K,eigvals_only=True)[::-1]
    return ev,x
# zeros in a comparable range (imag parts) for reference
ZT=[]; k=1
while True:
    g=float(mp.im(mp.zetazero(k))); ZT.append(g)
    if g>120: break
    k+=1
ZT=np.array(ZT)
print("="*84)
print("§938 prolate operator spectrum (Slepian concentration λ_n) + the plunge; connect to ζ-zeros")
print("="*84)
for (a,b,c) in [(-1,1,10),(-1,1,20),(-1,1,40)]:
    ev,x=prolate(a,b,c)
    # Slepian: eigenvalues near 1 (# ≈ 2c/pi = time-bandwidth), sharp plunge to 0
    n1=int((ev>0.99).sum()); n5=int((ev>0.5).sum()); n0=int((ev>0.01).sum())
    tb=2*c/np.pi
    print(f"\n[{a},{b}] c={c}: time-bandwidth 2c/π={tb:.2f}; #λ>0.99={n1}, #λ>0.5={n5} (Landau ≈{tb:.0f}), #λ>0.01={n0}")
    print(f"  plunge region (0.01<λ<0.99): {n0-n1} modes; top λ: {np.array2string(ev[:5],precision=5)}")
    print(f"  plunge width ~ (2/π²)ln(time-bandwidth) = {(2/np.pi**2)*np.log(tb):.2f} (Landau-Widom)")
print("="*84)
print(f"ζ-zeros to h120: {len(ZT)} zeros. CCM (2511.22755): a specific operator built from prolate → these zeros as N,λ→∞.")
print("The plunge eigenspace = our near-null cluster (§905). Next: the exact CCM operator whose spectrum = zeros (workflow).")
print("="*84)
