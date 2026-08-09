#!/usr/bin/env python3
"""§931 (function-field comparison, where RH is PROVEN): RH for a curve/F_q ⟺ Frobenius eigenvalues α_i on |α|=√q
⟺ the normalized power-sum Toeplitz matrix T_N=[c_{m-n}/q^{(m-n)/2}], c_k=Σα_i^k, is PSD (Bochner/moment positivity).
On-circle → T is a Gram matrix VV* (PSD, PROVEN by Hodge index). Off-circle → non-PSD. Demonstrate: (1) positivity ⟺
zeros-on-circle, (2) the razor (min eig → 0 as N grows to interpolate the 2g zeros), (3) off-line → negative margin.
Compare to our number-field razor: SAME sampling razor, but here positivity is a clean GRAM fact given RH."""
import numpy as np
np.set_printoptions(precision=4,suppress=True)
def toeplitz_margin(radii,angles,q,N):
    # normalized power sums s_k = Σ_i (r_i/√q)^k e^{i k θ_i}, k=-(N-1)..(N-1); Hermitian Toeplitz T[m,n]=s_{m-n}
    def s(k):
        return sum((r/np.sqrt(q))**abs(k)*np.exp(1j*k*t) for r,t in zip(radii,angles))
    T=np.array([[s(m-n) for n in range(N)] for m in range(N)])
    T=0.5*(T+T.conj().T)
    return float(np.linalg.eigvalsh(T)[0])  # min eigenvalue = the Weil margin
q=25; g=4  # genus-4 curve over F_25: 2g=8 Frobenius eigenvalues on |α|=5
np.random.seed=None
angles=[0.4,1.1,1.9,2.7]  # 4 angles; eigenvalues come in conjugate pairs θ,−θ → 8 total
ang8=angles+[-a for a in angles]
print("="*82)
print(f"§931 function-field Weil positivity (curve/F_{q}, genus {g}, 2g={2*g} Frobenius eigenvalues). RH ⟺ PSD")
print("="*82)
print(f"{'N (band-limit)':>15} {'min eig (on-circle=RH)':>22} {'min eig (ONE off: r=5.6)':>24}")
for N in [3,5,7,8,9,11]:
    on=toeplitz_margin([np.sqrt(q)]*(2*g), ang8, q, N)
    off_r=[np.sqrt(q)]*(2*g); off_r[0]=5.6  # move ONE eigenvalue OFF the circle (RH violation)
    off=toeplitz_margin(off_r, ang8, q, N)
    print(f"{N:>15} {on:>22.3e} {off:>24.3e}")
print("="*82)
print("on-circle: min eig ≥ 0 (PSD, RH holds — a GRAM matrix VV*), → RAZOR (→0) as N passes 2g=8 (interpolation).")
print("one off-circle: min eig goes NEGATIVE (RH violated ⟺ Weil positivity fails). The razor DETECTS off-line zeros.")
print("KEY: in function fields this positivity is PROVEN (Hodge index ⇒ eigenvalues on circle ⇒ Gram). Number fields")
print("lack the surface/Hodge-index (Deninger no-go) — so our arch-vs-prime razor is the analytic stand-in with no proof.")
print("="*82)
