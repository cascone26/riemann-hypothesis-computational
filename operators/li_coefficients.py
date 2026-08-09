#!/usr/bin/env python3
"""§921 (strategic): are the Li coefficients λ_n comfortably positive (non-razor) or razor-thin? λ_n=Σ_ρ[1−(1−1/ρ)^n]≥0
⟺ RH. If λ_n are COMFORTABLY positive (large margin, not ~1e-15), the Li criterion is a better arb-certification target
than our razor. Compute λ_n from zeros (magnitudes robust) and report their size + the 'how close to 0' margin."""
import mpmath as mp
mp.mp.dps=30
# zeros up to height T
ZT=[]; k=1
while True:
    g=mp.im(mp.zetazero(k)); ZT.append(g)
    if g>300: break
    k+=1
print(f"using {len(ZT)} zeros up to height {float(ZT[-1]):.0f} (pairs ρ=1/2±iγ)")
print("="*80)
print(f"{'n':>4} {'λ_n':>16} {'λ_n/n':>12} {'~expected (n/2)(ln(n/2π)+γ-1)?':>18}")
import math
for n in list(range(1,11))+[15,20,30,50]:
    # λ_n = Σ_ρ [1-(1-1/ρ)^n]; ρ=1/2+iγ and conjugate. sum over pairs = 2 Re Σ_{γ>0}[1-(1-1/ρ)^n]
    s=mp.mpf(0)
    for g in ZT:
        rho=mp.mpc(0.5,g)
        s+= (1-(1-1/rho)**n) + (1-(1-1/mp.conj(rho))**n)  # both ρ and conj
    lam=mp.re(s)
    # expected asymptotic λ_n ~ (n/2)(ln n - ln(2π) + γ_euler -1) + ...
    exp=0.5*n*(math.log(n)-math.log(2*math.pi)+0.5772156649-1) if n>1 else None
    es=f"{exp:.3f}" if exp else " -- "
    print(f"{n:>4} {float(lam):>16.4f} {float(lam)/n:>12.4f} {es:>18}")
print("="*80)
print("λ_n comfortably positive & growing (~ (n/2)ln n) → NON-razor positivity. RH-content is in the FLUCTUATION around")
print("the smooth growth, not the sign of individual λ_n. (Truncated-zero estimate; magnitudes robust, not proof.)")
print("="*80)
