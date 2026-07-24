#!/usr/bin/env python3
"""
lens_g_offline_damage.py — LEDGER §127 (Lens G core): quantify the DAMAGE a single off-line zero does to the localized
Weil form λ_a, as a function of (a, δ, γ₀), and verify the e^{2δa} scaling law. δ=β−½ is the off-line distance.
de Branges/Kreĭn rep: Q_W(v)=Σ_ρ v̂(z_ρ)² (z_ρ = zero in the z-variable, ρ=½+iz_ρ). On-line pair ±γ₀ ⟹ 2v̂(γ₀)²≥0.
Off-line rectangle {±γ₀±iδ} ⟹ 4Re[v̂(γ₀−iδ)²]. DAMAGE(v) = 2v̂(γ₀)² − 4Re[v̂(γ₀−iδ)²] (weighted 2/(γ₀²+δ²)); the max
over unit-norm v = the largest eigenvalue of the damage form = how far one off-line zero can push λ_a NEGATIVE.
Since v̂(γ₀−iδ)=∫v e^{iγ₀x}e^{δx}dx ~ e^{δa}, expect max-damage ~ e^{2δa}/γ₀².
"""
import numpy as np

def offline_max_damage(a, delta, gamma0, K=40, M=4000):
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]; w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    kk=np.arange(1,K+1); B=np.sin(np.pi*np.outer(x+a,kk)/(2*a))        # (M+1,K) sine basis, vanish at ±a
    def shat(z):                                                       # ŝ_k(z)=∫ sin_k(x) e^{izx} dx  (complex z ok)
        return (B*np.exp(1j*z*x)[:,None]*w[:,None]).sum(0)             # (K,)
    s0=shat(gamma0)                    # real-ish (v̂ at real γ0), K-vector
    u =shat(gamma0-1j*delta)           # complex, grows ~e^{δa}
    wt=2/(gamma0**2+delta**2)
    # damage form D = wt*(2 Re(s0 s0^T) - 4 Re(u u^T)) ; mass Mmat = <v_k,v_l>
    D = wt*(2*np.real(np.outer(s0,s0)) - 4*np.real(np.outer(u,u)))
    D=(D+D.T)/2
    Mmat = B.T@(w[:,None]*B); Mmat=(Mmat+Mmat.T)/2
    # max generalized eigenvalue of D rel Mmat = max damage per unit ||v||^2
    L=np.linalg.cholesky(Mmat); Li=np.linalg.inv(L)
    ev=np.linalg.eigvalsh(Li@D@Li.T)
    return ev[-1]                       # largest = max damage (λ_a can be pushed down by this)

print("Lens G: max damage to λ_a from ONE off-line zero at ½+δ+iγ₀, vs a. (verify e^{2δa}/γ₀² law)")
print(f"{'δ':>6} {'γ₀':>6} " + " ".join(f"{'a='+str(a):>11}" for a in [1.0,1.5,2.0,2.5,3.0]))
for delta,gamma0 in [(0.1,14.13),(0.2,14.13),(0.1,25.0),(0.05,14.13)]:
    row=[]
    for a in [1.0,1.5,2.0,2.5,3.0]:
        d=offline_max_damage(a,delta,gamma0)
        row.append(d)
    print(f"{delta:>6.2f} {gamma0:>6.1f} " + " ".join(f"{r:>11.4e}" for r in row), flush=True)
    # check scaling: damage(a)/damage(a-0.5) should be ~ e^{2δ·0.5}=e^{δ}
    ratios=[row[i]/row[i-1] for i in range(1,len(row)) if row[i-1]>0]
    print(f"       ratios a→a+0.5: {[round(r,3) for r in ratios]}  (predicted e^{{δ}}={np.exp(delta):.3f})",flush=True)
print("\nREAD: max-damage should GROW ~e^{2δa} (ratio per Δa=0.5 ≈ e^{δ}) ⟹ one off-line zero breaks λ_a≥0 at large a,")
print("confirming δ(a)↛0 (density alone can't give RH) BUT giving the explicit per-zero damage law for the δ(a) bound +")
print("the off-line-zero DETECTOR (measured λ_a < −damage(a,δ,γ) ⟹ off-line zero in that (δ,γ) window).")
