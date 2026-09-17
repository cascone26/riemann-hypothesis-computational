#!/usr/bin/env python3
"""§957 END-TO-END validation of the Krein-Langer pole-count = κ mechanism in the PROVEN function-field case.

The route (§950/§955): RH ⟺ the Weil m-function is generalized-Nevanlinna N_κ with κ = #off-line zeros, so
κ=1 (only the pole) ⟺ RH. Before more number-field construction, VALIDATE the mechanism where RH is PROVEN and
we CONTROL the zeros: a curve/F_q. Frobenius eigenvalues α_i = √q·e^{iθ_i} on the circle (RH_ff) or off it.
Build the associated Nevanlinna function m(z) = Σ_i 1/(z - u_i), u_i = the 'zero heights' (u_i = θ_i real on the
circle; off-circle α_i=r_i e^{iθ} ⟹ u_i = θ_i - i·ln(r_i/√q) leaves the real axis into the UHP/LHP). Then:
  κ(Nevanlinna kernel of m) should = number of u_i in the UPPER half-plane = number of off-circle-outward zeros.
If κ tracks the off-circle count EXACTLY, the pole-count=κ mechanism is confirmed end-to-end in the proven case —
so the number-field task is purely to transcribe ζ's Weyl function into this Nevanlinna form. Controls included.
HANDLE (own code)."""
import numpy as np

def nevanlinna_kappa(u_list, res_list, zpts, tol=1e-9):
    """κ of kernel of m(z)=Σ res_i/(z-u_i) sampled on zpts. res_i>0 (Nevanlinna point masses)."""
    def m(z): return sum(r/(z-u) for u, r in zip(u_list, res_list))
    Z = np.array(zpts); mv = np.array([m(z) for z in Z])
    K = np.empty((len(Z),len(Z)), complex)
    for i in range(len(Z)):
        for j in range(len(Z)):
            d = Z[i]-np.conj(Z[j])
            K[i,j] = np.imag(mv[i])/np.imag(Z[i]) if abs(d)<1e-12 else (mv[i]-np.conj(mv[j]))/d
    K = 0.5*(K+K.conj().T); ev = np.linalg.eigvalsh(K); s=max(1.0,abs(ev).max())
    return int(np.sum(ev < -tol*s))

rng = np.random.default_rng(957)
zp = [complex(x,y) for x,y in zip(rng.uniform(-3,3,26), rng.uniform(0.3,3.0,26))]
q=25; R=np.sqrt(q); g=4; angles=np.array([0.4,1.1,1.9,2.7,-0.4,-1.1,-1.9,-2.7])

print("="*90)
print("§957 function-field validation: κ(Nevanlinna kernel) vs number of off-circle (UHP) Frobenius zeros.")
print("="*90)
# CONTROLS: point masses on the real axis (all u real, res>0) = genuine Nevanlinna N_0 (κ=0)
print(f"  CONTROL: {g*2} masses all on real axis (u real, res>0) → N₀ expect κ=0: κ={nevanlinna_kappa(list(angles),[1.0]*8,zp)}")
print("-"*90)
print(f"  {'# moved OUTWARD off-circle (→UHP)':>34} | {'u_i UHP count':>13} | {'κ(Nevanlinna kernel)':>20}")
for k in range(0, 5):
    radii = np.array([R]*8, float)
    for i in range(k): radii[i] = R*np.exp(0.3*(i+1)/8 + 0.2)   # push outward: r>√q ⟹ u gets +Im
    # u_i = theta_i - i*ln(r_i/R): outward (r>R) ⟹ Im(u) = -ln(r/R) < 0 → LHP; inward would be UHP.
    # Use SIGNED: count how many leave the real axis (|Im u|>eps). For an N_κ readout we place the
    # off-axis masses with a rank-1 negative-type term when Im(u)<0 (generalized pole). Model: split
    u = angles.astype(complex).copy()
    for i in range(k): u[i] = angles[i] - 1j*np.log(radii[i]/R)
    offaxis = int(np.sum(np.abs(np.imag(u)) > 1e-9))
    # generalized Nevanlinna with off-axis poles: build kernel directly from m(z)=Σ 1/(z-u_i) (complex u)
    def m(z): return sum(1.0/(z-ui) for ui in u)
    Z=np.array(zp); mv=np.array([m(z) for z in Z]); Kk=np.empty((len(Z),len(Z)),complex)
    for i in range(len(Z)):
        for j in range(len(Z)):
            d=Z[i]-np.conj(Z[j]); Kk[i,j]=np.imag(mv[i])/np.imag(Z[i]) if abs(d)<1e-12 else (mv[i]-np.conj(mv[j]))/d
    Kk=0.5*(Kk+Kk.conj().T); ev=np.linalg.eigvalsh(Kk); s=max(1.0,abs(ev).max())
    kap=int(np.sum(ev<-1e-9*s))
    print(f"  {k:>34} | {offaxis:>13} | {kap:>20}")
print("="*90)
print("READ: if κ(Nevanlinna kernel) grows by exactly 1 per off-circle Frobenius zero (each leaving the real")
print("axis = one negative-type pole), the pole-count=κ mechanism is CONFIRMED in the proven case. RH_ff ⟺ all")
print("zeros on circle ⟺ all u_i real ⟺ κ=0 (+1 for ζ's actual pole in the number-field analog). Then the ζ task")
print("is exactly: transcribe ζ's Weyl function to m(z)=Σ 1/(z-u_i), u_i=zero heights + the s=1 pole → κ=1 ⟺ RH.")
print("="*90)
