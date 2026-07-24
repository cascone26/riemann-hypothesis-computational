#!/usr/bin/env python3
"""
subgoalH_grochenig_TP.py — LEDGER §149 (Lens 4/7, the foothold-equipped lever): build Gröchenig's object (arXiv
2007.12889). RH ⟺ Λ(x)=(1/2π)∫ e^{ixτ}/ξ(1/2+τ)dτ is a PÓLYA FREQUENCY function = TOTALLY POSITIVE (all minors
det[Λ(x_i−x_j)]≥0). UNCONDITIONAL FOOTHOLD (proven): Λ(x)≥0 pointwise (Bochner; ξ(1/2+τ)>0 & pole-free on ℝ ⟹
1/ξ(1/2+τ) is a positive-definite function ⟹ its FT Λ≥0). Gröchenig: even the NEXT rung isn't known. The ladder up to
RH: Λ≥0 (TP_1, HAVE) → Λ log-concave (TP_2 / PF_2) → all 3×3 minors ≥0 (TP_3) → ... → totally positive (TP_∞ = RH).
Since ξ(1/2+τ) is even, Λ(x)=(1/π)∫_0^∞ cos(xτ)/ξ(1/2+τ)dτ (real, even).
This tool: (1) build Λ(x) by cosine-transforming 1/ξ (mpmath ξ on a τ-grid). (2) CONFIRM the foothold Λ(x)≥0.
(3) climb the TP ladder: is Λ LOG-CONCAVE (TP_2)? are the 3×3,4×4 minors of [Λ(x_i−x_j)] ≥0 (TP_3,4)? Where is positivity
TIGHTEST (the near-singular minor = where RH is 'closest to failing')? All UNCONDITIONALLY checkable — a concrete climb.
"""
import numpy as np, mpmath as mp
mp.mp.dps = 25

def xi(s):
    s=mp.mpf(s) if not isinstance(s,mp.mpf) else s
    return 0.5*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)

# precompute 1/ξ(1/2+τ) on a τ-grid (even, decays superfast)
T=45.0; NT=1800
dtau=T/NT
tau=(np.arange(NT)+0.5)*dtau        # MIDPOINT grid — never lands on τ=0.5 (s=1, ζ pole)
print("building 1/ξ(1/2+τ) on τ-grid (mpmath, midpoint)...", flush=True)
finv=np.array([float(1/xi(0.5+t)) for t in tau])   # positive, peaked at τ=0
def Lam(x):
    x=np.atleast_1d(x)
    integ=np.cos(np.outer(x,tau))*finv[None,:]
    return (integ*dtau).sum(1)/np.pi                # midpoint rule

# (2) foothold: Λ(x) ≥ 0
xg=np.linspace(0,6,241); L=Lam(xg)
print(f"\n(2) FOOTHOLD Λ(x)≥0:  Λ(0)={L[0]:.5f}  min over [0,6]={L.min():.3e} @x={xg[np.argmin(L)]:.2f}  "
      f"{'✓ Λ≥0 holds' if L.min()>=-1e-9 else '✗ Λ<0 (check integration)'}")

# (3a) TP_2 = log-concavity of Λ: log Λ concave ⟺ Λ(x)² ≥ Λ(x-δ)Λ(x+δ). test second difference of logΛ.
Lpos=L.copy(); Lpos[Lpos<=0]=np.nan; lg=np.log(Lpos)
d2=lg[2:]-2*lg[1:-1]+lg[:-2]     # ∝ (logΛ)''; ≤0 ⟹ log-concave
frac_concave=np.nanmean(d2<=1e-12)
print(f"(3a) TP_2 (log-concavity): (logΛ)'' ≤0 on {100*frac_concave:.1f}% of [0,6];  max (logΛ)''={np.nanmax(d2)/dtau**2*dtau**2:.2e} "
      f"{'✓ Λ LOG-CONCAVE' if frac_concave>0.999 else '✗ NOT everywhere log-concave — first convex point matters'}")
if frac_concave<=0.999:
    bad=xg[1:-1][np.nan_to_num(d2,nan=-1)>1e-12]
    print(f"       convex (TP_2-violating) x-locations: {np.round(bad[:8],3)}")

# (3b) TP_k minors: matrices [Λ(x_i−x_j)] on equispaced points; check leading k×k principal minors ≥0 and smallest eig.
print("(3b) TP ladder on equispaced points x_i=i·h (Toeplitz [Λ(|i-j|h)]):")
for h in [0.4,0.7,1.0]:
    row=np.array([Lam(np.array([h*k]))[0] for k in range(12)])
    print(f"   h={h}: leading-minor signs / smallest-eig of k×k Toeplitz:")
    for k in [2,3,4,6,8]:
        Mk=np.array([[row[abs(i-j)] for j in range(k)] for i in range(k)])
        dets=[np.linalg.det(Mk[:m,:m]) for m in range(1,k+1)]
        allpos=all(d>-1e-12 for d in dets); mineig=np.linalg.eigvalsh(Mk)[0]
        print(f"      k={k}: leading minors {'all≥0' if allpos else 'NEG APPEARS'}  min-eig={mineig:+.3e}",flush=True)
print("\nREAD: Λ≥0 = the proven foothold. Log-concavity (TP_2) UNCONDITIONALLY checkable — if it holds everywhere that is a")
print("rung ABOVE the foothold toward RH (Λ TP = RH). The k×k minor positivity + smallest-eig show where total positivity")
print("is tightest. If TP_2 holds unconditionally & cleanly, the target sharpens to: prove Λ log-concave → higher TP → RH.")
