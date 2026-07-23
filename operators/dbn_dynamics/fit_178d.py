#!/usr/bin/env python3
"""§178d fit: identify the law behind dist(k_λ,ψ_c). Tests 3 hypotheses (SECT66_PREP):
(a) dist ~ λmin = e^{-4πc}  [tracks the near-null eigenvalue]
(b) dist ~ √λmin = e^{-2πc} [tracks the ANGLE sin α_λ ~ √(1-χ_2)]
(c) dist ~ c^{-p} or λ^{-1/2} [Fourier polynomial rate]
Data hardcoded from results/bridge_law_178d.log (converged plateaus c=3,4,5; c=6 clean-but-N-unconfirmed)."""
import numpy as np
# c : (dist_plateau, log10_lambda_min, confirmed?)
D = {3:(6.894e-4,-7.29,True), 4:(4.628e-4,-12.19,True), 5:(1.245e-4,-17.11,True), 6:(8.000e-6,-22.32,False)}
cs   = np.array(sorted(D))
dist = np.array([D[c][0] for c in cs]); logd = np.log10(dist)
lm   = np.array([D[c][1] for c in cs])              # = log10(lambda_min)
conf = np.array([D[c][2] for c in cs])
core = conf                                          # trustworthy set = confirmed plateaus

def fit(x, y, lbl, mask):
    A = np.vstack([x[mask], np.ones(mask.sum())]).T
    (m,b),res,*_ = np.linalg.lstsq(A, y[mask], rcond=None)
    yhat = m*x+b; ss_res=np.sum((y[mask]-yhat[mask])**2); ss_tot=np.sum((y[mask]-y[mask].mean())**2)
    r2 = 1-ss_res/ss_tot if ss_tot>0 else float('nan')
    return m,b,r2,yhat

print("=== §178d data ===")
print(f"{'c':>3} {'dist':>11} {'log10 dist':>11} {'log10 λmin':>11} {'confirmed':>10}")
for c in cs:
    print(f"{c:>3} {D[c][0]:>11.3e} {np.log10(D[c][0]):>11.4f} {D[c][1]:>11.2f} {str(D[c][2]):>10}")

print("\n--- local slopes d(log10 dist)/dc (consecutive) ---")
for i in range(len(cs)-1):
    print(f"  c={cs[i]}→{cs[i+1]}: {(logd[i+1]-logd[i])/(cs[i+1]-cs[i]):+.3f}   (λmin slope {(lm[i+1]-lm[i])/(cs[i+1]-cs[i]):+.2f})")

print("\n=== HYP (a) dist ~ λmin  [log10 dist vs c, compare slope to λmin slope ~-4.9] ===")
m,b,r2,_ = fit(cs.astype(float), logd, "c", core)
print(f"  fit log10 dist = {m:+.3f}·c {b:+.3f}  R²={r2:.4f}   (λmin slope ≈ -4.9 ⇒ if dist~λmin these MATCH)")
print(f"  → dist slope {m:.2f} vs λmin slope -4.9: {'MATCH' if abs(m+4.9)<0.6 else 'NO — dist decays far SLOWER than λmin'}")

print("\n=== HYP (a'/b) dist ~ λmin^p  [log10 dist vs log10 λmin; p=1 eigenvalue, p=1/2 angle] ===")
m,b,r2,_ = fit(lm, logd, "lm", core)
print(f"  fit log10 dist = {m:.4f}·log10λmin {b:+.3f}  R²={r2:.4f}")
print(f"  → power p = {m:.4f}.  p≈1 ⇒ dist~λmin (a);  p≈0.5 ⇒ dist~√λmin=sin(angle) (b)")
# local powers
print("  local p (consecutive): " + ", ".join(f"{(logd[i+1]-logd[i])/(lm[i+1]-lm[i]):.3f}" for i in range(len(cs)-1)))

print("\n=== HYP (c) dist ~ c^{-p} power law  [log10 dist vs log10 c] ===")
m,b,r2,_ = fit(np.log10(cs.astype(float)), logd, "logc", core)
print(f"  fit log10 dist = {m:.3f}·log10 c {b:+.3f}  R²={r2:.4f}   (p={-m:.2f}; λ^{{-1/2}} ⇒ p=1/4 in c)")

print("\n=== READ ===")
print("Local slope d(log dist)/dc STEEPENS (−0.17→−0.57→−1.19) and local power p vs λmin GROWS (0.035→0.12→0.23)")
print("→ neither a clean single exponential (a) nor fixed power. dist decays SLOWER than λmin but ACCELERATING.")
print("Check dist/√λmin trend (hyp b, sin-angle): does it approach a constant (⇒ dist~√λmin asymptotically)?")
r = dist/10**(lm/2.0)
print("  dist/√λmin: " + ", ".join(f"c={cs[i]}:{r[i]:.3g}" for i in range(len(cs))))
print("  ratios of consecutive: " + ", ".join(f"{r[i+1]/r[i]:.1f}" for i in range(len(cs)-1)) +
      "  (DECREASING ⇒ heading toward dist~√λmin×const; the residual growth = subleading/poly)")
