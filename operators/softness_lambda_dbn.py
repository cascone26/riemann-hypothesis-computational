import numpy as np, os

# SUPERPRISM #5: quantify "how barely true" — connect the operator softness to Lambda_dBN.
# de Bruijn-Newman: H_t(x) with zeros; Lambda = inf t s.t. all zeros real. RH <=> Lambda<=0;
# Rodgers-Tao Lambda>=0 => Lambda=0. The heat flow moves zeros: gamma_k(t) with
#   d gamma_k/dt = sum_{j!=k} 2/(gamma_k - gamma_j)   (the DBN/Calogero-type flow, real-line).
# Backward flow (t<0) can push neighboring zeros to collide -> complex (off-line). The MARGIN to
# collision at t=0 = how barely RH holds locally. Measure: min over k of the "collision time"
# tau_k ~ (gamma_{k+1}-gamma_k)^2 / (relative repulsion) and its distribution/scaling with height.
# This is a concrete, computable softness law from the actual zeros.

DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros10k.npy"))
g = Z.astype(float)
# local repulsion velocity difference between adjacent zeros under the flow
# v_k = sum_{j!=k} 2/(g_k-g_j); collision of (k,k+1) when gap closes: gap' = v_{k+1}-v_k
sp = np.diff(g)
# approximate v_k with a windowed sum (dominant near neighbors)
def vel(i, W=60):
    lo, hi = max(0, i - W), min(len(g), i + W + 1)
    idx = np.arange(lo, hi); idx = idx[idx != i]
    return np.sum(2.0 / (g[i] - g[idx]))
V = np.array([vel(i) for i in range(len(g))])
gap = sp
dgap_dt = V[1:] - V[:-1]   # rate of change of gap under forward flow
# forward flow: if dgap/dt>0 gap widens (stable); collision needs dgap/dt<0 backward
# "collision time" going backward: tau = gap / max(dgap_dt, 0) (time to close going backward = -t)
tau = np.where(dgap_dt > 0, gap / dgap_dt, np.inf)
h = 0.5 * (g[1:] + g[:-1])
print("softness law: backward-collision time tau_k = gap/(d gap/dt) vs height (the local RH margin):")
for lo, hi in [(0, 500), (500, 1500), (1500, 3000), (3000, 6000), (6000, 10000)]:
    m = (h >= lo) & (h < hi)
    finite = tau[m][np.isfinite(tau[m])]
    if len(finite):
        print(f"  height [{lo:>5},{hi:>5}): median tau={np.median(finite):.3f}, min tau={finite.min():.4f}, "
              f"mean gap={gap[m].mean():.3f}")
print("\nmin tau across all = the tightest local margin = closest approach to a DBN collision (Lambda proxy).")
mn = np.nanmin(tau[np.isfinite(tau)])
print(f"global min backward-collision time = {mn:.4f} (at height {h[np.argmin(np.where(np.isfinite(tau),tau,np.inf))]:.1f})")
print("does min-tau SHRINK with height (margin -> 0, Lambda=0 saturation) or stay bounded?")
