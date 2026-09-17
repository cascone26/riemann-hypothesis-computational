import numpy as np
z = np.load('operators/zeros10k.npy')

# Spacing ratio r_n = min(s_n,s_{n-1})/max(...) is UNFOLDING-FREE (Oganesyan-Huse / Atas)
# Reference means: GUE=0.5996, GSE=0.6744, GOE=0.5307, Poisson=0.3863
d = np.diff(z)
r = np.minimum(d[1:], d[:-1]) / np.maximum(d[1:], d[:-1])
print(f"N zeros = {len(z)},  <r> (all) = {r.mean():.4f}")

# by window to check stability / running to high T
for lo,hi in [(0,2000),(2000,5000),(5000,10000)]:
    rr = r[lo:hi]
    print(f"  window [{lo}:{hi}]  <r> = {rr.mean():.4f}  (n={len(rr)})")

print("\nReference: Poisson=0.3863  GOE=0.5307  GUE=0.5996  GSE=0.6744")

# Number variance / rigidity is the other discriminant; do the cheap unfolded NN-spacing dist too
def Nsmooth(T):
    return (T/(2*np.pi))*np.log(T/(2*np.pi)) - T/(2*np.pi) + 7/8
u = Nsmooth(z)                      # unfolded positions (mean spacing -> 1)
s = np.diff(u)
print(f"\nUnfolded mean spacing = {s.mean():.4f} (should be ~1)")
# small-s behavior: fraction of spacings below thresholds (level repulsion order)
for thr in [0.1,0.2,0.3,0.5]:
    print(f"  P(s<{thr}) = {(s<thr).mean():.4f}")
# GUE predicts P(s)~ (32/pi^2)s^2 small s  -> beta=2 repulsion; GSE ~ s^4 (beta=4)
# fit exponent via log-log of cumulative small-s
mask = s < 0.4
import numpy as np
sv = np.sort(s[mask])
cdf = np.arange(1,len(sv)+1)/len(s)
# CDF ~ s^(beta+1); slope of log cdf vs log s
lo = (sv>0.05)&(sv<0.35)
A = np.polyfit(np.log(sv[lo]), np.log(cdf[lo]), 1)
print(f"\nSmall-s CDF slope = {A[0]:.2f}  => repulsion exponent beta ~ {A[0]-1:.2f}  (GUE beta=2, GSE beta=4)")
