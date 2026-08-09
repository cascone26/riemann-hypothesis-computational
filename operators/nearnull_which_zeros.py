#!/usr/bin/env python3
"""§858: WHICH zeros does the RH-critical near-null mode struggle to vanish at? The zero side (validated)
is Q_zero(v)=2 sum_gamma F_v(gamma)^2 (test function at the zero heights, a manifest sum of squares). The
near-null modes vanish at the window's zeros; the per-zero contribution 2*F_v(gamma)^2 shows which zeros
carry the residual. Extract the margin eigenvector at [14,44], show per-zero contributions."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import l2_fejer as lf
from scipy.linalg import eigh
a, b = 14.0, 44.0
N = 15
CEN = np.linspace(a, b, N)
d = lf.build(CEN)
full, G = d['full'], d['overlap']
w, V = eigh(full, G)               # ascending; V columns G-orthonormal
vmin = V[:, 0]                     # the margin eigenvector
print(f"§858 [14,44] 4W={4*lf.W:.1f}, N={N}: margin(min eig)={w[0]:+.4e}; next eigs {w[1]:.3e},{w[2]:.3e},{w[3]:.3e}")
# test function at zeros: F_v(gamma) = sum_c v_c F_c(gamma)
FZ = np.array([lf.Frow(c, lf.ZEROS) for c in CEN])   # N x nzeros
Fv_at_zeros = vmin @ FZ                               # value at each zero
contrib = 2.0 * Fv_at_zeros**2                        # per-zero contribution to Q_zero(v)
Qzero = contrib.sum()
print(f"Q_zero(vmin) = 2*sum F_v(gamma)^2 = {Qzero:+.4e}  (should ~= margin {w[0]:+.4e}, explicit-formula check)")
inwin = (lf.ZEROS >= a) & (lf.ZEROS <= b)
print(f"in-window zeros: {int(inwin.sum())}, out-window: {int((~inwin).sum())}")
print(f"{'gamma':>10} {'in?':>4} {'F_v(gamma)':>12} {'contribution':>13} {'% of Qzero':>10}")
order = np.argsort(-contrib)
for i in order[:12]:
    g = lf.ZEROS[i]
    print(f"{g:>10.3f} {'IN' if inwin[i] else 'out':>4} {Fv_at_zeros[i]:>+12.4e} {contrib[i]:>13.4e} {100*contrib[i]/Qzero:>9.1f}%")
print(f"in-window share of residual: {100*contrib[inwin].sum()/Qzero:.1f}%")
