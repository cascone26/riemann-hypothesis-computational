import numpy as np, os
from krein_fit import lanczos_jacobi, design

# Multiplicity detector: a DOUBLE zero of zeta at gamma_1 makes |zeta|^2 vanish to
# 4th order. Generalized-Jacobi: oscillation amplitude scales with the weight
# exponent => amplitude should ~2x. Simulate by multiplying the weight by the
# extra (t-g)^2(t+g)^2 factor (smooth normalization only shifts low-n terms).

A, M, K = 30.0, 8000, 70
G1 = 14.134725
GAMMAS = [14.134725, 21.022040, 25.010858]
th = (np.arange(M) + 0.5) * np.pi / M
ts = A * np.cos(th)
quad = A * np.sin(th) * (np.pi / M)
wz = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"zw_A{int(A)}_M{M}.npy"))
freqs = [2 * np.arccos(g / A) for g in GAMMAS]
ns = np.arange(1, K - 1)[4:]
_, bf = lanczos_jacobi(ts, quad, K)

def amps(w):
    _, bz = lanczos_jacobi(ts, w * quad, K)
    y = (bz[:-1] - bf[:-1])[5:]
    X = design(freqs, ns, False)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return [float(np.hypot(beta[1 + 2 * i], beta[2 + 2 * i])) for i in range(3)]

extra = ((ts - G1) ** 2) * ((ts + G1) ** 2) / G1 ** 4     # order 2 -> 4 at +-gamma_1
a1 = amps(wz)
a2 = amps(wz * extra)
print("amplitudes at [g1, g2, g3]:")
print(f"  simple zeros (real zeta): {[round(x, 4) for x in a1]}")
print(f"  g1 doubled  (order 4):    {[round(x, 4) for x in a2]}")
print(f"  ratio at g1: {a2[0] / a1[0]:.3f}  (prediction ~2.0)   g2,g3 control ratios: "
      f"{a2[1] / a1[1]:.3f}, {a2[2] / a1[2]:.3f} (prediction ~1.0)")
