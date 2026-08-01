import numpy as np, os
from krein_fit import lanczos_jacobi, design

# Operator-side razor: an off-line zero at beta=1/2+delta lifts the weight's double
# zero at t=gamma to |c|^2*((t-gamma)^2+delta^2). Simulate on gamma_1=14.1347 and
# measure the collapse of the Jacobi oscillation amplitude at 2*arccos(gamma_1/A).
# amplitude(delta)/amplitude(0) = the string's detection stiffness.

A, M, K = 30.0, 8000, 70
G1 = 14.134725
GAMMAS = [14.134725, 21.022040, 25.010858]

th = (np.arange(M) + 0.5) * np.pi / M
ts = A * np.cos(th)
quad = A * np.sin(th) * (np.pi / M)
wz = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"zw_A{int(A)}_M{M}.npy"))
_, bf = lanczos_jacobi(ts, quad, K)

freqs = [2 * np.arccos(g / A) for g in GAMMAS]
ns = np.arange(1, K - 1)[4:]

def amp_at_zero1(w):
    _, bz = lanczos_jacobi(ts, w * quad, K)
    y = (bz[:-1] - bf[:-1])[5:]
    X = design(freqs, ns, False)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return np.hypot(beta[1], beta[2])       # amplitude at gamma_1 frequency

# off-line lift acts on BOTH t=+-gamma_1 (weight even, zeros come in pairs)
def lift(delta):
    f = np.ones_like(ts)
    for g in (G1, -G1):
        f *= ((ts - g) ** 2 + delta ** 2) / np.maximum((ts - g) ** 2, 1e-18)
    return wz * f

a0 = amp_at_zero1(wz)
print(f"A={A}: baseline amplitude at 2*arccos(g1/A) = {a0:.4f}")
print(f"{'delta':>7} {'amp':>8} {'amp/amp0':>9}")
for d in [0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]:
    a = amp_at_zero1(lift(d))
    print(f"{d:>7.3f} {a:>8.4f} {a / a0:>9.4f}", flush=True)
print("\ncollapse curve = the Krein string's off-line detection stiffness at gamma_1;")
print("compare shape vs the Weil-form razor's ~delta^2 indefiniteness (474b).")
