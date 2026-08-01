import numpy as np, os
import mpmath as mp
from krein_fit import lanczos_jacobi
from scipy.fft import dct

# OUT-OF-SAMPLE validation of the string dictionary at untouched window A=50:
#   db_pred(n) = sum_k C*sin(theta_k)*cos(2n theta_k - 2 psi(theta_k) + c0)
# with C=0.6777 and scale s=-2 FROZEN from other windows; only the single window
# offset c0 is scanned. Metric: max_c0 corr(db_pred, db_meas) on n in [10, K-2],
# null = same with 500 random phase sets.

mp.mp.dps = 15
A, M, K = 50.0, 12000, 115
C_FROZEN, S_FROZEN = 0.6777, -2.0
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
gam = [float(g) for g in Z if g < A - 0.5]
print(f"A={A}: {len(gam)} zeros (window untouched by any fit)", flush=True)

th = (np.arange(M) + 0.5) * np.pi / M
ts = A * np.cos(th)
quad = A * np.sin(th) * (np.pi / M)
fn = os.path.join(DIR, f"zw_A50_M{M}.npy")
if os.path.exists(fn):
    wz = np.load(fn)
else:
    wz = np.empty(M)
    for i, t in enumerate(ts):
        wz[i] = float(abs(mp.zeta(mp.mpc(0.5, t)))) ** 2
        if i % 2000 == 0:
            print(f"  zeta grid {i}/{M}", flush=True)
    np.save(fn, wz)

_, bz = lanczos_jacobi(ts, wz * quad, K)
_, bf = lanczos_jacobi(ts, quad, K)
db = (bz[:-1] - bf[:-1])
thk = np.array([np.arccos(g / A) for g in gam])
L = 0.5 * np.log(np.maximum(wz, 1e-300))
a = dct(L, type=2, norm=None) / M
KF = min(3000, M - 1)
ks = np.arange(1, KF)
psi = np.array([np.sum(a[1:KF] * np.sin(ks * t)) for t in thk])

ns = np.arange(len(db))
sel = (ns >= 10) & (ns <= K - 3)
y = db[sel]

c0s = np.linspace(-np.pi, np.pi, 361)

def best_r(phis):
    # pred(c0) = cos(c0)*U - sin(c0)*V; identical c0 treatment for real and null
    P = (C_FROZEN * np.sin(thk)[None, :] *
         np.exp(1j * (2 * ns[sel][:, None] * thk[None, :] - phis[None, :]))).sum(axis=1)
    U, V = P.real, P.imag
    preds = np.cos(c0s)[:, None] * U[None, :] - np.sin(c0s)[:, None] * V[None, :]
    pc = preds - preds.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    rs = pc @ yc / (np.linalg.norm(pc, axis=1) * np.linalg.norm(yc))
    i = int(np.argmax(rs))
    return float(rs[i]), float(c0s[i])

phi0 = S_FROZEN * psi
r_best, c0_best = best_r(phi0)
rng = np.random.default_rng(9)
null = np.array([best_r(rng.uniform(-np.pi, np.pi, len(gam)))[0] for _ in range(2000)])
z = (r_best - null.mean()) / null.std()
pemp = float(np.mean(null >= r_best))
print(f"\nfrozen-dictionary prediction: r = {r_best:.4f} at c0={c0_best:+.2f}")
print(f"random-phase null (same amp/freq structure, same c0 scan, n=2000): mean={null.mean():.4f} "
      f"sd={null.std():.4f} max={null.max():.4f} | z={z:+.2f} empirical p={pemp:.4f}", flush=True)
