import numpy as np, os
import mpmath as mp
from krein_fit import lanczos_jacobi, design, r2
import krein_fit

# A->inf ladder rung: A=60 (13 zeros inside), K=140. Do the pinned zero-frequency
# amplitudes persist as the window grows? (continuum limit toward the de Branges
# canonical system of zeta)

mp.mp.dps = 15
A, M, K = 100.0, 20000, 220
Z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "zeros500.npy"))
gam = [float(g) for g in Z if g < A - 0.5]
print(f"A={A}: {len(gam)} zeros in window: {[round(g,2) for g in gam]}", flush=True)

th = (np.arange(M) + 0.5) * np.pi / M
ts = A * np.cos(th)
quad = A * np.sin(th) * (np.pi / M)
fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"zw_A100_M{M}.npy")
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
db = bz[:-1] - bf[:-1]
ns = np.arange(1, len(db))[4:]
y = db[5:]
freqs = [2 * np.arccos(g / A) for g in gam]
nf = len(freqs)
rng = np.random.default_rng(11)
R, beta = r2(freqs, ns, y, False)
ctrl = np.array([r2(rng.uniform(0.1, np.pi - 0.1, nf), ns, y, False)[0] for _ in range(200)])
z = (R - ctrl.mean()) / ctrl.std()
amps = [float(np.hypot(beta[1 + 2 * i], beta[2 + 2 * i])) for i in range(nf)]
print(f"\nR^2(zeros)={R:.4f} | ctrl mean={ctrl.mean():.4f} sd={ctrl.std():.4f} max={ctrl.max():.4f} | z={z:+.2f}", flush=True)
print("per-zero amplitudes:")
for g, a in zip(gam, amps):
    print(f"  gamma={g:8.4f}  amp={a:.4f}")
print("\ncompare vs A=30/40/60 ladders:")
print("persisting amplitudes as A grows => the infinite string carries ALL zeros.", flush=True)
