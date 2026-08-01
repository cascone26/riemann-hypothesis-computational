import numpy as np, os
import mpmath as mp
from scipy.fft import dct

# A=200 string (79 zeros): full parameter-free dictionary verification at 2.7x scale.
# (i) build weight grid, (ii) Lanczos K=440, (iii) pinned-frequency fit,
# (iv) phases vs Vanlessen (1.6) [counting + slope + DCT conj of log h], parameter-free,
# (v) amplitudes vs theory scale (late half), (vi) frequency z-score vs random controls.

mp.mp.dps = 15
A, M, K = 200.0, 30000, 440
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros10k.npy"))
gam = [float(g) for g in Z if g < A - 0.5]
print(f"A={A}: {len(gam)} zeros, M={M}, K={K}", flush=True)

th = (np.arange(M) + 0.5) * np.pi / M
ts = A * np.cos(th)
x = np.cos(th)
quad = A * np.sin(th) * (np.pi / M)
fn = os.path.join(DIR, f"zw_A200_M{M}.npy")
if os.path.exists(fn):
    wz = np.load(fn)
else:
    wz = np.empty(M)
    for i, t in enumerate(ts):
        wz[i] = float(abs(mp.zeta(mp.mpc(0.5, t)))) ** 2
        if i % 2000 == 0:
            print(f"  zeta grid {i}/{M}", flush=True)
    np.save(fn, wz)

def lanczos(mu, KK):
    v = np.sqrt(mu); v = v / np.linalg.norm(v)
    V = [v]; b = []
    w = ts * v
    for k in range(KK):
        ak = v @ w
        w = w - ak * v - (b[-1] * V[-2] if b else 0)
        for u in V: w -= (u @ w) * u
        for u in V: w -= (u @ w) * u
        bk = np.linalg.norm(w); b.append(bk)
        vn = w / bk; V.append(vn); v = vn; w = ts * v
        if k % 50 == 0:
            print(f"  lanczos {k}/{KK}", flush=True)
    return np.array(b)

bz = lanczos(wz * quad, K)
bf = lanczos(quad, K)
db = bz[:-1] - bf[:-1]
np.save(os.path.join(DIR, "db_A200.npy"), db)

thk = np.array([np.arccos(g / A) for g in gam])
idx = np.arange(len(db)); sel = (idx >= K // 2) & (idx <= K - 2)
y = db[sel]; n = (idx[sel] + 1).astype(float)
cols = [np.ones(len(n))]
for t_ in thk:
    amp = (A / n) * np.sin(t_)
    cols += [amp * np.cos(2 * n * t_), amp * np.sin(2 * n * t_)]
X = np.column_stack(cols)
beta, *_ = np.linalg.lstsq(X, y, rcond=None)
res = y - X @ beta
R2 = 1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2)
scales = np.array([np.hypot(beta[1 + 2 * i], beta[2 + 2 * i]) for i in range(len(gam))])
print(f"\nlate-half dictionary fit: R2={R2:.4f}", flush=True)
print(f"per-zero scales: mean={scales.mean():.3f} sd={scales.std():.3f} min={scales.min():.3f} max={scales.max():.3f} (theory->1, finite-n below)")
print(f"simplicity check: any scale > 1.5x median? {np.sum(scales > 1.5 * np.median(scales))} zeros (multiplicity-2 would read ~2x)")

# parameter-free phase closure vs Vanlessen (1.6)
xs = np.array(sorted([g / A for g in gam] + [-g / A for g in gam]))
n0 = len(xs)
logh = np.log(np.maximum(wz, 1e-300))
for s in xs: logh -= 2 * np.log(np.abs(x - s))
a = dct(logh, type=2, norm=None) / M
KF = min(6000, M - 1); ks = np.arange(1, KF)
def conj_h(t_): return float(np.sum(a[1:KF] * np.sin(ks * t_)))
def Phi(xv):
    t_ = np.arccos(xv)
    above = int(np.sum(xs > xv + 1e-12))
    return (1 + 2 * above) * np.pi - (2 * n0) * t_ - conj_h(t_)
resid = []
for i, g in enumerate(gam):
    phi_fit = np.arctan2(beta[2 + 2 * i], beta[1 + 2 * i])
    eff = (Phi(g / A) - Phi(-g / A)) / 2
    resid.append(np.angle(np.exp(1j * (phi_fit - eff - np.pi))))
resid = np.array(resid)
R = float(np.abs(np.mean(np.exp(1j * resid)))); mn = float(np.angle(np.mean(np.exp(1j * resid))))
print(f"\nPHASE CLOSURE (parameter-free): resultant R={R:.3f} mean={mn:+.3f} sd={resid.std():.3f} over {len(gam)} zeros")
print("DONE", flush=True)
