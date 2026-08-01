import numpy as np, os
import mpmath as mp

# THE ZERO-DRUM: Krein/Jacobi string of the discrete measure sum_k delta_{+-gamma_k}
# (unit masses) vs the smooth Weyl-density control (positions t_k solving N_bar = k).
# This is the explicit finite section of the (Suzuki-abstract) Hilbert-Polya string.
# Duality prediction: its coefficient deviations db_n carry the PRIMES.

mp.mp.dps = 15
T = 1000.0
K = 500
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros10k.npy"))
gam = np.array([g for g in Z if g < T])
N = len(gam)
print(f"T={T}: {N} zeros -> {2*N} support points, K={K}", flush=True)

# smooth control positions: Nbar(t) = theta_H(t)/pi + 1 ; solve Nbar(tbar_k) = k
fn = os.path.join(DIR, f"weylpos_half_T{int(T)}.npy")
if os.path.exists(fn):
    tbar = np.load(fn)
else:
    tbar = np.empty(N)
    t = 14.0
    for k in range(1, N + 1):
        # Newton on f(t) = theta(t)/pi + 1 - k
        for _ in range(60):
            f = float(mp.siegeltheta(t)) / np.pi + 1 - (k - 0.5)
            fp = float(mp.log(t / (2 * np.pi))) / (2 * np.pi)
            step = f / fp
            t -= step
            if abs(step) < 1e-10:
                break
        tbar[k - 1] = t
        if k % 100 == 0:
            print(f"  weyl {k}/{N}", flush=True)
    np.save(fn, tbar)
print(f"weyl positions: first {tbar[0]:.3f} (gamma_1={gam[0]:.3f}), last {tbar[-1]:.3f} (gamma_N={gam[-1]:.3f})", flush=True)

def jacobi_of_points(pts, K):
    x = np.concatenate([-pts[::-1], pts])
    mu = np.ones(len(x)) / len(x)
    v = np.sqrt(mu); v = v / np.linalg.norm(v)
    V = [v]; b = []
    w = x * v
    for k in range(K):
        ak = v @ w
        w = w - ak * v - (b[-1] * V[-2] if b else 0)
        for u in V: w -= (u @ w) * u
        for u in V: w -= (u @ w) * u
        bk = np.linalg.norm(w); b.append(bk)
        vn = w / bk; V.append(vn); v = vn; w = x * v
    return np.array(b)

bz = jacobi_of_points(gam, K)
bs = jacobi_of_points(tbar, K)
db = bz[:-1] - bs[:-1]
np.save(os.path.join(DIR, "zerodrum_db2.npy"), db)
print(f"db: rms={np.sqrt(np.mean(db**2)):.4f} max|.|={np.abs(db).max():.4f}", flush=True)

# Fourier analysis of db_n: peaks vs candidate prime-comb scalings
ns = np.arange(len(db))
sel = ns >= 10
y = db[sel] - db[sel].mean()
n = ns[sel].astype(float)
# scan oscillation frequency Omega in db ~ cos(Omega n): power spectrum
Oms = np.linspace(0.005, np.pi, 6000)
P = np.array([np.abs(np.sum(y * np.exp(1j * Om * n))) for Om in Oms])
from scipy.signal import argrelmax
ix = argrelmax(P, order=15)[0]
top = ix[np.argsort(P[ix])[::-1][:12]]
print("\ntop-12 spectral peaks of db_n (frequency Omega, power):", flush=True)
def isp(q): return q > 1 and all(q % d for d in range(2, int(q**0.5) + 1))
pps = []
for p in [q for q in range(2, 60) if isp(q)]:
    pm = p
    while pm < 60:
        pps.append(pm); pm *= p
pps = sorted(set(pps))
for i in np.sort(top):
    Om = Oms[i]
    # candidate identifications of Omega with log p under scalings
    cands = []
    for pm in pps:
        for name, val in [("logp/pi*?", np.log(pm))]:
            pass
    print(f"  Omega={Om:.4f}  power={P[i]:.3f}  | Omega*K/pi={Om*K/np.pi:.2f}  T*?: {Om:.4f}")
print("\n(identification pass follows — compare Omega against log(p^m)/scalings in next step)", flush=True)
print("DONE", flush=True)
