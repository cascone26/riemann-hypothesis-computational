import numpy as np, os

# Does the Krein string's off-line detection stiffness GROW with string length K?
# If deficit(delta, K) scales like (K*delta)^2 before saturating, the infinite
# string detects any delta>0 — operator-side "RH detectable in the limit".

A, M = 30.0, 8000
G1 = 14.134725
GAMMAS = [14.134725, 21.022040, 25.010858]
th = (np.arange(M) + 0.5) * np.pi / M
ts = A * np.cos(th)
quad = A * np.sin(th) * (np.pi / M)
wz = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"zw_A{int(A)}_M{M}.npy"))
freqs = [2 * np.arccos(g / A) for g in GAMMAS]

def lanczos(mu, K):
    v = np.sqrt(mu); v = v / np.linalg.norm(v)
    V = [v]; b = []
    w = ts * v
    for k in range(K):
        ak = v @ w
        w = w - ak * v - (b[-1] * V[-2] if b else 0)
        for u in V: w -= (u @ w) * u
        for u in V: w -= (u @ w) * u
        bk = np.linalg.norm(w); b.append(bk)
        vn = w / bk; V.append(vn); v = vn; w = ts * v
    return np.array(b)

def lift(delta):
    f = np.ones_like(ts)
    for g in (G1, -G1):
        f *= ((ts - g) ** 2 + delta ** 2) / np.maximum((ts - g) ** 2, 1e-18)
    return wz * f

def amp(w, K, bf):
    bz = lanczos(w * quad, K)
    y = (bz[:-1] - bf[:-1])[5:]
    ns = np.arange(1, K - 1)[4:]
    cols = [np.ones_like(ns, dtype=float)]
    for f in freqs:
        cols.append(np.cos(f * ns)); cols.append(np.sin(f * ns))
    X = np.column_stack(cols)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return np.hypot(beta[1], beta[2])

print(f"{'K':>5} {'amp0':>8} " + " ".join(f"def(d={d})" for d in [0.003, 0.01, 0.03]), flush=True)
for K in [70, 150, 300, 600]:
    bf = lanczos(quad, K)
    a0 = amp(wz, K, bf)
    defs = []
    for d in [0.003, 0.01, 0.03]:
        a = amp(lift(d), K, bf)
        defs.append(1 - a / a0)
    print(f"{K:>5} {a0:>8.4f} " + " ".join(f"{x:>9.2e}" for x in defs), flush=True)
print("\nif deficit ~ K^2*delta^2: each K doubling -> ~4x deficit. saturation K ~ 1/delta.", flush=True)
