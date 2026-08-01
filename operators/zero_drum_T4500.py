import numpy as np, os
import mpmath as mp

# Big drum: T=3000 (~2470 zeros), K=2200 (valid window: K << 2N ~ 4930).
# Prime edges n* = T log(p)/pi: p=2:662, 3:1049, 4:1324, 5:1537, 7:1858, 8:1986,
# 9:2098, 11:2290(marginal). Upgraded filter: equilibrium kernel + low-t cleaning.

mp.mp.dps = 15
T, K = 4500.0, 3200
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros10k.npy"))
gam = np.array([g for g in Z if g < T])
N = len(gam)
print(f"T={T}: {N} zeros, K={K}", flush=True)

fn = os.path.join(DIR, f"weylpos_half_T{int(T)}.npy")
if os.path.exists(fn):
    tbar = np.load(fn)
else:
    tbar = np.empty(N)
    t = 14.0
    for k in range(1, N + 1):
        for _ in range(60):
            f = float(mp.siegeltheta(t)) / np.pi + 1 - (k - 0.5)
            fp = float(mp.log(t / (2 * np.pi))) / (2 * np.pi)
            step = f / fp
            t -= step
            if abs(step) < 1e-10:
                break
        tbar[k - 1] = t
        if k % 250 == 0:
            print(f"  weyl {k}/{N}", flush=True)
    np.save(fn, tbar)

def jac(pts, KK):
    x = np.concatenate([-pts[::-1], pts])
    mu = np.ones(len(x)) / len(x)
    v = np.sqrt(mu); v = v / np.linalg.norm(v)
    V = [v]; b = []
    w = x * v
    for k in range(KK):
        ak = v @ w
        w = w - ak * v - (b[-1] * V[-2] if b else 0)
        for u in V: w -= (u @ w) * u
        for u in V: w -= (u @ w) * u
        bk = np.linalg.norm(w); b.append(bk)
        vn = w / bk; V.append(vn); v = vn; w = x * v
        if k % 200 == 0:
            print(f"  lanczos {k}/{KK}", flush=True)
    return np.array(b)

db = jac(gam, K)[:-1] - jac(tbar, K)[:-1]
np.save(os.path.join(DIR, "zerodrum_db_T4500.npy"), db)
print(f"db rms={db.std():.4f}", flush=True)

ns = np.arange(len(db)); sel = ns >= 10
y = db[sel] - db[sel].mean(); nn = ns[sel].astype(float)
lowt = []
for g in gam[:30]:
    th = np.arccos(g / T)
    lowt += [np.cos(2 * nn * th), np.sin(2 * nn * th)]
Xl = np.column_stack(lowt)
yc = y - Xl @ np.linalg.lstsq(Xl, y, rcond=None)[0]
print(f"low-t projection removed {1 - np.var(yc) / np.var(y):.1%}", flush=True)
tq = np.linspace(0, T, 12001)[1:-1]
eqw = 1.0 / np.sqrt(T ** 2 - tq ** 2)
CM = np.cos(np.outer(2 * nn, np.arccos(tq / T)))
def score(om):
    t_ = CM @ (np.cos(om * tq) * eqw); t_ = t_ - t_.mean()
    t_ = t_ - Xl @ np.linalg.lstsq(Xl, t_, rcond=None)[0]
    return float(np.abs(t_ @ yc) / (np.linalg.norm(t_) * np.linalg.norm(yc)))

targets = [2, 3, 4, 5, 7, 8, 9, 11]
marks = np.array([np.log(t) for t in targets])
rng = np.random.default_rng(37)
scs = {}
for pm in targets:
    scs[pm] = score(np.log(pm))
    print(f"  p^m={pm:>2} n*={T*np.log(pm)/np.pi:5.0f}: |corr|={scs[pm]:.4f}", flush=True)
ctrl = []
for _ in range(200):
    while True:
        om = rng.uniform(0.4, 2.4)
        if np.min(np.abs(marks - om)) > 0.03: break
    ctrl.append(score(om))
ctrl = np.array(ctrl)
print(f"ctrl {ctrl.mean():.4f}±{ctrl.std():.4f} max {ctrl.max():.4f}")
print("z per target:", {pm: round((scs[pm] - ctrl.mean()) / ctrl.std(), 1) for pm in targets})
mp_ = np.mean(list(scs.values()))
print(f"prime-mean z = {(mp_ - ctrl.mean()) / (ctrl.std() / np.sqrt(len(targets))):+.2f}")
print("DONE", flush=True)
