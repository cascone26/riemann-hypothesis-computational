import numpy as np, os
from krein_fit import lanczos_jacobi, design

# Phase spectroscopy: extract per-zero phases phi_k from the pinned-frequency fit
# db_n ~ sum_k w_k cos(2*theta_k*n - phi_k). If phi_k is LINEAR in theta_k
# (phi = c1*theta + c0), then db_n ~ sum w_k cos(theta_k*(2n - c1) - c0): the
# prime comb survives with n shifted by c1/2 — rerun the comb at the corrected
# offset. If phases are unstructured, the comb is genuinely washed out.

DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
RUNS = [(30.0, 8000, 70, "A30_M8000"), (40.0, 10000, 95, "A40_M10000"),
        (60.0, 14000, 140, "A60_M14000"), (100.0, 20000, 220, "A100_M20000")]

def primepowers(lo, hi):
    def isp(q): return q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))
    out = []
    for p in [q for q in range(2, 200) if isp(q)]:
        pm = p
        while np.log(pm) < hi:
            if np.log(pm) > lo: out.append(pm)
            pm *= p
    return sorted(set(out))

rng = np.random.default_rng(5)
combz = []
for A, M, K, tag in RUNS:
    th = (np.arange(M) + 0.5) * np.pi / M
    ts = A * np.cos(th)
    quad = A * np.sin(th) * (np.pi / M)
    wz = np.load(os.path.join(DIR, f"zw_{tag}.npy"))
    _, bz = lanczos_jacobi(ts, wz * quad, K)
    _, bf = lanczos_jacobi(ts, quad, K)
    db = bz[:-1] - bf[:-1]
    gam = [float(g) for g in Z if g < A - 0.5]
    thk = np.array([np.arccos(g / A) for g in gam])
    ns = np.arange(1, len(db))[4:]
    y = db[5:]
    X = design(2 * thk, ns, False)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    cc = beta[1::2]; ss = beta[2::2]
    phi = np.arctan2(ss, cc)          # db ~ sum |w| cos(2 th n - phi)
    # circular-linear fit phi ~ c1*theta + c0: scan c1, maximize resultant length
    c1s = np.linspace(-8, 8, 3201)
    Rs = np.array([np.abs(np.mean(np.exp(1j * (phi - c1 * thk)))) for c1 in c1s])
    i = int(np.argmax(Rs)); c1 = float(c1s[i]); R = float(Rs[i])
    c0 = float(np.angle(np.mean(np.exp(1j * (phi - c1 * thk)))))
    print(f"A={A:5.0f}: {len(gam):2d} zeros | phase fit phi ~ {c1:+.3f}*theta {c0:+.3f} | "
          f"resultant R={R:.3f} (R~1 tight, R~0 unstructured)")
    for g, t_, p_ in list(zip(gam, thk, phi))[:6]:
        print(f"    g={g:7.3f} theta={t_:.4f} phi={p_:+.4f} pred={(c1 * t_ + c0 + np.pi) % (2 * np.pi) - np.pi:+.4f}")
    # corrected comb: db ~ sum w cos(theta*(2n - c1) - c0) => spikes where
    # (2n - c1)/A = log(p^m) i.e. n = (A*log(p^m) + c1)/2
    x = np.abs(db)
    med = np.array([np.median(x[max(0, i2 - 3):i2 + 4]) for i2 in range(len(x))])
    ex = x - med
    nlo, nhi = 8, len(ex) - 3
    pms = primepowers((2 * nlo - c1) / A, (2 * nhi - c1) / A)
    marks = np.array([(A * np.log(pm) + c1) / 2 for pm in pms])
    marks = marks[(marks > nlo) & (marks < nhi)]
    sc = lambda mk: float(np.mean(np.interp(mk, np.arange(len(ex)), ex)))
    s0 = sc(marks)
    span = nhi - nlo
    ctrl = np.array([sc(nlo + (marks - nlo + rng.uniform(2, span - 2)) % span) for _ in range(2000)])
    zz = (s0 - ctrl.mean()) / ctrl.std()
    combz.append(zz)
    print(f"    corrected comb ({len(marks)} marks, offset {c1 / 2:+.2f}): score={s0:+.4f} "
          f"ctrl={ctrl.mean():+.4f}±{ctrl.std():.4f} z={zz:+.2f}")
print(f"\ncombined corrected-comb Stouffer z = {sum(combz) / np.sqrt(len(combz)):+.2f}")
