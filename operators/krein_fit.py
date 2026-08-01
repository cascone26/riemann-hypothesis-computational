import numpy as np, os
rng = np.random.default_rng(7)

# Inverse-spectral test: db_n (Jacobi off-diagonal imprint of the |zeta|^2 measure)
# should be a superposition of oscillations at FIXED frequencies 2*arccos(gamma_k/A),
# one per zeta zero gamma_k < A (interior double zero of the weight -> generalized
# Jacobi oscillation). Frequencies are pinned by known zeros; only amplitudes/phases
# are fit. Control: same fit with random frequency sets.

GAMMAS = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062, 37.586178]

def lanczos_jacobi(x, mu, K):
    v = np.sqrt(mu); v = v / np.linalg.norm(v)
    V = [v]; a = []; b = []
    w = x * v
    for k in range(K):
        ak = v @ w; a.append(ak)
        w = w - ak * v - (b[-1] * V[-2] if b else 0)
        for u in V: w -= (u @ w) * u
        for u in V: w -= (u @ w) * u
        bk = np.linalg.norm(w); b.append(bk)
        vn = w / bk; V.append(vn); v = vn; w = x * v
    return np.array(a), np.array(b)

def get_db(A, M, K, tag):
    th = (np.arange(M) + 0.5) * np.pi / M
    ts = A * np.cos(th)
    quad = A * np.sin(th) * (np.pi / M)
    wz = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"zw_{tag}.npy"))
    _, bz = lanczos_jacobi(ts, wz * quad, K)
    _, bf = lanczos_jacobi(ts, quad, K)
    return bz[:-1] - bf[:-1]

def design(freqs, ns, envelope):
    cols = [np.ones_like(ns, dtype=float)]
    for f in freqs:
        env = (ns[0] / ns) if envelope else 1.0
        cols.append(np.cos(f * ns) * env)
        cols.append(np.sin(f * ns) * env)
    return np.column_stack(cols)

def r2(freqs, ns, y, envelope):
    X = design(freqs, ns, envelope)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    return 1 - np.sum(res ** 2) / np.sum((y - y.mean()) ** 2), beta

if __name__ != "__main__":
    pass
_RUNS = [(30.0, 8000, 70, "A30_M8000"), (40.0, 10000, 95, "A40_M10000")] if __name__ == "__main__" else []
for A, M, K, tag in _RUNS:
    db = get_db(A, M, K, tag)
    ns = np.arange(1, len(db))[4:]          # drop n<=4 edge transients
    y = db[5:]
    gam = [g for g in GAMMAS if g < A]
    freqs = [2 * np.arccos(g / A) for g in gam]
    nf = len(freqs)
    print(f"\n=== A={A}: {nf} zeros below A, {len(y)} samples, {1+2*nf} fit params ===")
    print("  pinned frequencies 2*arccos(g/A):", [round(f, 4) for f in freqs])
    for env in (False, True):
        R, beta = r2(freqs, ns, y, env)
        amps = [np.hypot(beta[1 + 2 * i], beta[2 + 2 * i]) for i in range(nf)]
        # control: random frequency sets, same count, uniform in (0.1, pi-0.1)
        ctrl = [r2(rng.uniform(0.1, np.pi - 0.1, nf), ns, y, env)[0] for _ in range(200)]
        ctrl = np.array(ctrl)
        z = (R - ctrl.mean()) / ctrl.std()
        print(f"  envelope={'1/n' if env else 'none'}: R^2(zeros)={R:.4f} | "
              f"random-freq R^2 mean={ctrl.mean():.4f} sd={ctrl.std():.4f} max={ctrl.max():.4f} | z={z:+.2f}")
        if not env:
            print("    per-zero amplitudes:", [f"{g:.2f}:{a:.3f}" for g, a in zip(gam, amps)])
print("\nverdict: z >> 3 with pinned zero-frequencies beating the random-control max => "
      "the zeta zeros are READABLE from the Jacobi coefficients (inverse-spectral imprint).")
