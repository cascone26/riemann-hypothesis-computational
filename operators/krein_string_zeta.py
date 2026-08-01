import numpy as np, os, sys
import mpmath as mp

# Krein-string / Szego test: does the Jacobi matrix of the zeta spectral measure
#   d mu(t) = |zeta(1/2+it)|^2 dt  on [-A, A]
# carry the PRIMES in its off-diagonal coefficients?
# Prediction (Szego/Verblunsky heuristics): under t = A cos(theta), a t-frequency
# log(p^m) in log|zeta|^2 lands at Jacobi index n ~ A*log(p^m). The killer control:
# rerun at a second A and check the bumps TRACK A*log p.

mp.mp.dps = 15

def zeta_weight(ts, tag):
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), f"zw_{tag}.npy")
    if os.path.exists(fn):
        w = np.load(fn)
        if len(w) == len(ts):
            return w
    w = np.empty(len(ts))
    for i, t in enumerate(ts):
        z = mp.zeta(mp.mpc(0.5, t))
        w[i] = float(abs(z)) ** 2
        if i % 1000 == 0:
            print(f"  [{tag}] zeta grid {i}/{len(ts)}", flush=True)
    np.save(fn, w)
    return w

def lanczos_jacobi(x, mu, K):
    # discretized-Stieltjes: Lanczos with full reorthogonalization on diag(x), weights mu
    v = np.sqrt(mu); v = v / np.linalg.norm(v)
    V = [v]
    a = []; b = []
    w = x * v
    for k in range(K):
        ak = v @ w
        a.append(ak)
        w = w - ak * v - (b[-1] * V[-2] if b else 0)
        for u in V:            # full reorthogonalization, twice
            w -= (u @ w) * u
        for u in V:
            w -= (u @ w) * u
        bk = np.linalg.norm(w)
        b.append(bk)
        vn = w / bk
        V.append(vn)
        v = vn
        w = x * v
    return np.array(a), np.array(b)

def run(A, M, K):
    print(f"\n=== A={A}  M={M} nodes  K={K} Jacobi coefficients ===", flush=True)
    th = (np.arange(M) + 0.5) * np.pi / M
    ts = A * np.cos(th)
    quad = A * np.sin(th) * (np.pi / M)        # dt weight for theta-uniform nodes
    wz = zeta_weight(ts, f"A{int(A)}_M{M}")
    mu_z = wz * quad
    mu_flat = quad                              # w==1 baseline -> semicircle, b_n -> A/2
    a_z, b_z = lanczos_jacobi(ts, mu_z, K)
    a_f, b_f = lanczos_jacobi(ts, mu_flat, K)
    db = b_z[:-1] - b_f[:-1]
    # Fourier side: cosine coefficients of log-weight in theta
    L = np.log(np.maximum(wz, 1e-300))
    ks = np.arange(K + 5)
    ch = np.array([(2.0 / M) * np.sum(L * np.cos(k * th)) for k in ks])
    marks = []
    for p in [2, 3, 5, 7, 11, 13]:
        pm = p
        while A * np.log(pm) < K:
            marks.append((pm, A * np.log(pm)))
            pm *= p
    marks.sort(key=lambda z: z[1])
    print("predicted prime-power bump locations n ~ A*log(p^m):", flush=True)
    print("  " + "  ".join(f"{pm}:{loc:.1f}" for pm, loc in marks), flush=True)
    print(f"\n{'n':>4} {'b_zeta':>10} {'b_flat':>10} {'db':>9} {'chat(logw)':>11}", flush=True)
    for n in range(1, len(db)):
        near = "".join(f" <-- {pm}^ ~A*log" for pm, loc in marks if abs(n - loc) < 0.55)
        print(f"{n:>4} {b_z[n]:>10.4f} {b_f[n]:>10.4f} {db[n]:>9.4f} {ch[n]:>11.4f}{near}", flush=True)
    # bump detection on |db| after removing smooth trend (centered moving median, w=7)
    A_ = np.abs(db[1:])
    med = np.array([np.median(A_[max(0, i - 3):i + 4]) for i in range(len(A_))])
    ex = A_ - med
    top = np.argsort(ex)[::-1][:12] + 1
    print(f"\ntop-12 |db| excess bumps at n = {sorted(top.tolist())}", flush=True)
    print(f"prime-power grid at A={A}: {[round(l, 1) for _, l in marks]}", flush=True)
    hits = sum(1 for n in top if any(abs(n - loc) <= 1.5 for _, loc in marks))
    print(f"bumps within 1.5 of an A*log(p^m) mark: {hits}/12", flush=True)
    return db, ch, marks

if __name__ == "__main__":
    db1, ch1, m1 = run(A=30.0, M=8000, K=70)
    db2, ch2, m2 = run(A=40.0, M=10000, K=95)
    print("\nCONTROL: if the bump family moved from the A=30 grid to the A=40 grid,")
    print("the primes are IN the operator coefficients (Krein string of zeta).", flush=True)
