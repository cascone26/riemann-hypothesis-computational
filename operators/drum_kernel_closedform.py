import numpy as np, os

# DIRECTION 2: drum kernel closed form via equilibrium-measure calculus.
# For a discrete measure of N points on [-T,T] with local density rho(t) (~Weyl),
# the equilibrium (arcsine-type) map is theta(t)=arccos(t/T). Perturbing one point
# x_j by dx generates, to leading order, the response in Jacobi off-diag:
#   R_n(x_j) = d b_n / d x_j.
# HYPOTHESIS (derivative-Chebyshev / equilibrium): the string coefficient b_n couples
# to the measure via the n-th Chebyshev moment; displacing a point moves that moment by
#   d/dx [ w_j T_n(x/T) ] with w_j the (equal) mass. In the orthonormal frame this gives
#   R_n(x) = (C0/N) * (n / (T sin theta))^? * sin(2 n theta) * <lattice form factor>.
# Rather than guess the power, EXTRACT it cleanly from the CENTER (theta~pi/2) where the
# single-frequency fit was R^2=1.000, across a ladder of T so the T-scaling is pinned.

def jacobis(pts, K):
    x = np.concatenate([-pts[::-1], pts]); mu = np.ones(len(x)) / len(x)
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

DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros10k.npy"))
print("CENTER kernel amplitude & T-scaling (theta=pi/2, x_j~0): A_center vs T, N")
print(f"{'T':>6} {'N':>5} {'K':>4} {'A_center':>12} {'A*N':>10} {'A*N*T':>10}")
rows = []
for T in [500.0, 1000.0, 2000.0]:
    gam = np.array([g for g in Z if g < T])
    N = len(gam)
    # smooth control = Weyl positions (k-1/2)
    import mpmath as mp
    mp.mp.dps = 12
    tb = np.empty(N); t = 14.0
    for k in range(1, N + 1):
        for _ in range(50):
            f = float(mp.siegeltheta(t)) / np.pi + 1 - (k - 0.5)
            fp = float(mp.log(t / (2 * np.pi))) / (2 * np.pi)
            t -= f / fp
            if abs(f / fp) < 1e-9: break
        tb[k - 1] = t
    K = min(int(0.45 * 2 * N), 400)
    b0 = jacobis(tb, K)
    # displace the point nearest 0 (there is none; smallest is ~14) -> use central index in x-order
    j = 0  # smallest positive point ~ gamma_1; theta near pi/2 only if x<<T
    # instead measure amplitude of response to the SMALLEST point (closest to center of [-T,T])
    eps = 1e-3
    tp = tb.copy(); tp[j] += eps
    R = (jacobis(tp, K) - b0)[:-1] / eps
    th = np.arccos(tb[j] / T)
    n = np.arange(len(R)).astype(float); sel = (n >= 20) & (n <= 0.7 * K)
    X = np.column_stack([n[sel] * np.cos(2 * th * n[sel]), n[sel] * np.sin(2 * th * n[sel])])
    beta, *_ = np.linalg.lstsq(X, R[sel], rcond=None)
    A = float(np.hypot(*beta)); R2 = 1 - np.var(R[sel] - X @ beta) / np.var(R[sel])
    rows.append((T, N, A))
    print(f"{T:>6.0f} {N:>5} {K:>4} {A:>12.3e} {A*N:>10.4f} {A*N*T:>10.3f}  (R2={R2:.3f}, theta={th:.3f})", flush=True)
print("\nif A*N ~ const: kernel amp = C/N (pure mass factor). if A*N*T ~ const: C/(N T).")
print("DONE", flush=True)
