import numpy as np, os
from krein_fit import lanczos_jacobi

# TWO tests on the Krein string, all four cached windows A=30/40/60/100.
#
# 1) PRIME COMB (corrected scale): in the A->inf limit 2n*arccos(g/A) ~ n*pi - 2ng/A,
#    so db_n ~ (-1)^n * sum_k w_k cos(2*(n/A)*g_k + ph) — the explicit-formula kernel
#    at u = 2n/A, which spikes at u = log(p^m). Prediction: |db| bumps at
#    n ~ (A/2)*log(p^m). Statistic: mean excess at the comb marks, p-value from
#    circular shifts of the comb.
# 2) AMPLITUDE LAW: amp(g,A) = C*sin(arccos(g/A)) with ONE constant C across all
#    windows — no gamma-dependence beyond geometry <=> equal weight (all zeros
#    simple) in the string.

DIR = os.path.dirname(os.path.abspath(__file__))
RUNS = [(30.0, 8000, 70, "A30_M8000"), (40.0, 10000, 95, "A40_M10000"),
        (60.0, 14000, 140, "A60_M14000"), (100.0, 20000, 220, "A100_M20000")]

def primepowers(lo, hi):
    def isp(q):
        return q > 1 and all(q % d for d in range(2, int(q ** 0.5) + 1))
    out = []
    for p in [q for q in range(2, 200) if isp(q)]:
        pm = p
        while np.log(pm) < hi:
            if np.log(pm) > lo:
                out.append(pm)
            pm *= p
    return sorted(set(out))

rng = np.random.default_rng(3)
print("== 1) PRIME COMB at n ~ (A/2)*log(p^m) ==")
allscores = []
for A, M, K, tag in RUNS:
    th = (np.arange(M) + 0.5) * np.pi / M
    ts = A * np.cos(th)
    quad = A * np.sin(th) * (np.pi / M)
    wz = np.load(os.path.join(DIR, f"zw_{tag}.npy"))
    _, bz = lanczos_jacobi(ts, wz * quad, K)
    _, bf = lanczos_jacobi(ts, quad, K)
    db = bz[:-1] - bf[:-1]
    x = np.abs(db)
    med = np.array([np.median(x[max(0, i - 3):i + 4]) for i in range(len(x))])
    ex = x - med
    nlo, nhi = 8, len(ex) - 3
    pms = primepowers(2 * nlo / A, 2 * nhi / A)
    marks = np.array([(A / 2) * np.log(pm) for pm in pms])
    def score(mk):
        return float(np.mean(np.interp(mk, np.arange(len(ex)), ex)))
    s0 = score(marks)
    span = nhi - nlo
    ctrl = np.array([score(nlo + (marks - nlo + rng.uniform(2, span - 2)) % span)
                     for _ in range(2000)])
    p = float(np.mean(ctrl >= s0))
    z = (s0 - ctrl.mean()) / ctrl.std()
    allscores.append((s0 - ctrl.mean()) / ctrl.std())
    print(f"  A={A:5.0f}: {len(marks):3d} marks in [{nlo},{nhi}] | score={s0:+.4f} "
          f"ctrl={ctrl.mean():+.4f}±{ctrl.std():.4f} | z={z:+.2f} p={p:.4f}")
print(f"  combined Stouffer z = {sum(allscores) / np.sqrt(len(allscores)):+.2f}")

print("\n== 2) AMPLITUDE LAW amp = C*sqrt(1-(g/A)^2) ==")
AMPS = {
    30.0: [(14.1347, 0.719), (21.0220, 0.484), (25.0109, 0.398)],
    40.0: [(14.1347, 0.679), (21.0220, 0.580), (25.0109, 0.527), (30.4249, 0.527),
           (32.9351, 0.473), (37.5862, 0.314)],
    60.0: [(14.1347, 0.6320), (21.0220, 0.6841), (25.0109, 0.6127), (30.4249, 0.6259),
           (32.9351, 0.6253), (37.5862, 0.4140), (40.9187, 0.4356), (43.3271, 0.4756),
           (48.0052, 0.4986), (49.7738, 0.4183), (52.9703, 0.3900), (56.4462, 0.2390),
           (59.3470, 0.1638)],
    100.0: [(14.1347, 0.6999), (21.0220, 0.6243), (25.0109, 0.5627), (30.4249, 0.6972),
            (32.9351, 0.6283), (37.5862, 0.5146), (40.9187, 0.5937), (43.3271, 0.5779),
            (48.0052, 0.6066), (49.7738, 0.5314), (52.9703, 0.5309), (56.4462, 0.4228),
            (59.3470, 0.4825), (60.8318, 0.5938), (65.1125, 0.5203), (67.0798, 0.3982),
            (69.5464, 0.3996), (72.0672, 0.3820), (75.7047, 0.5082), (77.1448, 0.3970),
            (79.3374, 0.4083), (82.9104, 0.3892), (84.7355, 0.3654), (87.4253, 0.4126),
            (88.8091, 0.4545), (92.4919, 0.3201), (94.6513, 0.2875), (95.8706, 0.2773),
            (98.8312, 0.1622)],
}
xs, ys = [], []
for A, rows in AMPS.items():
    for g, a in rows:
        xs.append(np.sqrt(1 - (g / A) ** 2)); ys.append(a)
xs = np.array(xs); ys = np.array(ys)
C = float(xs @ ys / (xs @ xs))
res = ys - C * xs
R2 = 1 - np.sum(res ** 2) / np.sum((ys - ys.mean()) ** 2)
print(f"  {len(ys)} (gamma,A) points | best C = {C:.4f} | R^2 through model = {R2:.4f}")
print(f"  residual sd = {res.std():.4f} vs amp sd = {ys.std():.4f}")
corr = np.corrcoef(xs, ys)[0, 1]
print(f"  corr(amp, sin_theta) = {corr:.4f}")
print("  law holding with ONE C across windows => equal per-zero weight (simplicity) in the string.")
