"""S3015: RIGOROUS (Arb ball arithmetic) off-line-zero detector certificate (lane 2).
Builds tau = W02 - W_R - W_p rigorously by re-using the S2880 Arb builder (hp_arb_tau_certify.py, eshift=0),
forms the Loewner generator from the centre row g_j = pi*(nu_c - nu_j)*tau[c,j] (nu_n = 2 pi n / L, c <-> n=0),
removes the CERTIFIED zero contributions pi*sum_{0<gamma<=V+1} (4/L) sin^2(gamma L/2) (1/(nu_j-gamma) + 1/(nu_j+gamma))
(zeros via acb.zeta_zeros, completeness via arb.zeta_nzeros = certified N(V+1)),
and certifies every step s_j = -(gt_{j+1}-gt_j)/pi > 0 as a ball (lower endpoint > 0) => detector CERTIFIED QUIET.
usage: python3 arb_offline_detector_cert.py lam2 N prec"""
import sys, os
lam2s, Ns, precs = sys.argv[1], sys.argv[2], sys.argv[3]
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'hp_arb_tau_certify.py')).read()
src = src[:src.index('# RIGOROUS certificate')]
sys.argv = ['x', lam2s, Ns, precs, '0']
g_ = {}
exec(src, g_)
T, L, PI, N, M, arb, acb = g_['T'], g_['L'], g_['PI'], g_['N'], g_['M'], g_['arb'], g_['acb']
nu = [2 * PI * n / L for n in range(-N, N + 1)]
c = N
g = [arb(0) if j == c else PI * (nu[c] - nu[j]) * T[c, j] for j in range(M)]
# Loewner residual (ball must contain 0 for every off-diagonal pair if the single-Loewner structure holds)
resid_ok = all(((g[j] - g[i]) - PI * (nu[i] - nu[j]) * T[i, j]).contains(0) for i in range(M) for j in range(M) if i != j)
resid_max = max(float(abs(((g[j] - g[i]) - PI * (nu[i] - nu[j]) * T[i, j])).upper()) for i in range(M) for j in range(M) if i != j)
V = nu[-1]
Vp = V + 1
nz = int(arb.zeta_nzeros(Vp).mid().str(5, radius=False).split('.')[0].replace(' ', '')) if True else 0
nzb = arb.zeta_nzeros(Vp)
assert nzb.rad() == 0 or float(nzb.rad()) < 0.5, nzb
nz = int(round(float(nzb.mid())))
zs = [z.imag for z in acb.zeta_zeros(1, nz)] if nz > 0 else []
assert all(z < Vp for z in zs) and nz == len(zs)
def ztail(nj):
    s = arb(0)
    for z in zs:
        s += (4 / L) * (z * L / 2).sin() ** 2 * (1 / (nj - z) + 1 / (nj + z))
    return PI * s
gt = [g[j] - ztail(nu[j]) for j in range(M)]
steps = [-(gt[j + 1] - gt[j]) / PI for j in range(M - 1)]
cert = all(s.lower() > 0 for s in steps)
mnlow = min(float(s.lower()) for s in steps)
maxrad = max(float(s.rad()) for s in steps)
print(f"lam2={lam2s} N={N} prec={precs}: certified N(V+1)={nz} zeros removed (V+1={float(Vp):.4f}); "
      f"Loewner residual balls contain 0: {resid_ok} (max |resid| upper={resid_max:.2e}); "
      f"min step lower bound={mnlow:.6e} max step radius={maxrad:.2e}; DETECTOR CERTIFIED QUIET (all steps > 0): {cert}", flush=True)
# S3021 (additive): optional per-step dump for the nu-normalised statistic. Set DUMP=path.
if os.environ.get('DUMP'):
    with open(os.environ['DUMP'], 'w') as fh:
        fh.write('# j nu_j nu_{j+1} step_mid step_lower\n')
        for j, s in enumerate(steps):
            fh.write(f"{j} {float(nu[j].mid()):.10f} {float(nu[j+1].mid()):.10f} {float(s.mid()):.12e} {float(s.lower()):.12e}\n")
# S3035 (additive): certified centre-row PARTS (W02, W_R, W_P separately) + Fejer-smoothed decomposition. Set PARTS=1.
if os.environ.get('PARTS'):
    G = g_
    sh2, archconst, Jsig, Dd, omega, pw = G['sh2'], G['archconst'], G['Jsig'], G['Dd'], G['omega'], G['pw']
    idx = list(range(-N, N + 1))
    def w02(n, m): return 32*L*sh2*(L*L-16*PI*PI*m*n)/((L*L+16*PI*PI*m*m)*(L*L+16*PI*PI*n*n))
    def wr(n, m): return (archconst + Dd[abs(n)]) if n == m else (Jsig(m) - Jsig(n)) / (PI * (n - m))
    def wp(n, m): return sum((lk * arb(k) ** arb(-0.5) * omega(n, m, arb(k).log()) for k, lk in pw), arb(0))
    def genp(f): return [arb(0) if j == c else PI * (nu[c] - nu[j]) * f(0, idx[j]) for j in range(M)]
    g02, gR, gP = genp(w02), genp(wr), genp(wp)
    gZ = [ztail(nu[j]) for j in range(M)]
    stp = lambda g: [-(g[j + 1] - g[j]) / PI for j in range(M - 1)]
    S02, SR, SP, SZ = stp(g02), [-x for x in stp(gR)], [-x for x in stp(gP)], [-x for x in stp(gZ)]
    def sm(a, h, j):
        w = [arb(h + 1 - abs(i)) for i in range(-h, h + 1)]
        return sum((wi * a[j + i] for wi, i in zip(w, range(-h, h + 1))), arb(0)) / sum(w, arb(0))
    for h in [0, 8, 16, 24]:
        for name, j in [("centre", N), ("8th", N + N // 4)]:
            if j + h >= M - 1 or j - h < 0: continue
            a, b, p, z = sm(S02, h, j), sm(SR, h, j), sm(SP, h, j), sm(SZ, h, j)
            arith = a + b + p
            print(f"PARTS h={h:2d} {name:6s} s02={a.str(6)} -sR={b.str(6)} -sP={p.str(6)} arith={arith.str(6)} -sZ={z.str(6)} "
                  f"total={(arith + z).str(6)}  CERT arith>0:{arith > 0}  CERT explicit>|prime|:{(a + b) > abs(p)}", flush=True)
# S3038 (additive): ALL-grid-point certified smoothed dominance. Set PARTSALL="h1,h2,..." (requires PARTS=1).
if os.environ.get('PARTS') and os.environ.get('PARTSALL'):
    for h in [int(x) for x in os.environ['PARTSALL'].split(',')]:
        js = [j for j in range(h, M - 1 - h)]
        nA = nE = 0; worst = None
        for j in js:
            a, b, p = sm(S02, h, j), sm(SR, h, j), sm(SP, h, j)
            if (a + b + p) > 0: nA += 1
            if (a + b) > abs(p): nE += 1
            r = float(abs(p).upper()) / max(float((a + b).lower()), 1e-300)
            if worst is None or r > worst[0]: worst = (r, j)
        print(f"PARTSALL h={h:2d}: grid points {len(js)}: CERT arith>0 at {nA}/{len(js)}; CERT explicit>|prime| at {nE}/{len(js)}; "
              f"worst |prime|/explicit = {worst[0]:.3f} at j={worst[1]} (nu={float(nu[worst[1]].mid()):.2f})", flush=True)
