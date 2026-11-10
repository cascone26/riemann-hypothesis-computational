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
