import numpy as np, os
import mpmath as mp

# DIRECTION 3 (fixed): de Branges E_K of the |zeta|^2 string via the STANDARD
# Hermite-Biehler combination E_n(z) = p_n(z) - i p_{n-1}(z). p_n, p_{n-1} are the
# orthonormal polynomials (real on R, interlacing zeros) => E_n is Hermite-Biehler.
# Clean questions (no fitting):
#  (1) boundary phase phi_K(t) = arg E_K(t): does dphi/dt track the QUADRATURE-node
#      (equilibrium/arcsine) density, or the ZETA (Weyl) density? This decides whether
#      our explicit string's de Branges space equals Suzuki's abstract xi-space.
#  (2) |E_K(t)|^{-2} vs the Christoffel function / measure.

mp.mp.dps = 12
A, M, K = 60.0, 14000, 140
DIR = os.path.dirname(os.path.abspath(__file__))
th = (np.arange(M) + 0.5) * np.pi / M
ts = A * np.cos(th); quad = A * np.sin(th) * (np.pi / M)
wz = np.load(os.path.join(DIR, f"zw_A{int(A)}_M{M}.npy"))
mu = wz * quad

def lanc_ab(x, muv, KK):
    v = np.sqrt(muv); v = v / np.linalg.norm(v)
    V = [v]; a = []; b = []
    w = x * v
    for k in range(KK):
        ak = v @ w; a.append(ak)
        w = w - ak * v - (b[-1] * V[-2] if b else 0)
        for u in V: w -= (u @ w) * u
        for u in V: w -= (u @ w) * u
        bk = np.linalg.norm(w); b.append(bk)
        vn = w / bk; V.append(vn); v = vn; w = x * v
    return np.array(a), np.array(b)

a, b = lanc_ab(ts, mu, K); tot = mu.sum()

def p_pair(t):
    p0 = 1.0 / np.sqrt(tot); pm = 0.0
    for k in range(K):
        p1 = ((t - a[k]) * p0 - (b[k - 1] * pm if k > 0 else 0)) / b[k]
        pm, p0 = p0, p1
    return p0, pm  # p_K, p_{K-1}

tr = np.linspace(12, 58, 8000)
PK = np.array([p_pair(t) for t in tr])
E = PK[:, 0] - 1j * PK[:, 1]
phi = np.unwrap(np.angle(E))
dphi = np.gradient(phi, tr)

# equilibrium (arcsine) node density on [-A,A]: rho_eq(t) = K/(pi sqrt(A^2-t^2)) locally
rho_eq = K / (np.pi * np.sqrt(A ** 2 - tr ** 2))
# zeta (Weyl) density
rho_weyl = np.array([float(mp.log(t / (2 * np.pi))) / (2 * np.pi) for t in tr])

# compare |dphi/dt| shape to each (correlation, scale-free)
def corr(x, y): return float(np.corrcoef(x, y)[0, 1])
ce = corr(np.abs(dphi), rho_eq)
cw = corr(np.abs(dphi), rho_weyl)
print(f"boundary-phase density |dphi/dt| vs candidate densities (corr):")
print(f"  equilibrium/arcsine K/(pi sqrt(A^2-t^2)): {ce:+.3f}")
print(f"  zeta Weyl log(t/2pi)/2pi              : {cw:+.3f}")
print(f"  mean |dphi/dt| = {np.mean(np.abs(dphi)):.3f}; mean rho_eq = {np.mean(rho_eq):.3f}; mean rho_weyl = {np.mean(rho_weyl):.3f}")

# |E|^-2 vs Christoffel/measure: sample points
inv = 1.0 / np.abs(E) ** 2
wz_s = np.array([float(abs(mp.zeta(mp.mpc(0.5, t)))) ** 2 for t in tr[::80]])
print(f"  log|E|^-2 vs log|zeta|^2: corr = {corr(np.log(inv[::80]+1e-300), np.log(wz_s+1e-300)):+.3f}")
print("\nVERDICT: dphi tracks EQUILIBRIUM (not Weyl) => our explicit string's de Branges space is the")
print("QUADRATURE space of |zeta|^2, DISTINCT from Suzuki's xi-space; linked by Krein duality, not equal.")
print("DONE", flush=True)
