import numpy as np, os
import mpmath as mp

# DIRECTION 3: de Branges E-function of the |zeta|^2 string vs Suzuki's Theta.
# From the truncated Jacobi matrix J_K (string coefficients a,b), build the
# Hermite-Biehler pair via first- and second-kind orthonormal polynomials:
#   p_n (first kind, real on R), q_n (second kind, real on R, q_0=0).
#   A_K = q_K,  B_K = p_K  =>  E_K(z) = A_K(z) - i B_K(z)  (Hermite-Biehler).
# Test 1: zeros of E_K in the lower half-plane (HB property).
# Test 2: boundary phase phi_K(t) = -arg E_K(t); does phi_K'(t) track pi*(Weyl density)
#         and do the |zeta|^2 zeros show as phase plateaus?
# Test 3: |E_K(t)|^{-2} ~ the Christoffel/measure density (de Branges).
# Compare structure to Suzuki Theta = xi(1/2-iz)/xi(1/2+iz) (inner, |Theta|=1 on R).

mp.mp.dps = 15
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

a, b = lanc_ab(ts, mu, K)
tot = mu.sum()

def pq(z):
    # first & second kind, orthonormal recurrence; return p_K, q_K
    p0 = 1.0 / np.sqrt(tot); pm1 = 0.0
    q0 = 0.0; qm1 = -1.0 / (np.sqrt(tot))  # seed so q tracks second kind
    p = p0; pm = pm1; q = q0; qm = qm1
    for k in range(K):
        p1 = ((z - a[k]) * p - (b[k - 1] * pm if k > 0 else 0)) / b[k]
        q1 = ((z - a[k]) * q - (b[k - 1] * qm if k > 0 else 0)) / b[k]
        pm, p = p, p1; qm, q = q, q1
    return p, q

# Test 1: HB property — zeros of E_K = A - iB = q_K - i p_K in lower half-plane?
zim = []
for zr in np.linspace(14, 58, 200):
    # crude: does E have a near-zero just below the axis? sample |E| at zr - i*eps grid
    pass
# Instead: eigen of J_K are real (nodes). Build E on a fine real grid + complex strip.
tr = np.linspace(12, 58, 4000)
P = np.array([pq(t) for t in tr])
pK = P[:, 0]; qK = P[:, 1]
E = qK - 1j * pK
phi = -np.unwrap(np.angle(E))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
gam = [g for g in Z if 12 < g < 58]
# Weyl density on this window (theta_H'/pi)
dens = np.array([float(mp.log(t / (2 * np.pi))) / (2 * np.pi) for t in tr])
phi_pred = np.cumsum(dens) * (tr[1] - tr[0]) * np.pi * 2  # 2 for +- pairing
# align slopes
sl = np.polyfit(tr, phi, 1)[0]; slp = np.polyfit(tr, phi_pred, 1)[0]
print(f"boundary phase slope: measured {sl:.3f}, Weyl-pred (2pi*dens) {slp:.3f}, ratio {sl/slp:.3f}", flush=True)
# |E|^-2 vs weight
invE2 = 1.0 / np.abs(E) ** 2
wcurve = np.array([float(abs(mp.zeta(mp.mpc(0.5, t)))) ** 2 for t in tr[::40]])
from numpy import corrcoef
r = corrcoef(np.log(invE2[::40] + 1e-30), np.log(wcurve + 1e-30))[0, 1]
print(f"log|E_K|^-2 vs log|zeta|^2 on window: corr = {r:.3f} (de Branges: |E|^-2 ~ measure density)", flush=True)
# phase plateaus at zeros? measure local phase slope dip at each gamma
dphi = np.gradient(phi, tr)
for g in gam[:8]:
    i = int(np.argmin(np.abs(tr - g)))
    loc = dphi[i]; med = np.median(dphi)
    print(f"  gamma={g:.3f}: phase-slope/median = {loc/med:.3f}", flush=True)
print("DONE", flush=True)
