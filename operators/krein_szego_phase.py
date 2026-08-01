import numpy as np, os
from krein_fit import lanczos_jacobi, design
from scipy.fft import dct

# Szego-phase test: theory (generalized-Jacobi / Riemann-Hilbert) says the per-zero
# oscillation phase is set by the Szego function's phase = harmonic conjugate of
# (1/2)log w. On the zeta line the conjugate of log|zeta|^2 is ~ 2*arg zeta = S(T)
# content. Model: phi_k ~ s * psi(theta_k) + c0(window), with ONE scale s shared
# across ALL windows (offsets per window absorb theory constants).
# 51 phases vs 1 shared scale + 4 offsets.

DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
RUNS = [(30.0, 8000, 70, "A30_M8000"), (40.0, 10000, 95, "A40_M10000"),
        (60.0, 14000, 140, "A60_M14000"), (100.0, 20000, 220, "A100_M20000")]

data = []
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
    X = design(2 * thk, ns, False)
    beta, *_ = np.linalg.lstsq(X, db[5:], rcond=None)
    phi = np.arctan2(beta[2::2], beta[1::2])
    # conjugate function of (1/2) log w on the theta-circle via DCT/sine resum
    L = 0.5 * np.log(np.maximum(wz, 1e-300))
    a = dct(L, type=2, norm=None) / M          # a[k] ~ cosine coefficients (a[0]=2*mean)
    KF = min(3000, M - 1)
    ks = np.arange(1, KF)
    psi = np.array([np.sum(a[1:KF] * np.sin(ks * t)) for t in thk])
    data.append((A, thk, phi, psi))
    print(f"A={A:5.0f}: psi(theta_k) = {np.round(psi[:6], 3)} | fitted phi = {np.round(phi[:6], 3)}", flush=True)

def combR(s):
    tot, n = 0.0, 0
    for A, thk, phi, psi in data:
        r = np.abs(np.sum(np.exp(1j * (phi - s * psi))))   # offset absorbed by |.|
        tot += r; n += len(phi)
    return tot / n

ss = np.linspace(-3, 3, 1201)
Rs = np.array([combR(s) for s in ss])
i = int(np.argmax(Rs))
print(f"\nbest shared scale s = {ss[i]:+.3f} with combined resultant R = {Rs[i]:.3f}")
print(f"(control: s=0 gives R = {combR(0.0):.3f}; random phases would give ~{1/np.sqrt(np.mean([len(d[2]) for d in data])):.2f} per window)")
for s in [-2.0, -1.0, 1.0, 2.0]:
    print(f"  s={s:+.1f}: R={combR(s):.3f}")
top = np.argsort(Rs)[::-1][:5]
print("top scales:", [(round(float(ss[j]), 3), round(float(Rs[j]), 3)) for j in top])
