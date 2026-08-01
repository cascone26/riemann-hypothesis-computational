import numpy as np, os
import mpmath as mp
from scipy.fft import dct

# Identify the intrinsic phase: the DCT-conjugate psi(theta_k) is ~A-independent
# per zero => it is a function of gamma alone. Candidates (computed by mpmath):
#   P1: arg zeta'(1/2 + i gamma_k)
#   P2: pi*S(gamma_k) fractional content, S = N(T)-1-theta_H(T)/pi at T=gamma_k
#   P3: -theta_H(gamma_k) mod pi
# Test: circular corr of psi_avg(gamma_k) (averaged over the 4 windows) vs each,
# scanning one shared scale in {+-1, +-2} and offset.

mp.mp.dps = 20
DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
RUNS = [(30.0, 8000, "A30_M8000"), (40.0, 10000, "A40_M10000"),
        (60.0, 14000, "A60_M14000"), (100.0, 20000, "A100_M20000")]

# psi from each window's weight (DCT conjugate), averaged per zero
psis = {}
for A, M, tag in RUNS:
    th = (np.arange(M) + 0.5) * np.pi / M
    wz = np.load(os.path.join(DIR, f"zw_{tag}.npy"))
    L = 0.5 * np.log(np.maximum(wz, 1e-300))
    a = dct(L, type=2, norm=None) / M
    KF = min(3000, M - 1); ks = np.arange(1, KF)
    for g in [float(x) for x in Z if x < A - 0.5]:
        t = np.arccos(g / A)
        psis.setdefault(round(g, 3), []).append(float(np.sum(a[1:KF] * np.sin(ks * t))))
gams = sorted(psis)
psi_avg = np.array([np.mean(psis[g]) for g in gams])
spread = np.array([np.std(psis[g]) for g in gams])
print(f"{len(gams)} zeros; psi window-spread median = {np.median(spread):.3f} (small => intrinsic)")

argzp, Sg, thH = [], [], []
for g in gams:
    zp = mp.zeta(mp.mpc(0.5, g), derivative=1)
    argzp.append(float(mp.arg(zp)))
    thH.append(float(mp.siegeltheta(g)))
argzp = np.array(argzp); thH = np.array(thH)
k = np.arange(1, len(gams) + 1)
Sg = k - 1 - thH / np.pi          # S(T) just below the k-th zero (up to integer conv)

def circfit(x, target):
    best = (0, 0, 0.0)
    for s in [-2, -1, 1, 2, -0.5, 0.5]:
        R = np.abs(np.mean(np.exp(1j * (target - s * x))))
        if R > best[2]:
            best = (s, float(np.angle(np.mean(np.exp(1j * (target - s * x))))), float(R))
    return best

for name, x in [("arg zeta'(rho)", argzp), ("pi*S(gamma)", np.pi * Sg),
                ("theta_H mod pi", np.mod(thH, np.pi)), ("arg zeta' + theta_H", argzp + thH)]:
    s, c0, R = circfit(x, psi_avg)
    print(f"  psi vs {name:>18}: best s={s:+.1f}  R={R:.3f}")
print(f"  (chance for {len(gams)} pts ~ {1/np.sqrt(len(gams)):.2f}; psi self-consistency R across windows would be the ceiling)")
print("\npsi_avg:", np.round(psi_avg, 3))
print("argzp:  ", np.round(argzp, 3))
