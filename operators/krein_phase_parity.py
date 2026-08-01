import numpy as np, os
from krein_fit import lanczos_jacobi, design

# Parity/counting test: is phi_k = pi*k + smooth(theta_k)? (Szego phase winds by pi
# per interior zero => string phases carry the zero-COUNTING function.)
# Compare circular-fit resultant R for models:
#   M0: phi ~ c0                    M1: phi ~ c1*theta + c0
#   M2: phi - pi*k ~ c0             M3: phi - pi*k ~ c1*theta + c0

DIR = os.path.dirname(os.path.abspath(__file__))
Z = np.load(os.path.join(DIR, "zeros500.npy"))
RUNS = [(30.0, 8000, 70, "A30_M8000"), (40.0, 10000, 95, "A40_M10000"),
        (60.0, 14000, 140, "A60_M14000"), (100.0, 20000, 220, "A100_M20000")]

def bestR(ph, thk, linear):
    if not linear:
        return float(np.abs(np.mean(np.exp(1j * ph)))), 0.0
    c1s = np.linspace(-8, 8, 3201)
    Rs = np.array([np.abs(np.mean(np.exp(1j * (ph - c1 * thk)))) for c1 in c1s])
    i = int(np.argmax(Rs))
    return float(Rs[i]), float(c1s[i])

print(f"{'A':>5} {'nz':>3} | {'M0':>6} {'M1':>6} {'M2:pi*k':>8} {'M3:pi*k+lin':>11} | c1(M3)")
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
    k = np.arange(1, len(gam) + 1)
    R0, _ = bestR(phi, thk, False)
    R1, _ = bestR(phi, thk, True)
    R2, _ = bestR(phi - np.pi * k, thk, False)
    R3, c1 = bestR(phi - np.pi * k, thk, True)
    print(f"{A:>5.0f} {len(gam):>3} | {R0:>6.3f} {R1:>6.3f} {R2:>8.3f} {R3:>11.3f} | {c1:+.2f}")
    resid = np.angle(np.exp(1j * (phi - np.pi * k - c1 * thk)))
    resid -= np.mean(resid)
    print(f"      M3 residual phases: {np.round(resid, 2)}")
print("\nR2/R3 >> R0/R1 => phases = pi*(zero index) + smooth: the string counts zeros (N(T) imprint).")
