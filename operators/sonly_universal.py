"""S-only universality: does the POSITION-ONLY phase reconstruct the true arg L'(rho)
across the full Selberg-class family set (degrees 1-4)?

Route A (exact theta, gamma factor known): S_pos(k) = (k-1/2) - theta_L(g_k)/pi + C
Route B (self-detrend, NO gamma knowledge): fit (k-1/2) ~ a*t*log t + b*t + c*log t + d
  (the universal main-term shape), S_pos = residual. Only 4 global params -> window-scale
  S fluctuations (the carrier, per LEDGER 602) survive intact.

Decisive metric: R = |circ-mean of exp(i*(argLp_k - pi*S_pos(k)))|.
R ~ 1  <=>  position phase = true phase up to ONE constant  <=>  every readout
in LEDGER 570-583 carries over verbatim with no L-derivative data.
"""
import numpy as np, mpmath as mp
mp.mp.dps = 30

def circR(d):
    z = np.exp(1j*d).mean()
    return abs(z), np.angle(z)

def routeB(g, ph):
    k = np.arange(1, len(g)+1) - 0.5
    X = np.column_stack([g*np.log(g), g, np.log(g), np.ones_like(g)])
    coef, *_ = np.linalg.lstsq(X, k, rcond=None)
    S = k - X@coef
    R, mu = circR(ph - np.pi*S)
    return R, np.std(np.mod(ph - np.pi*S - mu + np.pi, 2*np.pi) - np.pi)

def routeA(g, ph, theta):
    k = np.arange(1, len(g)+1) - 0.5
    th = np.array([float(theta(t)) for t in g])
    S = k - th/np.pi
    R, mu = circR(ph - np.pi*S)
    return R, np.std(np.mod(ph - np.pi*S - mu + np.pi, 2*np.pi) - np.pi)

def th_dirichlet(q, a):
    return lambda t: mp.im(mp.loggamma((0.5+1j*t+a)/2)) - (t/2)*mp.log(mp.pi/q)

def th_delta(t):  # Gamma_C(s+11/2), analytic norm
    return mp.im(mp.loggamma(0.5+1j*t+5.5)) - t*mp.log(2*mp.pi)

def th_sym2(t):  # Gamma_R(s+1) Gamma_C(s+11)
    return (mp.im(mp.loggamma((0.5+1j*t+1)/2)) - (t/2)*mp.log(mp.pi)
            + mp.im(mp.loggamma(0.5+1j*t+11)) - t*mp.log(2*mp.pi))

FAM = [
    ("zeta      d=1", "zeros10k.npy", "argzp10k.npy", lambda t: mp.siegeltheta(t)),
    ("chi3      d=1", "Lchi3_zeros_complete.npy", "Lchi3_argLp_complete.npy", th_dirichlet(3,1)),
    ("chi4      d=1", "Lchi4_zeros_complete.npy", "Lchi4_argLp_complete.npy", th_dirichlet(4,1)),
    ("Delta     d=2", "Delta_zeros.npy", "Delta_argLp.npy", th_delta),
    ("sym2Delta d=3", "sym2_zeros.npy", "sym2_argLp.npy", th_sym2),
    ("sym3Delta d=4", "sym3_zeros.npy", "sym3_argLp.npy", None),
]

print(f"{'family':<14} {'N':>5}  {'R_exact':>8} {'sd_ex':>6}   {'R_detrend':>9} {'sd_dt':>6}")
for name, zf, pf, theta in FAM:
    g = np.load(zf); ph = np.load(pf)
    n = min(len(g), len(ph)); g, ph = g[:n], ph[:n]
    if name.startswith("zeta"):
        gz = g[:2000]; phz = ph[:2000]  # cap exact route cost; detrend uses all
        Ra, sda = routeA(gz, phz, theta)
    elif theta is not None:
        Ra, sda = routeA(g, ph, theta)
    else:
        Ra, sda = float('nan'), float('nan')
    Rb, sdb = routeB(g, ph)
    print(f"{name:<14} {n:>5}  {Ra:8.4f} {sda:6.3f}   {Rb:9.4f} {sdb:6.3f}")
print("DONE")
