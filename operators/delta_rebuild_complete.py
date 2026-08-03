"""Rebuild the COMPLETE Delta (Ramanujan tau) zero list up to T=260 + arg L' at each zero.

Old Delta_zeros.npy (212 zeros) misses ~16 zeros starting gamma~69 (S-drift audit).
Method:
  scan: RS-style main sum Z(t) ~ 2 sum_{n<=M(t)} lam(n)/sqrt(n) cos(theta(t) - t log n),
        tapered cutoff at M = t/2pi + taper; grid 0.02; bisect sign changes.
  verify: EXACT incomplete-gamma Lambda(6+it) (classical s, weight 12, eps=+1):
        Lambda(s) = sum_n tau(n) [ (2pi n)^{-s} Gamma(s,2pi n) + (2pi n)^{s-12} Gamma(12-s,2pi n) ]
        real on the critical line; re-bisect each zero with exact Lambda to 1e-8.
  gate: S(k) = (k-1/2) - theta(g_k)/pi must stay in a bounded O(1) band over the WHOLE range.
  phases: L_cl'(6+ig) via central FD of Lambda/gammafactor along the line (dL/ds = -i dL/dt).
Outputs: Delta_zeros_complete.npy, Delta_argLp_complete.npy
"""
import numpy as np, mpmath as mp
mp.mp.dps = 25
TWO_PI = 2*np.pi
TMAX = 260.0

# ---- tau(n) via eta^24: q prod (1-q^n)^24 ----
N = 120
poly = [1] + [0]*N
for n in range(1, N+1):
    for _ in range(24):
        new = poly[:]
        for i in range(N, n-1, -1):
            new[i] -= poly[i-n]
        poly = new
tau = np.array([poly[n-1] for n in range(1, N+1)], dtype=float)  # tau(n) = coeff of q^n in q*prod
lam = tau / np.arange(1, N+1)**5.5  # analytic normalization
print("tau check: tau(2)=%d tau(3)=%d tau(5)=%d" % (poly[1], poly[2], poly[4]))
assert poly[1] == -24 and poly[2] == 252 and poly[4] == 4830, "tau FAILED"

# ---- theta and RS main sum (vectorized scan) ----
from scipy.special import loggamma as slg
def theta_v(t):
    return np.imag(slg(6.0 + 1j*t)) - t*np.log(TWO_PI)
def Z_rs(t):
    t = np.atleast_1d(t)
    th = theta_v(t)
    M = t/TWO_PI
    out = np.zeros_like(t)
    ns = np.arange(1, N+1)
    for i, (ti, thi, Mi) in enumerate(zip(t, th, M)):
        cut = ns <= Mi + 3.0
        n = ns[cut]
        w = np.ones(len(n))
        edge = n > Mi - 3.0
        w[edge] = 0.5*(1 + np.cos(np.pi*np.clip((n[edge]-(Mi-3.0))/6.0, 0, 1)))
        out[i] = 2*np.sum(lam[cut]*w*np.cos(thi - ti*np.log(n))/1.0)
    return out

# ---- exact Lambda ----
def Lam(t):
    s = mp.mpf(6) + 1j*mp.mpf(t)
    tot = mp.mpf(0)
    for n in range(1, 61):
        x = TWO_PI*n
        a = mp.gammainc(s, x, mp.inf) * mp.power(x, -s)
        b = mp.gammainc(12-s, x, mp.inf) * mp.power(x, s-12)
        tot += int(poly[n-1]) * (a + b)
    assert abs(mp.im(tot)) < 1e-10*max(1, abs(mp.re(tot))), "Lambda not real!"
    return float(mp.re(tot))

# ---- scan ----
grid = np.arange(6.0, TMAX, 0.02)
zv = Z_rs(grid)
idx = np.where(np.sign(zv[:-1]) != np.sign(zv[1:]))[0]
print("RS scan sign changes:", len(idx))

zeros = []
for i in idx:
    a, b = grid[i], grid[i+1]
    fa = Lam(a); fb = Lam(b)
    if np.sign(fa) == np.sign(fb):
        # RS/exact disagree near here: fine-scan exact +-0.06
        tt = np.arange(a-0.06, b+0.06, 0.01)
        vv = [Lam(x) for x in tt]
        jj = [j for j in range(len(vv)-1) if np.sign(vv[j]) != np.sign(vv[j+1])]
        if not jj:
            print("  DROPPED spurious RS crossing near %.3f" % a); continue
        a, b = tt[jj[0]], tt[jj[0]+1]; fa, fb = vv[jj[0]], vv[jj[0]+1]
    for _ in range(40):
        m = 0.5*(a+b); fm = Lam(m)
        if np.sign(fm) == np.sign(fa): a, fa = m, fm
        else: b, fb = m, fm
        if b-a < 1e-9: break
    zeros.append(0.5*(a+b))
zeros = np.array(sorted(zeros))
print("exact-verified zeros:", len(zeros), "first:", zeros[:3], "last:", zeros[-1])

# ---- completeness gate: S bounded ----
k = np.arange(1, len(zeros)+1) - 0.5
S = k - theta_v(zeros)/np.pi
print("S band: min %.3f max %.3f  (PASS if width < ~2 and no drift)" % (S.min(), S.max()))
drift = abs(S[-10:].mean() - S[:10].mean())
print("S drift end-start: %.3f  (PASS if < 0.7)" % drift)

# ---- old-list containment ----
old = np.load("Delta_zeros.npy")
d = np.array([np.min(np.abs(zeros - og)) for og in old])
print("old zeros matched within: max %.4f (n=%d)" % (d.max(), len(old)))
missing_new = len(zeros) - len(old)
print("net new zeros: %d" % missing_new)

# ---- phases: arg L'(1/2+i g) analytic = arg L_cl'(6+ig) ----
def LamC(t):
    s = mp.mpf(6) + 1j*mp.mpf(t)
    tot = mp.mpc(0)
    for n in range(1, 61):
        x = TWO_PI*n
        tot += int(poly[n-1]) * (mp.gammainc(s, x, mp.inf)*mp.power(x, -s)
                                 + mp.gammainc(12-s, x, mp.inf)*mp.power(x, s-12))
    return tot
def argLp(t):
    h = 1e-5
    gam = lambda tt: mp.power(TWO_PI, -(6+1j*mp.mpf(tt))) * mp.gamma(6+1j*mp.mpf(tt))
    Lp = (LamC(t+h)/gam(t+h) - LamC(t-h)/gam(t-h)) / (2*h)   # dL/dt
    Lps = -1j*Lp                                              # dL/ds = -i dL/dt... s=sigma+it => dL/ds along line: dL/dt * (1/i)
    return float(mp.arg(Lps))
ph = np.array([argLp(g) for g in zeros])

np.save("Delta_zeros_complete.npy", zeros)
np.save("Delta_argLp_complete.npy", ph)

# ---- immediate S-only retest ----
d541 = ph - np.pi*(k - theta_v(zeros)/np.pi)
z = np.exp(1j*d541).mean()
R = abs(z)
sd = np.std(np.mod(d541 - np.angle(z) + np.pi, 2*np.pi) - np.pi)
print("POSITION-ONLY retest on rebuilt list: R=%.4f  sd=%.4f  const=%.4f" % (R, sd, np.angle(z)))
print("DONE")
