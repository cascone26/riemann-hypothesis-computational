"""Delta complete zero rebuild v2.
RS main sum (corrected: lam(n)/sqrt(n)) with two tapers; S-band completeness gate;
adaptive-dps exact Lambda spot checks; phases via sgn Z'; readout vs b(p)=-lam(p)."""
import numpy as np, mpmath as mp
from scipy.special import loggamma as slg
TWO_PI = 2*np.pi
TMAX = 260.0

N = 120
poly = [1] + [0]*N
for n in range(1, N+1):
    for _ in range(24):
        new = poly[:]
        for i in range(N, n-1, -1):
            new[i] -= poly[i-n]
        poly = new
tau = np.array([poly[n-1] for n in range(1, N+1)], dtype=float)
lam = tau / np.arange(1, N+1)**5.5
assert poly[1] == -24 and poly[2] == 252 and poly[4] == 4830

def theta_v(t):
    return np.imag(slg(6.0 + 1j*np.atleast_1d(t))) - np.atleast_1d(t)*np.log(TWO_PI)

def Z_rs(t, taper):
    t = np.atleast_1d(t)
    th = theta_v(t)
    out = np.zeros_like(t)
    ns = np.arange(1, N+1)
    for i, (ti, thi) in enumerate(zip(t, th)):
        Mi = ti/TWO_PI
        cut = ns <= Mi + taper
        n = ns[cut]
        w = np.ones(len(n))
        edge = n > Mi - taper
        w[edge] = 0.5*(1 + np.cos(np.pi*np.clip((n[edge]-(Mi-taper))/(2*taper), 0, 1)))
        out[i] = 2*np.sum(lam[cut]*w*np.cos(thi - ti*np.log(n))/np.sqrt(n))
    return out

grid = np.arange(8.0, TMAX, 0.01)
z3 = Z_rs(grid, 3.0)
z8 = Z_rs(grid, 8.0)
print("taper agreement: rms diff %.4f max %.4f" % (np.sqrt(np.mean((z3-z8)**2)), np.abs(z3-z8).max()))
idx = np.where(np.sign(z3[:-1]) != np.sign(z3[1:]))[0]
print("sign changes (taper3):", len(idx))

def bisect(f, a, b, tol=1e-7):
    fa = f(a)
    for _ in range(50):
        m = 0.5*(a+b); fm = f(m)
        if np.sign(fm) == np.sign(fa): a, fa = m, fm
        else: b = m
        if b-a < tol: break
    return 0.5*(a+b)

f3 = lambda x: Z_rs(np.array([x]), 3.0)[0]
f8 = lambda x: Z_rs(np.array([x]), 8.0)[0]
zeros, flags = [], []
for i in idx:
    r3 = bisect(f3, grid[i], grid[i+1])
    # taper-8 must also change sign nearby
    a8, b8 = r3-0.05, r3+0.05
    if np.sign(f8(a8)) == np.sign(f8(b8)):
        flags.append(r3); continue
    r8 = bisect(f8, a8, b8)
    zeros.append(0.5*(r3+r8))
    if abs(r3-r8) > 0.01: flags.append(r3)
zeros = np.array(sorted(zeros))
print("zeros: %d | taper-disagree flags: %d %s" % (len(zeros), len(flags), np.round(flags,2)[:8]))

k = np.arange(1, len(zeros)+1) - 0.5
S = k - theta_v(zeros)/np.pi
print("S band: [%.3f, %.3f] drift %.3f" % (S.min(), S.max(), abs(S[-10:].mean()-S[:10].mean())))

# exact Lambda spot check, adaptive dps, on 10 sampled zeros: sign change across [g-.004, g+.004]
def Lam_exact(t):
    mp.mp.dps = 45 + int(0.70*t)
    s = mp.mpf(6) + 1j*mp.mpf(t)
    tot = mp.mpf(0)
    for n in range(1, 61):
        x = TWO_PI*n
        tot += int(poly[n-1])*mp.re(mp.gammainc(s, x, mp.inf)*mp.power(x, -s)
                                    + mp.gammainc(12-s, x, mp.inf)*mp.power(x, s-12))
    return mp.re(tot)
samp = zeros[np.linspace(0, len(zeros)-1, 10).astype(int)]
ok = 0
for gz in samp:
    va, vb = Lam_exact(gz-0.004), Lam_exact(gz+0.004)
    good = mp.sign(va) != mp.sign(vb)
    ok += bool(good)
    print("  spot t=%.4f exact-sign-change=%s" % (gz, bool(good)))
print("exact spot checks: %d/10 PASS" % ok)

# phases via sgn Z' (from taper-8 sum, FD)
h = 1e-4
Zp = np.array([(f8(g+h)-f8(g-h))/(2*h) for g in zeros])
ph = -np.pi/2 - theta_v(zeros) + np.pi*(Zp < 0)
ph = np.mod(ph + np.pi, 2*np.pi) - np.pi
# parity check (non-circular completeness cross-check): sgn Z' should strictly alternate
alt = np.all(np.sign(Zp[1:]) == -np.sign(Zp[:-1]))
print("sgn Z' strict alternation:", alt)

np.save("Delta_zeros_complete.npy", zeros)
np.save("Delta_argLp_complete.npy", ph)

# READOUT vs b(n): b(p)=-lam(p), b(p^2)=+1, b(else composite)=0-ish per inverse Euler
kk = np.arange(1, len(zeros)+1)
w = (-1)**kk * np.exp(1j*theta_v(zeros))
c = None
print("elementary readout on rebuilt COMPLETE Delta list:")
lamn = lam  # analytic lambda(n)
tests = [(2, -lamn[1]), (3, -lamn[2]), (4, 1.0), (5, -lamn[4]), (7, -lamn[6]), (9, 1.0), (6, lamn[1]*lamn[2])]
for n, btrue in tests:
    G = (w*np.exp(1j*np.log(n)*zeros)).mean()
    if c is None: c = np.angle(G) - (0 if btrue > 0 else np.pi)  # anchor global const at n=2
    v = (G*np.exp(-1j*c)).real
    match = "OK" if np.sign(v) == np.sign(btrue) else "MISS"
    print(f"  n={n:<2} b_true={btrue:+.3f}  readout={v:+.4f}  {match}")
print("DONE")
