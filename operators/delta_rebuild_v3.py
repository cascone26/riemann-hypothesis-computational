"""Delta rebuild v3: RS sum = bracket finder ONLY; every zero refined with the EXACT
incomplete-gamma Lambda at adaptive dps (45 + 0.70*t digits, beats the e^{-pi t/2} scale),
multiprocessing over zeros. Gates: S-band bounded, first zeros vs LMFDB, readout vs b(n)."""
import numpy as np
from scipy.special import loggamma as slg
from multiprocessing import Pool
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
TAU = [int(poly[n-1]) for n in range(1, 61)]
tau = np.array(poly[:N], dtype=float)
lam = tau / np.arange(1, N+1)**5.5

def theta_v(t):
    return np.imag(slg(6.0 + 1j*np.asarray(t, dtype=float))) - np.asarray(t, dtype=float)*np.log(TWO_PI)

def Z_rs(t, taper=8.0):
    t = np.atleast_1d(t).astype(float)
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

def refine(bracket):
    a, b = bracket
    import mpmath as mp
    def Lam(t):
        mp.mp.dps = 45 + int(0.70*t)
        s = mp.mpf(6) + 1j*mp.mpf(t)
        tot = mp.mpf(0)
        for n, tn in enumerate(TAU, start=1):
            x = TWO_PI*n
            tot += tn*mp.re(mp.gammainc(s, x, mp.inf)*mp.power(x, -s)
                            + mp.gammainc(12-s, x, mp.inf)*mp.power(x, s-12))
        return tot
    fa, fb = Lam(a), Lam(b)
    if (fa > 0) == (fb > 0):
        # widen once
        a2, b2 = a-0.10, b+0.10
        fa, fb = Lam(a2), Lam(b2)
        if (fa > 0) == (fb > 0):
            return None
        a, b = a2, b2
    for _ in range(24):
        m = 0.5*(a+b); fm = Lam(m)
        if (fm > 0) == (fa > 0): a, fa = m, fm
        else: b, fb = m, fm
        if b - a < 2e-6: break
    return 0.5*(a+b)

if __name__ == "__main__":
    grid = np.arange(8.0, TMAX, 0.01)
    z = Z_rs(grid)
    idx = np.where(np.sign(z[:-1]) != np.sign(z[1:]))[0]
    brackets = [(grid[i]-0.04, grid[i+1]+0.04) for i in idx]
    print("RS brackets:", len(brackets), flush=True)
    with Pool(8) as p:
        res = p.map(refine, brackets)
    dropped = sum(r is None for r in res)
    zeros = np.array(sorted(r for r in res if r is not None))
    # dedup (widened brackets can converge to same zero)
    zeros = zeros[np.concatenate([[True], np.diff(zeros) > 1e-3])]
    print("exact zeros: %d (dropped %d)" % (len(zeros), dropped), flush=True)
    print("first five:", np.round(zeros[:5], 5), "| LMFDB truth: 9.22238 13.90755 17.44277 19.65651 22.33610")

    k = np.arange(1, len(zeros)+1) - 0.5
    S = k - theta_v(zeros)/np.pi
    print("S band: [%.3f, %.3f] drift %.3f  (PASS if width<2, drift<0.7)" %
          (S.min(), S.max(), abs(S[-10:].mean() - S[:10].mean())))

    # phases via exact Z' sign: Z(t) = Lam(6+it)/|gamma-part| real; sgn Z' = sgn dLam/dt
    import mpmath as mp
    def Lam_s(t):
        mp.mp.dps = 45 + int(0.70*t)
        s = mp.mpf(6) + 1j*mp.mpf(t)
        return sum(tn*mp.re(mp.gammainc(s, x, mp.inf)*mp.power(x, -s)
                            + mp.gammainc(12-s, x, mp.inf)*mp.power(x, s-12))
                   for n, tn in enumerate(TAU, start=1) for x in [TWO_PI*n])
    def sgnzp(gz):
        return 1.0 if Lam_s(gz+0.003) > Lam_s(gz-0.003) else -1.0
    with Pool(8) as p:
        sg = np.array(p.map(sgnzp, list(zeros)))
    alt = np.all(sg[1:] == -sg[:-1])
    print("sgn Lambda' strict alternation:", bool(alt))
    ph = -np.pi/2 - theta_v(zeros) + np.pi*(sg < 0)
    ph = np.mod(ph + np.pi, 2*np.pi) - np.pi

    np.save("Delta_zeros_complete.npy", zeros)
    np.save("Delta_argLp_complete.npy", ph)

    kk = np.arange(1, len(zeros)+1)
    w = (-1)**kk * np.exp(1j*theta_v(zeros))
    c = None
    print("elementary readout, COMPLETE exact Delta list:")
    tests = [(2, -lam[1]), (3, -lam[2]), (4, 1.0), (5, -lam[4]), (7, -lam[6]),
             (9, 1.0), (6, lam[1]*lam[2]), (11, -lam[10]), (13, -lam[12])]
    okn = 0
    for n, btrue in tests:
        G = (w*np.exp(1j*np.log(n)*zeros)).mean()
        if c is None:
            c = np.angle(G) - (0 if btrue > 0 else np.pi)
        v = (G*np.exp(-1j*c)).real
        ok = np.sign(v) == np.sign(btrue)
        okn += bool(ok)
        print(f"  n={n:<2} b_true={btrue:+.3f}  readout={v:+.4f}  {'OK' if ok else 'MISS'}")
    print(f"readout: {okn}/{len(tests)}")
    print("DONE")
