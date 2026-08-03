"""Post-rebuild validation: (1) position-only R test on the exact Delta list -> 6/6 table;
(2) Delta J-comb from POSITIONS ONLY: square wave (-1)^{N(t)} e^{i theta(t)} built from the
zero list alone, plain Fourier transform vs b(n) signs."""
import numpy as np
from scipy.special import loggamma as slg

def th_delta(t):
    return np.imag(slg(6.0 + 1j*np.asarray(t))) - np.asarray(t)*np.log(2*np.pi)

g = np.load("Delta_zeros_complete.npy")
ph = np.load("Delta_argLp_complete.npy")
N = len(g)
k = np.arange(1, N+1) - 0.5
S = k - th_delta(g)/np.pi
d = ph - np.pi*S
z = np.exp(1j*d).mean()
sd = np.std(np.mod(d - np.angle(z) + np.pi, 2*np.pi) - np.pi)
print("Delta exact list: N=%d  position-only R=%.4f sd=%.4f" % (N, abs(z), sd))

# J-comb from positions only
T0, T, dt = float(g[0]) + 0.5, float(g[-1]) - 0.5, 0.002
t = np.arange(T0, T, dt)
Ncount = np.searchsorted(g, t)          # N(t)
th = th_delta(t)
w = (-1.0)**Ncount * np.exp(1j*th)
K = (th[-1] - th[0])/np.pi
Ntau = 20
poly = [1] + [0]*Ntau
for n in range(1, Ntau+1):
    for _ in range(24):
        new = poly[:]
        for i in range(Ntau, n-1, -1):
            new[i] -= poly[i-n]
        poly = new
lam = np.array(poly[:Ntau], dtype=float)/np.arange(1, Ntau+1)**5.5
tests = [(2, -lam[1]), (3, -lam[2]), (4, 1.0), (5, -lam[4]), (7, -lam[6]),
         (9, 1.0), (6, lam[1]*lam[2]), (11, -lam[10]), (13, -lam[12])]
c = None
okn = 0
print("Delta J-comb (positions only):")
for n, btrue in tests:
    J = np.trapezoid(w*np.exp(1j*np.log(n)*t), dx=dt)/K
    if c is None:
        c = np.angle(J) - (0 if btrue > 0 else np.pi)
    v = (J*np.exp(-1j*c)).real
    ok = np.sign(v) == np.sign(btrue)
    okn += bool(ok)
    print(f"  n={n:<2} b_true={btrue:+.3f}  Re J={v:+.4f}  {'OK' if ok else 'MISS'}")
Joff = np.trapezoid(w*np.exp(1j*1.25*t), dx=dt)/K
print(f"  off-support u=1.25: |J|={abs(Joff):.4f}")
print(f"J-comb: {okn}/{len(tests)}")
print("DONE")
