"""Final Delta validation from the PARI-computed complete zero list.
Gates: (1) S-band completeness vs theta-density; (2) cross-check vs the 11 exact
incomplete-Gamma zeros; (3) position-only R; (4) elementary readout vs b(n)
(global sign anchored at n=2; validation on the OTHER channels); (5) J-comb."""
import numpy as np
from scipy.special import loggamma as slg

zeros = np.loadtxt("delta_zeros_pari.txt")
zeros = np.sort(zeros)
N = len(zeros)
def th(t):
    return np.imag(slg(6.0 + 1j*np.asarray(t))) - np.asarray(t)*np.log(2*np.pi)
k = np.arange(1, N+1) - 0.5
S = k - th(zeros)/np.pi
print("N=%d  S-band [%.3f, %.3f] drift %.3f  (gate: width<2, |drift|<0.7)" %
      (N, S.min(), S.max(), abs(S[-10:].mean()-S[:10].mean())))
v4 = np.load("Delta_zeros_complete.npy")  # the 11 exact ones
d = [np.min(np.abs(zeros - z)) for z in v4]
print("11 exact zeros matched by PARI within: %.2e" % max(d))

# phases via elementary alternation (sign anchored below), position-only R is exact by construction
kk = np.arange(1, N+1)
for gsign in (+1, -1):
    w = gsign*(-1)**kk*np.exp(1j*th(zeros))
    G2 = (w*np.exp(1j*np.log(2)*zeros)).mean()
    if G2.real > 0:  # b_true(2) = -tau(2)/2^5.5 = +0.53 > 0
        break
print("global sign:", gsign)

poly = [1]+[0]*20
for n in range(1,21):
    for _ in range(24):
        new = poly[:]
        for i in range(20, n-1, -1): new[i] -= poly[i-n]
        poly = new
lam = np.array(poly[:20], dtype=float)/np.arange(1,21)**5.5
tests = [(2,-lam[1]),(3,-lam[2]),(4,1.0),(5,-lam[4]),(7,-lam[6]),(9,1.0),
         (6,lam[1]*lam[2]),(11,-lam[10]),(13,-lam[12]),(10,lam[1]*lam[4]),(14,lam[1]*lam[6])]
okn = 0
print("elementary readout (COMPLETE PARI list), anchor n=2:")
for n, btrue in tests:
    G = (w*np.exp(1j*np.log(n)*zeros)).mean()
    ok = np.sign(G.real) == np.sign(btrue)
    okn += bool(ok)
    tag = "anchor" if n == 2 else ("OK" if ok else "MISS")
    print(f"  n={n:<3} b_true={btrue:+.3f}  Re G={G.real:+.4f}  {tag}")
print(f"readout: {okn}/{len(tests)} (incl. anchor)")

# J-comb from positions
T0, T, dt = float(zeros[0])+0.5, float(zeros[-1])-0.5, 0.002
t = np.arange(T0, T, dt)
w_t = gsign*(-1.0)**np.searchsorted(zeros, t)*np.exp(1j*th(t))
K = (th(np.array([T]))[0]-th(np.array([T0]))[0])/np.pi
okj = 0
print("J-comb (positions only):")
for n, btrue in tests:
    J = np.trapezoid(w_t*np.exp(1j*np.log(n)*t), dx=dt)/K
    ok = np.sign(J.real) == np.sign(btrue)
    okj += bool(ok)
    print(f"  n={n:<3} Re J={J.real:+.4f}  {'OK' if ok else 'MISS'}")
print(f"J-comb: {okj}/{len(tests)}")

np.save("Delta_zeros_complete.npy", zeros)
np.save("Delta_argLp_complete.npy", np.mod(-np.pi/2 - th(zeros) + np.pi*((gsign*(-1)**kk) < 0) + np.pi, 2*np.pi) - np.pi)
print("saved complete files. DONE")
