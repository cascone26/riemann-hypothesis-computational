"""Rigidity probe: displace one zero by 0.5 spacing; run v3 fine-delta solve from the
perturbed configuration. Does the prime field pull it back? Measure per-zero recovery."""
import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
g = np.load("zeros10k.npy")[:5000]
N = len(g); k = np.arange(1, N+1)
th_of = lambda tt: np.imag(slg(0.25+0.5j*np.asarray(tt))) - np.asarray(tt)/2*np.log(np.pi)
sp = np.gradient(g)
osc = []
for p in primerange(2, 500000):
    lp = np.log(p); m = 1
    while p**m < 500000:
        osc.append((lp*m, 1.0/(np.pi*np.sqrt(p**m)*m))); m += 1
Ls = np.array([o[0] for o in osc]); As = np.array([o[1] for o in osc])
def S_Sp(t, delta):
    w = As*np.exp(-0.5*(Ls*delta)**2)
    keep = w > 1e-8
    S = np.zeros_like(t); Sp = np.zeros_like(t)
    for L, a in zip(Ls[keep], w[keep]):
        S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
    return S, Sp
IDX = 2500
t = g.copy()
t[IDX] += 0.5*sp[IDX]        # displace one zero by half a spacing
print("displaced zero %d by +0.5 spacing" % IDX)
for delta in (0.25, 0.18, 0.12):
    for it in range(3):
        S, Sp = S_Sp(t, delta)
        thp = 0.5*np.log(t/(2*np.pi))/np.pi
        J = np.maximum(thp + Sp, 0.3*thp)
        t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/J, -1.5, 1.5)
err = (t - g)/sp
print("after solve: displaced zero error = %+.4f spacings (was +0.5)" % err[IDX])
print("neighbors err: %s" % np.round(err[IDX-3:IDX+4], 4))
print("global rms (excl. displaced +-3): %.4f" % np.sqrt(np.mean(np.delete(err, range(IDX-3,IDX+4))**2)))
print("DONE")
