"""Off-line (Lehmer-collision) probe: zeros k=2500,2501 merged to their midpoint (the
double-ordinate signature of an off-line pair). Run v3 at fine delta ONLY from the merged
configuration: repaired (no spurious fixed point) or stuck (spurious solution exists)?"""
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
I, J2 = 2500, 2501
t = g.copy()
mid = 0.5*(g[I]+g[J2])
t[I] = t[J2] = mid
print("merged zeros %d,%d at midpoint (displacements %+.3f/%+.3f spacings)" %
      (I, J2, (mid-g[I])/sp[I], (mid-g[J2])/sp[J2]))
def resid_norm(t, delta):
    S, _ = S_Sp(t, delta)
    return np.sqrt(np.mean(((th_of(t)/np.pi + 1 + S) - (k-0.5))**2))
print("initial residual norm (delta=0.12): %.5f" % resid_norm(t, 0.12))
for delta in (0.12,):
    for it in range(8):
        S, Sp = S_Sp(t, delta)
        thp = 0.5*np.log(t/(2*np.pi))/np.pi
        Jc = np.maximum(thp + Sp, 0.3*thp)
        t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/Jc, -1.5, 1.5)
        e = (t - g)/sp
        print("iter %d: merged-pair errors %+.4f %+.4f | global rms %.4f | resid %.5f" %
              (it, e[I], e[J2], np.sqrt(np.mean(e**2)), resid_norm(t, delta)), flush=True)
print("DONE")
