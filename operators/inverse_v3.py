"""Inverse map v3: damped Newton with full Jacobian (theta'/pi + S'_delta), trust region.
Target: beat 0.059 spacings; test the plateau law."""
import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
g = np.load("zeros10k.npy")[:5000]
N = len(g); k = np.arange(1, N+1)
th_of = lambda tt: np.imag(slg(0.25+0.5j*np.asarray(tt))) - np.asarray(tt)/2*np.log(np.pi)
tt = np.linspace(5, g[-1]+20, 400000)
NN = th_of(tt)/np.pi + 1
weyl = np.interp(k - 0.5, NN, tt)
sp = np.gradient(g)
osc = []
for p in primerange(2, 50000):
    lp = np.log(p); m = 1
    while p**m < 50000:
        osc.append((lp*m, 1.0/(np.pi*np.sqrt(p**m)*m))); m += 1
Ls = np.array([o[0] for o in osc]); As = np.array([o[1] for o in osc])
def S_and_Sp(t, delta):
    w = As*np.exp(-0.5*(Ls*delta)**2)
    keep = w > 1e-8
    S = np.zeros_like(t); Sp = np.zeros_like(t)
    for L, a in zip(Ls[keep], w[keep]):
        S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
    return S, Sp
t = weyl.copy()
best = (1.0, None)
for delta in (0.5, 0.35, 0.25, 0.18, 0.12, 0.08, 0.05):
    for it in range(3):
        S, Sp = S_and_Sp(t, delta)
        thp = 0.5*np.log(t/(2*np.pi))/np.pi
        J = thp + Sp
        J = np.where(J < 0.3*thp, 0.3*thp, J)      # trust region: floor the Jacobian
        resid = (th_of(t)/np.pi + 1 + S) - (k - 0.5)
        step = np.clip(resid/J, -1.5, 1.5)
        t = t - step
    rms = np.sqrt(np.mean(((t-g)/sp)**2))
    print(f"delta={delta:.2f}: rms = {rms:.4f} spacings", flush=True)
    if rms < best[0]: best = (rms, delta)
    else: pass
print("best: rms %.4f at delta %.2f" % best)
print("DONE")
