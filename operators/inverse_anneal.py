"""Annealed inverse map: solve theta/pi + 1 + S_delta(t_k) = k-1/2 with delta-schedule,
local warm-started solves. Target: beat the 0.059-spacing plateau of LEDGER 623."""
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
# oscillator bank to X=50000
osc = []
for p in primerange(2, 50000):
    lp = np.log(p); m = 1
    while p**m < 50000:
        osc.append((np.log(p)*m, 1.0/(np.pi*np.sqrt(p**m)*m))); m += 1
Ls = np.array([o[0] for o in osc]); As = np.array([o[1] for o in osc])
def S_delta(t, delta):
    S = np.zeros_like(t)
    w = As*np.exp(-0.5*(Ls*delta)**2)
    keep = w > 1e-8
    for L, a in zip(Ls[keep], w[keep]):
        S += -a*np.sin(t*L)
    return S
t = weyl.copy()
print("start rms: %.4f spacings" % np.sqrt(np.mean(((t-g)/sp)**2)), flush=True)
for delta in (0.5, 0.35, 0.25, 0.18, 0.12, 0.08):
    for it in range(2):
        S = S_delta(t, delta)
        # local solve: theta(t)/pi + 1 + S = k - 1/2  ->  t_new via one Newton step
        # d/dt [theta/pi + S] = theta'/pi + S'; approximate S' small vs theta' for stability
        thp = 0.5*np.log(t/(2*np.pi))/np.pi
        resid = (th_of(t)/np.pi + 1 + S) - (k - 0.5)
        t = t - np.clip(resid/thp, -2.5, 2.5)
    rms = np.sqrt(np.mean(((t-g)/sp)**2))
    print(f"delta={delta:.2f}: rms = {rms:.4f} spacings", flush=True)
print("DONE")
