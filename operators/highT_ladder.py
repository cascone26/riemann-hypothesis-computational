import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
z = np.load("zeros100k.npy")
g = z[80000:]; k0 = 80000
N = len(g); k = np.arange(1, N+1) + k0
th_of = lambda tt: np.imag(slg(0.25+0.5j*np.asarray(tt))) - np.asarray(tt)/2*np.log(np.pi)
tt = np.linspace(g[0]-5, g[-1]+5, 300000)
t = np.interp(k - 0.5, th_of(tt)/np.pi + 1, tt)
sp = np.gradient(g)
osc = []
for p in primerange(2, 50000):
    lp = np.log(p); m = 1
    while p**m < 50000:
        osc.append((lp*m, 1.0/(np.pi*np.sqrt(p**m)*m))); m += 1
Ls = np.array([o[0] for o in osc]); As = np.array([o[1] for o in osc])
best = 1.0
for delta in (0.5, 0.35, 0.25, 0.18, 0.12, 0.085, 0.06, 0.042):
    w = As*np.exp(-0.5*(Ls*delta)**2); keep = w > 1e-8
    for _ in range(3):
        S = np.zeros_like(t); Sp = np.zeros_like(t)
        for L, a in zip(Ls[keep], w[keep]):
            S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
        thp = 0.5*np.log(t/(2*np.pi))/np.pi
        J = np.maximum(thp + Sp, 0.3*thp)
        t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/J, -1.5, 1.5)
    r = np.sqrt(np.mean(((t-g)/sp)**2))
    best = min(best, r)
    print(f"delta={delta:.3f}: rms {r:.4f} sp", flush=True)
print("best %.4f -> const %.2f at log(T/2pi)=9.29" % (best, best*np.log(50000)**2))
print("DONE")
