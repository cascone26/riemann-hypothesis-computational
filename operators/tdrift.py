import numpy as np
from scipy.special import loggamma as slg
from sympy import primerange
g_all = np.load("zeros10k.npy")
th_of = lambda tt: np.imag(slg(0.25+0.5j*np.asarray(tt))) - np.asarray(tt)/2*np.log(np.pi)
osc = []
for p in primerange(2, 50000):
    lp = np.log(p); m = 1
    while p**m < 50000:
        osc.append((lp*m, 1.0/(np.pi*np.sqrt(p**m)*m))); m += 1
Ls = np.array([o[0] for o in osc]); As = np.array([o[1] for o in osc])
def run(g, k0):
    N = len(g); k = np.arange(1, N+1) + k0
    tt = np.linspace(g[0]-5, g[-1]+5, 150000)
    t = np.interp(k - 0.5, th_of(tt)/np.pi + 1, tt)
    sp = np.gradient(g)
    for delta in (0.5, 0.35, 0.25, 0.18, 0.12):
        w = As*np.exp(-0.5*(Ls*delta)**2); keep = w > 1e-8
        for _ in range(3):
            S = np.zeros_like(t); Sp = np.zeros_like(t)
            for L, a in zip(Ls[keep], w[keep]):
                S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
            thp = 0.5*np.log(t/(2*np.pi))/np.pi
            J = np.maximum(thp + Sp, 0.3*thp)
            t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/J, -1.5, 1.5)
    return np.sqrt(np.mean(((t-g)/sp)**2))
lx2 = np.log(50000)**2
print("windows: const = rms*(logX)^2 vs log(T/2pi):")
for lo, hi in ((100,1600),(1600,3200),(3200,4800),(4800,6900),(6900,9900)):
    gz = g_all[lo:hi]
    r = run(gz, lo)
    Tm = 0.5*(gz[0]+gz[-1])
    print(f"  Tmid={Tm:6.0f} log(T/2pi)={np.log(Tm/6.283):.2f}: const={r*lx2:.2f}", flush=True)
print("DONE")
