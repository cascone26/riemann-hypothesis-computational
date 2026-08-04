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
def run(X):
    osc = []
    for p in primerange(2, X):
        lp = np.log(p); m = 1
        while p**m < X:
            osc.append((lp*m, 1.0/(np.pi*np.sqrt(p**m)*m))); m += 1
    Ls = np.array([o[0] for o in osc]); As = np.array([o[1] for o in osc])
    def S_Sp(t, delta):
        w = As*np.exp(-0.5*(Ls*delta)**2)
        keep = w > 1e-8
        S = np.zeros_like(t); Sp = np.zeros_like(t)
        for L, a in zip(Ls[keep], w[keep]):
            S += -a*np.sin(t*L); Sp += -a*L*np.cos(t*L)
        return S, Sp
    t = weyl.copy(); best = 1.0
    for delta in (0.5, 0.35, 0.25, 0.18, 0.12, 0.08, 0.05):
        for it in range(3):
            S, Sp = S_Sp(t, delta)
            thp = 0.5*np.log(t/(2*np.pi))/np.pi
            J = np.maximum(thp + Sp, 0.3*thp)
            t = t - np.clip(((th_of(t)/np.pi + 1 + S) - (k-0.5))/J, -1.5, 1.5)
        best = min(best, np.sqrt(np.mean(((t-g)/sp)**2)))
    return best
for X in (500, 5000, 50000, 500000):
    print(f"X={X:>7}: best rms = {run(X):.4f}  (1/log X = {1/np.log(X):.3f})", flush=True)
print("DONE")
