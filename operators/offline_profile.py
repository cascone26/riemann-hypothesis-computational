"""R(sigma): how does the phase identity decay off the critical line?
For each sigma, compute arg zeta'(sigma+i gamma_k) at the first M zeros and
R(sigma) = |circ mean e^{i(argzp - pi S_k)}|. R(1/2)=1 exactly; profile shape = rigidity."""
import numpy as np, mpmath as mp
from scipy.special import loggamma as slg
mp.mp.dps = 15
g = np.load("zeros10k.npy")[:1500]
k = np.arange(1, len(g)+1)
th = np.imag(slg(0.25+0.5j*g)) - g/2*np.log(np.pi)
S = (k - 0.5) - th/np.pi - 1.0
for sig in (0.5, 0.52, 0.55, 0.6, 0.7, 0.85, 1.0, 1.25, 1.5, 2.0, 3.0):
    ph = np.array([float(mp.arg(mp.zeta(mp.mpc(sig, t), derivative=1))) for t in g])
    d = ph - np.pi*S
    z = np.exp(1j*d).mean()
    print("sigma=%.2f  R=%.4f  c=%+.3f" % (sig, abs(z), np.angle(z)), flush=True)
print("DONE")
