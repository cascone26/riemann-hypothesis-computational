"""
COORDINATE-FREE step-(b) test: the finite Weil minimizer θ_x (lowest eigenvector of the screw-g form on
[-a,a]) should, IF Connes' θ_x→k_λ holds, have a Mellin/Fourier transform whose real-axis zeros approach
the ζ-zeros γ_n (since M[k_λ]=Ξ, whose zeros are the γ).  Compute Z(t)=∫θ(x) e^{(1/2+it)x} dx on the window
and locate |Z(t)| minima; compare to γ_n as a grows.  Also try the plain cosine transform (no e^{x/2}) as a
control.  If the |Z| minima sharpen toward γ_n with a, step(b) is confirmed with NO basis/coordinate choice.
"""
import numpy as np
import mpmath as mp
from suzuki_screw_lambda_a import screw_g, vonmangoldt_upto

def weil_min(a, M):
    Nmax = int(np.floor(np.exp(2 * a))) + 1
    Lam = vonmangoldt_upto(min(Nmax, 30_000_000))
    x = np.linspace(-a, a, M + 1); hstep = x[1] - x[0]
    diffs = np.arange(-M, M + 1) * hstep
    gv = np.array([screw_g(d, Lam, Nmax) for d in diffs])
    idx = np.arange(M + 1); G = gv[(idx[:, None] - idx[None, :]) + M]
    w = np.sqrt(np.full(M + 1, hstep)); w[0] = w[-1] = np.sqrt(0.5 * hstep)
    Gs = w[:, None] * G * w[None, :]; Gs = (Gs + Gs.T) / 2
    ev, V = np.linalg.eigh(Gs); e = V[:, 0] / w
    return x, hstep, e / np.sqrt(np.sum(e * e) * hstep), ev[0]

def local_minima(t, y):
    mins = []
    for i in range(1, len(y) - 1):
        if y[i] < y[i - 1] and y[i] < y[i + 1]:
            mins.append(t[i])
    return np.array(mins)

GAMMA = np.array([float(mp.im(mp.zetazero(n))) for n in range(1, 13)])

print("ζ-zeros γ_1..6:", np.round(GAMMA[:6], 3))
print("\n|Z(t)| minima of the Weil-minimizer transform vs ζ-zeros, per window a")
print("(Mellin: Z=∫θ e^{(1/2+it)x}dx ; control: plain ∫θ cos(tx)dx)\n")
for a in [2.0, 3.0, 4.0, 5.0]:
    M = min(3000, max(1200, int(1400 / a) * 4))
    x, hstep, theta, eps = weil_min(a, M)
    tg = np.linspace(10, 45, 7000)
    # Mellin transform on critical line
    Zm = np.array([np.abs(np.sum(theta * np.exp((0.5 + 1j * t) * x)) * hstep) for t in tg])
    # plain cosine control
    Zc = np.array([abs(np.sum(theta * np.cos(t * x)) * hstep) for t in tg])
    mins_m = local_minima(tg, Zm)
    mins_c = local_minima(tg, Zc)
    # match first 4 γ to nearest minimum
    def match(mins):
        return [f"{mins[np.argmin(np.abs(mins-g))]:.2f}" if len(mins) else "--" for g in GAMMA[:4]]
    def err(mins):
        if not len(mins): return 9.99
        return np.mean([abs(mins[np.argmin(np.abs(mins-g))]-g) for g in GAMMA[:4]])
    print(f"a={a:.1f} eps={eps:.3f}  #minima(Mellin)={len(mins_m)} #minima(cos)={len(mins_c)}")
    print(f"   γ(1..4)   = {[f'{g:.2f}' for g in GAMMA[:4]]}")
    print(f"   Mellin →  = {match(mins_m)}   mean|err|={err(mins_m):.3f}")
    print(f"   cos    →  = {match(mins_c)}   mean|err|={err(mins_c):.3f}")
