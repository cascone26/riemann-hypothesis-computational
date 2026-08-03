"""Does order-3 near-resonance close the X-drift? Full W+Landau formula at X=1000
with singles + all pairs + all W-band triples."""
import numpy as np
from scipy.special import jv
from sympy import primerange, factorint
from itertools import combinations, product
g = np.load("zeros10k.npy"); N = len(g); T = float(g[-1])
land_c = T/(2*np.pi*N)
tg = np.arange(float(g[0]), T, 0.1)
thp = 0.5*np.log(tg/(2*np.pi))/np.pi
norm = np.trapezoid(thp, tg)
vg = np.arange(0, 0.8, 2e-4)
Wg = np.empty(len(vg), dtype=complex)
for i in range(0, len(vg), 200):
    sl = vg[i:i+200]
    Wg[i:i+200] = np.trapezoid(thp[None,:]*np.exp(1j*np.outer(sl, tg)), tg, axis=1)/norm
WgR, WgI = Wg.real.copy(), Wg.imag.copy()
def Wv_arr(v):
    a = np.abs(v)
    W = np.interp(a, vg, WgR) + 1j*np.interp(a, vg, WgI)
    return np.where(v >= 0, W, np.conj(W))*(a < 0.8)
def landau_scalar(v):
    x = np.exp(abs(v)); xr = round(x)
    if abs(x-xr) < 1e-9 and xr > 1 and len(factorint(xr)) == 1:
        return -land_c*np.log(list(factorint(xr))[0])/np.sqrt(xr)
    return 0.0
X = 3000
osc = []
for p in primerange(2, X):
    m = 1
    while p**m < X:
        osc.append((p**m, np.log(p)*m, 1.0/(m*np.sqrt(p**m)))); m += 1
M = len(osc)
Ls = np.array([o[1] for o in osc]); cs = np.array([o[2] for o in osc])
J0 = jv(0, cs); J1 = jv(1, cs)
B = np.prod(J0)
u = np.log(2)
tot = B*(Wv_arr(np.array([u]))[0] + landau_scalar(u))
# singles
for i in range(M):
    for k in (1,-1,2,-2):
        v = u + k*Ls[i]
        me = Wv_arr(np.array([v]))[0] + landau_scalar(v)
        if abs(me) > 1e-9:
            tot += B/J0[i]*jv(k, cs[i])*me
# pairs (vectorized over j)
for i in range(M):
    for s1, s2 in product((1,-1), repeat=2):
        vv = u + s1*Ls[i] + s2*Ls[i+1:]
        me = Wv_arr(vv)
        # add Landau at exact prime powers
        for jj in np.where(np.abs(me) > 0)[0] if False else []:
            pass
        # Landau hits: check integer frequencies
        xs = np.exp(np.abs(vv)); xr = np.round(xs)
        hit = np.abs(xs - xr) < 1e-9
        for jj in np.where(hit)[0]:
            me[jj] += landau_scalar(vv[jj])
        contrib = B/(J0[i]*J0[i+1:])*jv(s1, cs[i])*jv(s2, cs[i+1:])*me
        tot += contrib.sum()
print("X=1000 order<=2 with full W-band pairs: %+.5f" % tot.real, flush=True)
# triples: W-band only (|v|<0.8), vectorize over 3rd index
tot3 = 0.0
for i in range(M):
    for j in range(i+1, M):
        pref2 = B/(J0[i]*J0[j])
        for s1, s2, s3 in product((1,-1), repeat=3):
            base = u + s1*Ls[i] + s2*Ls[j]
            vv = base + s3*Ls[j+1:]
            a = np.abs(vv)
            mask = a < 0.8
            if not mask.any(): continue
            me = Wv_arr(vv[mask])
            contrib = pref2/(J0[j+1:][mask])*jv(s1,cs[i])*jv(s2,cs[j])*jv(s3,cs[j+1:][mask])*me
            tot3 += contrib.sum()
print("order-3 W-band triple total: %+.5f%+.5fi" % (tot3.real, tot3.imag))
print("grand total: %+.5f  (true -0.335)" % (tot.real + tot3.real))
print("DONE")
