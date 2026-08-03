import numpy as np
from scipy.special import jv
from sympy import primerange, factorint
from itertools import combinations
g = np.load("zeros10k.npy"); N = len(g); T = float(g[-1])
land_c = T/(2*np.pi*N)
tg = np.arange(float(g[0]), T, 0.02)
thp = 0.5*np.log(tg/(2*np.pi))/np.pi
norm = np.trapezoid(thp, tg)
# W on grid |v| <= 0.6 (near-resonance band), resolution 2e-4
vg = np.arange(0, 0.6, 2e-4)
Wg = np.empty(len(vg), dtype=complex)
for i in range(0, len(vg), 200):
    sl = vg[i:i+200]
    Wg[i:i+200] = np.trapezoid(thp[None,:]*np.exp(1j*np.outer(sl, tg)), tg, axis=1)/norm
def Wv(v):
    a = abs(v)
    if a >= 0.6: return 0.0j
    W = np.interp(a, vg, Wg.real) + 1j*np.interp(a, vg, Wg.imag)
    return W if v >= 0 else np.conj(W)
def landau(v):
    x = np.exp(abs(v)); xr = round(x)
    if abs(x-xr) < 1e-9 and xr > 1 and len(factorint(xr)) == 1:
        return -land_c*np.log(list(factorint(xr))[0])/np.sqrt(xr)
    return 0.0
u = np.log(2)
for X in (200, 1000, 3000):
    osc = []
    for p in primerange(2, X):
        m = 1
        while p**m < X:
            osc.append((p**m, np.log(p)*m, 1.0/(m*np.sqrt(p**m)))); m += 1
    B = np.prod([jv(0,c) for _,_,c in osc])
    tot = B*(Wv(u)+landau(u))
    for n,L,c in osc:
        Jr = B/jv(0,c)
        for k in (1,-1,2,-2):
            v = u+k*L
            me = Wv(v)+landau(v)
            if abs(me) > 1e-9: tot += Jr*jv(k,c)*me
    Ls = np.array([o[1] for o in osc]); cs = np.array([o[2] for o in osc])
    J1s = jv(1, cs); J0s = jv(0, cs)
    M = len(osc)
    for i in range(M):
        vi = u + Ls[i]      # k_i=+1
        wj = u - Ls[i]      # k_i=-1 base handled in loop below
        # vector over j>i for both sign patterns hitting near-zero: v = u + Li - Lj and u - Li + Lj
        for s1, s2 in ((1,-1),(-1,1),(1,1),(-1,-1)):
            vv = u + s1*Ls[i] + s2*Ls[i+1:]
            mask = np.abs(vv) < 0.6
            idx = np.where(mask)[0]
            for jj in idx:
                j = i+1+jj
                me = Wv(vv[jj]) + landau(vv[jj])
                if abs(me) > 1e-9:
                    tot += B/(J0s[i]*J0s[j])*jv(s1,cs[i])*jv(s2,cs[j])*me
    print(f"X={X:>5}: formula(W+Landau, full pairs) = {tot.real:+.4f}  (true -0.335)", flush=True)
print("DONE")
