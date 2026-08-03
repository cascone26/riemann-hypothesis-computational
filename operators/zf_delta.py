"""Zero-free reconstruction for Delta (degree 2): oscillators from Lambda_Delta(n) =
(alpha^m+beta^m) log p (analytic norm: alpha beta = 1, alpha+beta = lam(p));
Landau: <e^{i gamma log x}> = -(T/2piN) a_Lam(x)/sqrt(x)."""
import numpy as np
from scipy.special import jv
from sympy import primerange, factorint
from itertools import combinations, product
zeros = np.load("Delta_zeros_600.npy"); ph = np.load("Delta_argLp_600.npy")
N = len(zeros); T = float(zeros[-1])
land_c = T/(2*np.pi*N)
poly = [1]+[0]*200
for n in range(1,201):
    for _ in range(24):
        new = poly[:]
        for i in range(200, n-1, -1): new[i] -= poly[i-n]
        poly = new
lam = {n: poly[n-1]/n**5.5 for n in range(1,201)}
def alphabeta_sum(p, m):  # alpha^m + beta^m via Chebyshev-like recursion, alpha+beta=lam(p), ab=1
    s1, s0 = lam[p], 2.0
    if m == 0: return s0
    if m == 1: return s1
    prev, cur = s0, s1
    for _ in range(m-1):
        prev, cur = cur, lam[p]*cur - prev
    return cur
def aLam(x):
    f = factorint(x)
    if len(f) != 1: return 0.0
    p, m = list(f.items())[0]
    if p > 200: return 0.0
    return alphabeta_sum(p, m)*np.log(p)
def mean_exp(v):
    if abs(v) < 1e-9: return 1.0
    x = np.exp(abs(v)); xr = round(x)
    if abs(x-xr) < 1e-9 and xr > 1 and len(factorint(xr)) == 1:
        return -land_c*aLam(xr)/np.sqrt(xr)
    return 0.0
X = 200
osc = []
for p in primerange(2, X):
    m = 1
    while p**m < X:
        n = p**m
        c = alphabeta_sum(p, m)/(m*np.sqrt(n))   # signed amplitude
        if abs(c) > 1e-12: osc.append((n, np.log(p)*m, c))
        m += 1
B = np.prod([jv(0, abs(c)) for _,_,c in osc])
oscA = [o for o in osc if o[0] <= 40]
def G_formula(u):
    tot = B*mean_exp(u)
    for n,L,c in osc:
        Jr = B/jv(0,abs(c))
        for k in (1,-1,2,-2):
            me = mean_exp(u+k*L)
            if me: tot += Jr*jv(k,abs(c))*(np.sign(c)**abs(k))*me
    for (n1,L1,c1),(n2,L2,c2) in combinations(oscA,2):
        Jr = B/(jv(0,abs(c1))*jv(0,abs(c2)))
        for k1,k2 in product((1,-1),repeat=2):
            me = mean_exp(u+k1*L1+k2*L2)
            if me: tot += Jr*jv(k1,abs(c1))*np.sign(c1)*jv(k2,abs(c2))*np.sign(c2)*me
    return tot
w = np.exp(-1j*ph)
G2m = (w*np.exp(1j*np.log(2)*zeros)).mean()
rot = np.exp(-1j*np.angle(G2m))  # measured rotation anchor
print(" n   b_true     formula    measured(rot)")
bs = {2:-lam[2], 3:-lam[3], 4:1.0, 5:-lam[5], 6:lam[2]*lam[3], 7:-lam[7], 9:1.0}
for n in (2,3,4,5,6,7,9):
    f = G_formula(np.log(n))
    meas = ((w*np.exp(1j*np.log(n)*zeros)).mean()*rot).real
    print(f"{n:>3}  {bs[n]:+.3f}    {f:+.4f}    {meas:+.4f}")
print("DONE")
