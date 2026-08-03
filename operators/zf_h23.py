import numpy as np
from scipy.special import jv
from sympy import primerange, factorint
from itertools import combinations, product
def ks(a, n):  # kronecker symbol via sympy fallback
    from sympy.ntheory.residue_ntheory import jacobi_symbol
    return jacobi_symbol(a % n if n % 2 else a, n) if n % 2 else 0
an = np.loadtxt("h23_an.txt")
zeros = np.array([float(l.replace(" E","e")) for l in open("h23_zeros.txt") if l.strip()])
ph = np.array([float(l.replace(" E","e")) for l in open("h23_argLp.txt") if l.strip()])
N = len(zeros); T = float(zeros[-1])
land_c = T/(2*np.pi*N)
from sympy import jacobi_symbol as jac
def chi(p):
    if p == 23: return 0
    if p == 2: return 1 if (-23) % 8 in (1, 7) else -1
    return int(jac(-23, p))
def satake_sum(p, m):
    s0, s1 = 2.0, float(an[p-1])
    if m == 0: return s0
    if m == 1: return s1
    prev, cur = s0, s1
    for _ in range(m-1):
        prev, cur = cur, float(an[p-1])*cur - chi(p)*prev
    return cur
def aLam(x):
    f = factorint(x)
    if len(f) != 1: return 0.0
    p, m = list(f.items())[0]
    if p > 200: return 0.0
    return satake_sum(p, m)*np.log(p)
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
        c = satake_sum(p, m)/(m*np.sqrt(n))
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
rot = np.exp(-1j*np.angle((w*np.exp(1j*np.log(2)*zeros)).mean()))
print(" n    formula    measured   (channel)")
for n, tag in [(2,"npr"),(3,"npr"),(4,"npr p2"),(13,"npr"),(25,"inert p2"),(59,"PRINCIPAL"),(6,"npr*npr")]:
    f = G_formula(np.log(n))
    meas = ((w*np.exp(1j*np.log(n)*zeros)).mean()*rot).real
    print(f"{n:>4}  {f:+.4f}    {meas:+.4f}   {tag}")
print("DONE")
