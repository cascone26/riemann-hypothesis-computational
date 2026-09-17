import sys, mpmath as mp
sys.path.insert(0,"/Users/scones/projects/riemann/operators")
from qw_lambda_faithful import build_QW
mp.mp.dps = 60
N=16; lam=float(mp.e)
Q, idx, L = build_QW(N, lam, verbose=False)
n=len(idx)
# high-precision symmetric eigensolve
E, V = mp.eigsy(mp.matrix([[Q[i,j] for j in range(n)] for i in range(n)]))
# ground state = smallest eigenvalue
order=sorted(range(n), key=lambda k: E[k])
i0=order[0]
eps1=E[order[0]]; eps2=E[order[1]]
v=[V[r,i0] for r in range(n)]
ns=[int(k) for k in idx]
vmap={ns[i]:v[i] for i in range(n)}
num=sum((vmap[k]+vmap.get(-k,0))**2 for k in ns); den=sum((2*vmap[k])**2 for k in ns) or mp.mpf(1)
evenfrac=num/den
w=[2*mp.pi*k/L for k in ns]
def F(x): return sum(v[i]*mp.cos(w[i]*x) for i in range(n))
def dF(x): return sum(-v[i]*w[i]*mp.sin(w[i]*x) for i in range(n))
def ddF(x): return sum(-v[i]*w[i]*w[i]*mp.cos(w[i]*x) for i in range(n))
mn=mp.mpf('1e9')
for j in range(400):
    x=L*j/400
    lag=dF(x)**2 - F(x)*ddF(x)
    if lag<mn: mn=lag
print(f"N={N} lam={lam:.4f} dps=60")
print(f"eps1 = {mp.nstr(eps1,6)}   eps2 = {mp.nstr(eps2,6)}   gap = {mp.nstr(eps2-eps1,4)}")
print(f"evenfrac = {mp.nstr(evenfrac,6)}   (1.0=perfectly even)")
print(f"min L[F] over grid = {mp.nstr(mn,6)}   (>=0 => Laguerre-inequality holds => consistent w/ real-rooted)")
