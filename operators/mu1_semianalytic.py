import mpmath as mp
mp.mp.dps=40
# Q_jl = ∫∫ |x-y|log|x-y| vj'(x) vl'(y) dxdy, vj'(x)=(jπ/2)cos(jπ(x+1)/2) on (-1,1)
# = ∫_{-2}^{2} k(u) R_jl(u) du,  k=|u|log|u|,  R_jl(u)=∫ vj'(x) vl'(x-u) dx (overlap)
# Compute R_jl(u) in closed form; then 1D adaptive quad (handles log-sing at 0). Machine precision.
def R(j,l,u):
    # u in (0,2): overlap x+1=w in (u,2). integrand (jπ/2)(lπ/2) cos(jπw/2)cos(lπ(w-u)/2)
    pj=mp.pi*j/2; pl=mp.pi*l/2; amp=pj*pl
    # ½[cos((j-l)πw/2 + lπu/2) + cos((j+l)πw/2 - lπu/2)]
    def integ(a,b):  # ∫_u^2 cos(a w + b) dw
        if abs(a)<mp.mpf('1e-30'): return (2-u)*mp.cos(b)
        return (mp.sin(a*2+b)-mp.sin(a*u+b))/a
    a1=mp.pi*(j-l)/2; b1=pl*u
    a2=mp.pi*(j+l)/2; b2=-pl*u
    return amp*mp.mpf('0.5')*(integ(a1,b1)+integ(a2,b2))
def k(u): 
    au=abs(u); return au*mp.log(au) if au>0 else mp.mpf(0)
def Qelem(j,l):
    # Q_jl = ∫_0^2 k(u)(R_jl(u)+R_jl(-u)) du ; R_jl(-u)=R_lj(u) by symmetry of the correlation
    f=lambda u: k(u)*(R(j,l,u)+R(l,j,u))
    return mp.quad(f,[0,mp.mpf('0.0001'),2])   # split near log-sing
K=24
print(f"semi-analytic λ_min^(1), K={K} (machine-precision matrix elements):",flush=True)
Q=mp.matrix(K,K)
for j in range(1,K+1):
    for l in range(j,K+1):
        v=Qelem(j,l); Q[j-1,l-1]=v; Q[l-1,j-1]=v
# ||v||^2 = sum c^2 (mass=I). min eigenvalue of Q
E=mp.eigsy(Q,eigvals_only=True)
lam=min(E)
print(f"  lambda_min^(1) = {mp.nstr(lam,15)}",flush=True)
print(f"  mu1 = 1 + lam/2 = {mp.nstr(1+lam/2,15)}",flush=True)
print(f"  identify(mu1): {mp.identify(1+lam/2,['euler','log(2)','log(pi)','pi','catalan','zeta(3)'])}",flush=True)
print(f"  identify(lam): {mp.identify(lam,['euler','log(2)','log(pi)','pi'])}",flush=True)
