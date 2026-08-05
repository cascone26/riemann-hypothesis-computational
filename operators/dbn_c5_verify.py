"""Cross-check C_5(u0) via the DIRECT multi-index formula (paper eq for C_5) vs h-extrapolation,
and pin u0* (sign change). Two independent methods agreeing = solid basis for certification."""
import mpmath as mp
from itertools import product
mp.mp.dps = 50
def Phi(u):
    u=mp.mpf(u); tot=mp.mpf(0)
    e4u=mp.e**(4*u); e9u=mp.e**(9*u); e5u=mp.e**(5*u)
    for n in range(1,80):
        term=(2*mp.pi**2*n**4*e9u-3*mp.pi*n**2*e5u)*mp.e**(-mp.pi*n**2*e4u)
        tot+=term
        if abs(term)<mp.mpf(10)**(-45) and n>3: break
    return tot
def K(u): return Phi(abs(mp.mpf(u)))
# Taylor coeffs a_k(u0)=K^{(k)}(u0)/k! for k=0..14 (u0>0 so K=Phi, smooth)
def a_coeffs(u0, kmax=14):
    return [mp.diff(lambda x: Phi(x), mp.mpf(u0), k)/mp.factorial(k) for k in range(kmax+1)]
# W(k0..k4)=det[(i-j)^{k_i}]_{i,j=0..4}  (generalized Vandermonde; (i-j)^0=1, incl (i-j)=0 -> 0^0=1)
def Wgv(ks):
    M=mp.matrix(5,5)
    for i in range(5):
        for j in range(5):
            d=i-j
            M[i,j]=mp.mpf(1) if ks[i]==0 else mp.mpf(d)**ks[i]
    return mp.det(M)
def C5_direct(u0):
    a=a_coeffs(u0)
    tot=mp.mpf(0)
    # sum over k0+..+k4=20, 0<=ki<=14
    for ks in product(range(15),repeat=5):
        if sum(ks)!=20: continue
        w=Wgv(ks)
        if w==0: continue
        tot+=mp.fprod([a[k] for k in ks])*w
    return tot
def D5(u0,h):
    M=mp.matrix(5,5)
    for i in range(5):
        for j in range(5): M[i,j]=K(mp.mpf(u0)+(i-j)*mp.mpf(h))
    return mp.det(M)
def C5_extrap(u0):
    hs=[mp.mpf('0.003'),mp.mpf('0.002')]
    v=[D5(u0,h)/h**20 for h in hs]
    b=(v[0]-v[1])/(hs[0]**2-hs[1]**2)
    return v[1]-b*hs[1]**2
for u0 in ('0.010','0.031','0.032'):
    cd=C5_direct(mp.mpf(u0)); ce=C5_extrap(mp.mpf(u0))
    print(f"u0={u0}: C5_direct={mp.nstr(cd,8)}  C5_extrap={mp.nstr(ce,8)}  rel.diff={mp.nstr(abs((cd-ce)/cd),3)}", flush=True)
# pin u0* by bisection on C5_direct
lo,hi=mp.mpf('0.031'),mp.mpf('0.032')
for _ in range(30):
    mid=(lo+hi)/2
    if mp.sign(C5_direct(mid))==mp.sign(C5_direct(lo)): lo=mid
    else: hi=mid
print(f"u0* (C5 sign change) = {mp.nstr((lo+hi)/2, 12)}")
print("DONE")
