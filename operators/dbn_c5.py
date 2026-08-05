"""OPEN PIECE (arXiv:2602.20313, withdrawn-as-uncertified): the C_5(u0) leading coefficient
= lim_{h->0} D_5(u0,h)/h^20. Author claims C_5<0 near origin, sign change near u0~0.0311.
Compute it at high precision across u0, pin the threshold. Step toward re-certification."""
import mpmath as mp
mp.mp.dps = 60
def Phi(u):
    u=mp.mpf(u); tot=mp.mpf(0)
    e4u=mp.e**(4*u); e9u=mp.e**(9*u); e5u=mp.e**(5*u)
    for n in range(1,80):
        term=(2*mp.pi**2*n**4*e9u-3*mp.pi*n**2*e5u)*mp.e**(-mp.pi*n**2*e4u)
        tot+=term
        if abs(term)<mp.mpf(10)**(-55) and n>3: break
    return tot
def K(u): return Phi(abs(mp.mpf(u)))
def D5(u0,h):
    M=mp.matrix(5,5)
    for i in range(5):
        for j in range(5): M[i,j]=K(mp.mpf(u0)+(i-j)*mp.mpf(h))
    return mp.det(M)
def C5(u0):
    # leading coeff: D5/h^20 as h->0; Richardson extrapolate over small h
    hs=[mp.mpf('0.004'),mp.mpf('0.003'),mp.mpf('0.002')]
    vals=[D5(u0,h)/h**20 for h in hs]
    # linear extrapolation in h^2 to h=0
    h2=[h**2 for h in hs]
    # fit v = C5 + b h^2 using first two
    b=(vals[0]-vals[1])/(h2[0]-h2[1])
    return vals[1]-b*h2[1]
print("C_5(u0) across u0 (find sign change ~0.0311):")
prev=None; thresh=None
for u0 in [mp.mpf(x)/1000 for x in [1,10,20,28,30,31,32,35,40,50]]:
    c=C5(u0)
    print(f"  u0={float(u0):.3f}: C5={mp.nstr(c,8)}  sign={'+' if c>0 else '-'}", flush=True)
    if prev is not None and mp.sign(prev)!=mp.sign(c) and thresh is None:
        thresh=(prev_u,u0)
    prev=c; prev_u=u0
if thresh: print(f"sign change between u0={float(thresh[0]):.3f} and {float(thresh[1]):.3f}")
print("DONE")
