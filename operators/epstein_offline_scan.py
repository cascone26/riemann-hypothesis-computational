import mpmath as mp, numpy as np
mp.mp.dps=22
def Lchi(s,q,chi): return mp.power(q,-s)*mp.fsum(chi(a)*mp.zeta(s,mp.mpf(a)/q) for a in range(1,q+1))
def chi_m4(n): r=n%4; return 1 if r==1 else(-1 if r==3 else 0)
def chi_5(n): r=n%5; return 0 if r==0 else(1 if r in(1,4) else -1)
def chi_m20(n): return chi_m4(n)*chi_5(n)
def F(s): s=mp.mpc(s); return mp.zeta(s)*Lchi(s,20,chi_m20)+Lchi(s,4,chi_m4)*Lchi(s,5,chi_5)
# scan OFF-LINE region Re in [0.6,1.25] (definitively off-line if zero here), t in [8,70], |F| minima -> findroot
offzeros=[]
for sig in [0.65,0.8,0.95,1.1]:
    ts=np.arange(8.0,70,0.35); av=[float(abs(F(mp.mpc(sig,t)))) for t in ts]
    for i in range(1,len(av)-1):
        if av[i]<av[i-1] and av[i]<av[i+1] and av[i]<0.4:
            try:
                z=mp.findroot(F, mp.mpc(sig,ts[i]))
                if abs(F(z))<1e-8 and mp.re(z)>0.55 and 5<mp.im(z)<72:
                    zc=complex(z)
                    if all(abs(zc-w)>1e-2 for w in offzeros): offzeros.append(zc)
            except: pass
offzeros.sort(key=lambda z:z.imag)
print(f"OFF-LINE zeros of zeta_Q1 (Re>0.55) found in t<70: {len(offzeros)}")
for z in offzeros: print(f"   s = {z.real:.4f} + {z.imag:.3f}i   (off-line by {z.real-0.5:+.4f})")
np.save('/tmp/epstein_offzeros.npy', np.array(offzeros))
