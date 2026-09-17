import mpmath as mp
mp.mp.dps=20
chiv={1:1,2:1j,3:-1j,4:-1,0:0}
def chi(n): return chiv[n%5]
def L(s,conj=False):
    tot=mp.mpc(0)
    for a in range(1,5):
        c=mp.conj(chi(a)) if conj else chi(a); tot+=c*mp.zeta(s,mp.mpf(a)/5)
    return 5**(-s)*tot
kappa=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
def F(s): return (1-1j*kappa)/2*L(s)+(1+1j*kappa)/2*L(s,conj=True)
# scan for OFF-line zeros: Re in (0.5,0.95], find where |F| has a deep minimum, then refine
print("Scanning D-H off-line zeros (Re s > 0.5) by height...")
found=[]
for T in range(5,160,1):  # height bands
    for sig in [x/100 for x in range(55,95,5)]:
        s0=mp.mpc(sig,T)
        try:
            z=mp.findroot(F,s0)
            if 0.51<mp.re(z)<0.99 and 1<mp.im(z)<160 and abs(F(z))<1e-10:
                zt=(round(float(mp.re(z)),4),round(float(mp.im(z)),3))
                if all(abs(zt[1]-f[1])>0.5 for f in found): found.append(zt)
        except: pass
found.sort(key=lambda z:z[1])
print(f"off-line zeros found (Re,Im), Re>0.51: {found}")
for H in [50,70,90,110,140]:
    print(f"  #off-line pairs with Im<={H}: {sum(1 for z in found if z[1]<=H)}")
print("Compare: our certified D-H n_neg = 1,4,8,14 at c=25,50,75,100 (cycle 215) -- Bombieri: n_neg = resolved off-line pairs.")
