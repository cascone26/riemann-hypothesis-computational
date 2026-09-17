import mpmath as mp
mp.mp.dps=15
chiv={1:1,2:1j,3:-1j,4:-1,0:0}
def chi(n): return chiv[n%5]
def Lc(s,conj=False):
    tot=mp.mpc(0)
    for a in range(1,5):
        c=mp.conj(chi(a)) if conj else chi(a); tot+=c*mp.zeta(s,mp.mpf(a)/5)
    return 5**(-s)*tot
kap=(mp.sqrt(10-2*mp.sqrt(5))-2)/(mp.sqrt(5)-1)
def F(s): return (1-1j*kap)/2*Lc(s)+(1+1j*kap)/2*Lc(s,conj=True)
# Count zeros on the line vs off, via argument principle proxy: scan |F| on sigma=0.5 (on-line) count sign changes,
# and detect off-line by |F| minima at sigma>0.5. Fast: sample |F(0.7+iT)| and look for deep dips.
print("Off-line detection: local minima of |F(sigma+iT)| for sigma in {0.6,0.7,0.8} vs T (deep dip => off-line zero)")
offs=[]
for sig in [0.62,0.70,0.78,0.86]:
    prev2=prev=None
    T=2.0
    while T<160:
        v=abs(F(mp.mpc(sig,T)))
        if prev is not None and prev2 is not None and prev<prev2 and prev<v and prev<0.05:
            # local min below 0.05 -> refine
            try:
                z=mp.findroot(F,mp.mpc(sig,T-0.5))
                if 0.53<mp.re(z)<0.98 and abs(F(z))<1e-8:
                    zt=(round(float(mp.re(z)),3),round(float(mp.im(z)),2))
                    if all(abs(zt[1]-o[1])>0.3 for o in offs): offs.append(zt)
            except: pass
        prev2=prev;prev=v;T+=0.5
    print(f"  sigma~{sig}: scanned")
offs.sort(key=lambda z:z[1])
print(f"OFF-LINE zeros (Re,Im), Re>0.53: {offs}")
for H in [50,70,90,120,150]:
    print(f"  #off-line with Im<={H}: {sum(1 for z in offs if z[1]<=H)}")
print("our certified D-H n_neg=1,4,8,14 (c=25,50,75,100). Bombieri: n_neg = resolved off-line pairs.")
