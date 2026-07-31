import numpy as np
from flint import arb, ctx
from math import gcd as _gcd
ctx.dps=25
DIV=[1,2,3,5,6,10,15,30]; PRIMES=[2,3,5]
G30=[_gcd(r,30) for r in range(30)]; DIVS_OF={g:[d for d in DIV if g%d==0] for g in set(G30)}
def Flo(t_str,y_str,Nm,Np,EXP,FAC=1.03):
    t=arb(t_str);y=arb(y_str);qt=t/4
    bp={p:(qt*arb(p).log()**2).exp() for p in PRIMES}; Lam={}
    for d in DIV:
        v=arb(1)
        for p in PRIMES:
            if d%p==0: v=v*(-bp[p])
        Lam[d]=v
    xN=4*arb.pi()*arb(Nm)**2-arb.pi()*t/4; sig=(1+y)/2+qt*(xN/(4*arb.pi())).log()
    c=arb(Nm)**(arb(str(-EXP))); a1=c*arb(str(FAC)); r=(1-a1)/(1+a1)
    DN=30*Np; SUM=arb(0); p02=arb(str(EXP)); fF=arb(str(FAC)); bc={}
    def bm(m):
        v=bc.get(m)
        if v is None:
            v=(qt*arb(m).log()**2).exp()
            if len(bc)<2000000: bc[m]=v
        return v
    for n in range(2,DN+1):
        ds=DIVS_OF[G30[n%30]]; bn_=None;an_=None
        for d in ds:
            m=n//d
            if n==d*m and m<=Np:
                b=bm(m); bn_=(Lam[d]*b if bn_ is None else bn_+Lam[d]*b)
                an_=(Lam[d]*(arb(m)**p02)*b if an_ is None else an_+Lam[d]*(arb(m)**p02)*b)
        if bn_ is None: continue
        an_=an_*fF; t1=abs(bn_-c*an_); t2=r*abs(bn_+c*an_); SUM=SUM+(t1 if t1>t2 else t2)*(arb(n)**(-sig))
    return float((1-a1-SUM).lower())
print("EXP scan at FAC=1.03 (safe), binding tile, checking if smaller EXP crosses F_lo>0 past the floor.")
print("First: does EXP improve at the INHERITED certified target (0.170,0.160)? (validity sanity: must stay >0)")
for EXP in [0.20,0.15,0.10,0.07,0.05]:
    v=Flo('0.170','0.160',690988,800000,EXP)
    print(f"  (0.170,0.160) EXP={EXP:.2f}: F_lo={v:.5f}",flush=True)
print("\nNow the target past the floor (0.163,0.152, Lambda=0.17455):")
best=(-9,None)
for EXP in [0.20,0.15,0.10,0.07,0.05,0.03,0.02]:
    v=Flo('0.163','0.152',690988,800000,EXP)
    if v>best[0]: best=(v,EXP)
    print(f"  EXP={EXP:.2f}: F_lo={v:.5f} {'>0!' if v>2.15e-3 else ''}",flush=True)
print(f"\nBEST at this target: F_lo={best[0]:.5f} at EXP={best[1]}")
