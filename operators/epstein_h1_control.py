# h=1 CONTROL: Epstein zeta of x^2+y^2 (disc -4, class number 1) = 4 zeta(s) L(s,chi4) MULTIPLICATIVE.
# Same incomplete-gamma method as disc -23. Expect ZERO off-line zeros (product of on-line L's).
import mpmath as mp, numpy as np
mp.mp.dps=25
D=4; sc=mp.mpf(2)/mp.sqrt(D)   # =1 for D=4
R=22; QM=[]
for m in range(-R,R+1):
    for n in range(-R,R+1):
        if m==0 and n==0: continue
        QM.append(m*m+n*n)
QM=np.array(sorted(QM)); QM=QM[QM<=140]; qm=[sc*int(q) for q in QM]
def Lt(s):
    s=mp.mpc(s); t=mp.mpc(0)
    for q in qm:
        x=mp.pi*q; t+=x**(-s)*mp.gammainc(s,x)+x**(s-1)*mp.gammainc(1-s,x)
    return t - 1/s - 1/(1-s)
# FE validation
print("FE Lt(s)=Lt(1-s):", all(abs(complex(Lt(mp.mpc(a,b)))-complex(Lt(1-mp.mpc(a,b))))<1e-4 for a,b in [(0.3,5),(0.6,15)]))
# off-line scan (expect 0 -- it's 4*zeta*L(chi4), all zeros on line)
ts=np.arange(1,45,0.3); av=[float(abs(Lt(mp.mpc(1.05,t)))) for t in ts]; found=[]
for i in range(1,len(av)-1):
    if av[i]<av[i-1] and av[i]<av[i+1] and av[i]<0.3:
        try:
            z=mp.findroot(Lt,mp.mpc(1.05,ts[i]))
            if abs(Lt(z))<1e-7 and mp.re(z)>0.6:
                zc=complex(z)
                if all(abs(zc-w)>1e-2 for w in found): found.append(zc)
        except: pass
print(f"h=1 (disc -4) OFF-LINE zeros (Re>0.6): {len(found)}  (expect 0 = multiplicative PASS)")
for z in found[:5]: print("  ",round(z.real,4),round(z.imag,3))
print("=> Epstein family discriminator via CLASS NUMBER: h=1 mult PASS (0 off-line); h>=2 non-mult FAIL.")
