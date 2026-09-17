import mpmath as mp, numpy as np
mp.mp.dps=25
D=23; sc=mp.mpf(2)/mp.sqrt(D)   # normalize lattice to covolume 1: q_m = Q_m * 2/sqrt(D)
R=20; QM=[]
for m in range(-R,R+1):
    for n in range(-R,R+1):
        if m==0 and n==0: continue
        QM.append(m*m+m*n+6*n*n)
QM=np.array(sorted(QM)); QM=QM[QM<=120]
qm=[sc*int(q) for q in QM]
def Lt(s):  # Lambda-tilde, self-dual: Lt(s)=Lt(1-s)
    s=mp.mpc(s); t=mp.mpc(0)
    for q in qm:
        x=mp.pi*q
        t+=x**(-s)*mp.gammainc(s,x)+x**(s-1)*mp.gammainc(1-s,x)
    return t - 1/s - 1/(1-s)
# Z_Q(s) = (2/sqrtD)^s pi^s Lt(s)/Gamma(s); zeros of Z_Q = zeros of Lt (prefactor zero-free)
def ZQ(s): s=mp.mpc(s); return (sc**s)*(mp.pi**s)*Lt(s)/mp.gamma(s)
direct=sum(mp.mpf(int(q))**(-2) for q in QM)
print(f"Z_Q(2): continuation={mp.nstr(ZQ(2),9)}  direct-lattice={mp.nstr(direct,9)}  match={abs(complex(ZQ(2))-complex(direct))<1e-3}")
if abs(complex(ZQ(2))-complex(direct))<1e-3:
    ts=np.arange(1,45,0.3); av=[float(abs(Lt(mp.mpc(1.05,t)))) for t in ts]; found=[]
    for i in range(1,len(av)-1):
        if av[i]<av[i-1] and av[i]<av[i+1] and av[i]<0.3:
            try:
                z=mp.findroot(Lt,mp.mpc(1.05,ts[i]))
                if abs(Lt(z))<1e-8 and mp.re(z)>0.6:
                    zc=complex(z)
                    if all(abs(zc-w)>1e-2 for w in found): found.append(zc)
            except: pass
    print(f"OFF-LINE zeros of Z_Q (disc -23, h=3) found (Re>0.6): {len(found)}")
    for z in found[:6]: print(f"   s={z.real:.4f}+{z.imag:.3f}i  (off-line by {z.real-0.5:+.4f})")

print("\nFE check Lt(s)=Lt(1-s) [truncation-independent validation]:")
for s in [mp.mpc(0.3,5), mp.mpc(0.7,12), mp.mpc(0.55,20)]:
    a=Lt(s); b=Lt(1-s)
    print(f"  Lt({mp.nstr(s,4)})={mp.nstr(a,6)}  Lt(1-s)={mp.nstr(b,6)}  match={abs(complex(a)-complex(b))<1e-4}")
print("\nsigma=1.05 off-line zero scan (formula validated by FE):")
ts=np.arange(1,45,0.3); av=[float(abs(Lt(mp.mpc(1.05,t)))) for t in ts]; found=[]
for i in range(1,len(av)-1):
    if av[i]<av[i-1] and av[i]<av[i+1] and av[i]<0.3:
        try:
            z=mp.findroot(Lt,mp.mpc(1.05,ts[i]))
            if abs(Lt(z))<1e-7 and mp.re(z)>0.6:
                zc=complex(z)
                if all(abs(zc-w)>1e-2 for w in found): found.append(zc)
        except: pass
print(f"OFF-LINE zeros of Z_Q (disc -23, h=3) found (Re>0.6): {len(found)}")
for z in found[:6]: print(f"   s={z.real:.4f}+{z.imag:.3f}i  (off-line by {z.real-0.5:+.4f})")
