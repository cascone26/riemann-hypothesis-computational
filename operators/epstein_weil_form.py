import mpmath as mp, numpy as np
mp.mp.dps=22
def Lchi(s,q,chi): return mp.power(q,-s)*mp.fsum(chi(a)*mp.zeta(s,mp.mpf(a)/q) for a in range(1,q+1))
def chi_m4(n): r=n%4; return 1 if r==1 else(-1 if r==3 else 0)
def chi_5(n): r=n%5; return 0 if r==0 else(1 if r in(1,4) else -1)
def chi_m20(n): return chi_m4(n)*chi_5(n)
def F(s): s=mp.mpc(s); return mp.zeta(s)*Lchi(s,20,chi_m20)+Lchi(s,4,chi_m4)*Lchi(s,5,chi_5)
# on-line zeros: |F(0.5+it)| dips -> findroot. scan coarse, refine dips.
ts=np.arange(1.0,32,0.2); av=[float(abs(F(mp.mpc(0.5,t)))) for t in ts]
onl=[]
for i in range(1,len(av)-1):
    if av[i]<av[i-1] and av[i]<av[i+1] and av[i]<1.0:
        try:
            z=mp.findroot(F, mp.mpc(0.5, ts[i]))
            if abs(F(z))<1e-9 and abs(mp.re(z)-0.5)<1e-3:
                zc=complex(z)
                if all(abs(zc.imag-w.imag)>1e-2 for w in onl): onl.append(zc)
        except Exception: pass
# off-line quadruple (found cycle 284): both positive-t members
off=[complex(0.93767,29.9834), complex(0.06233,29.9834)]  # FE partners; Re=0.5+-0.4377
print(f"zeta_Q1: {len(onl)} on-line zeros (t<32), 1 off-line QUADRUPLE at Re=0.938/0.062, t=29.98")
SQ=np.sqrt(2*np.pi)
def phih(a,z): return SQ/2*(np.exp(-(z-a)**2/2)+np.exp(-(z+a)**2/2))
A=np.array([22.,25.,28.,30.,32.]); n=len(A)
def buildQ(zeros):
    Q=np.zeros((n,n))
    for z in zeros:
        g=z.imag; d=z.real-0.5; w=np.array([phih(a, g-1j*d) for a in A])
        Q+=2*np.real(np.outer(w,w))
    return 0.5*(Q+Q.T)
Qon=buildQ(onl); Qall=buildQ(onl+off)
evon=np.linalg.eigvalsh(Qon); evall=np.linalg.eigvalsh(Qall)
print(f"\nGram Weil form (a={A.tolist()}):")
print(f"  ON-LINE zeros only:  n_neg={(evon<-1e-9).sum()}  min-eig={evon.min():+.4e}  (PSD expected)")
print(f"  + OFF-LINE quadruple: n_neg={(evall<-1e-9).sum()}  min-eig={evall.min():+.4e}  (INDEFINITE = FAIL)")
print(f"  eigenvalues with off-line: {np.round(evall,4).tolist()}")
