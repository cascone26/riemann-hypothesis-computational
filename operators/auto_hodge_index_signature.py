import numpy as np, math
from scipy.special import digamma

def von_mangoldt(N):
    Lam=np.zeros(N+1); sieve=np.ones(N+1,bool); sieve[:2]=False
    for p in range(2,N+1):
        if sieve[p]:
            lp=math.log(p); pk=p
            while pk<=N: Lam[pk]=lp; pk*=p
            for m in range(p*p,N+1,p): sieve[m]=False
    return Lam

N=20000; Lam=von_mangoldt(N)
ns=np.nonzero(Lam>0)[0]; logn=np.log(ns.astype(float))
weight=Lam[ns]/np.sqrt(ns.astype(float))

w=0.35; centers=np.linspace(0.3,4.8,14); J=len(centers)
Rmax=45.0; Nr=9001; r=np.linspace(-Rmax,Rmax,Nr); dr=r[1]-r[0]
Omega=digamma(0.25+1j*r/2).real-math.log(math.pi)
pref=2*w*math.sqrt(2*math.pi)
Gr=pref*np.exp(-0.5*w*w*r*r)[None,:]*np.cos(np.outer(centers,r))
A=(Gr*Omega[None,:])@Gr.T*(dr/(2*math.pi))
gpole=pref*math.exp(w*w/8.0)*np.cosh(centers/2.0)
Pole=2.0*np.outer(gpole,gpole)
cs=math.sqrt(math.pi)*w; Pr=np.zeros((J,J))
for j in range(J):
    for k in range(J):
        acc=0.0
        for s1 in (1.0,-1.0):
            for s2 in (1.0,-1.0):
                d=s1*centers[j]-s2*centers[k]-logn
                acc+=np.sum(weight*cs*np.exp(-(d*d)/(4*w*w)))
        Pr[j,k]=2.0*acc
M=Pole+A-Pr
ev=np.linalg.eigvalsh(M); Mfin=Pole-Pr; evf=np.linalg.eigvalsh(Mfin)
tol=1e-9*max(1.0,abs(ev).max())
print('Weil-form M eigenvalues min=%.4e  #neg=%d/%d'%(ev.min(),(ev<-tol).sum(),J))
print('Signature(M): npos=%d nneg=%d'%((ev>tol).sum(),(ev<-tol).sum()))
print('Finite Weil form signature: npos=%d nneg=%d'%((evf>tol).sum(),(evf<-tol).sum()))
if (evf>tol).sum()==1:
    print('HODGE_INDEX_LORENTZIAN: finite Weil form has signature (1,%d) -- matches surface intersection form'%((evf<-tol).sum()))
    defc=(-evf[evf<-tol]).sum(); ain=(ev[ev>tol]).sum()
    print('archimedean must dominate Lorentzian deficit: neg-mass(finite)=%.4e , A pos-mass=%.4e'%(defc,ain))
if ev.min()>=-tol:
    print('WEIL_FORM_PSD_ON_BASIS: min eig %.4e >=0 (consistent with RH; A fills the Lorentzian deficit)'%ev.min())
else:
    print('NEGATIVE_DIRECTION: min eig %.4e (candidate obstruction test function; probe basis widths)'%ev.min())
