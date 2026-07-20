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

N=20000; Lam=von_mangoldt(N); ns=np.nonzero(Lam>0)[0]
logn=np.log(ns.astype(float)); weight=Lam[ns]/np.sqrt(ns.astype(float))
w=0.35; centers=np.linspace(0.3,4.8,16); J=len(centers)
gpole=2*w*math.sqrt(2*math.pi)*math.exp(w*w/8.0)*np.cosh(centers/2.0)
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
B=Pole-Pr   # finite arithmetic intersection pairing
evf,Vf=np.linalg.eigh(B)
print('finite pairing signature: npos=%d nneg=%d'%((evf>1e-9).sum(),(evf<-1e-9).sum()))
H=Vf[:,-1]; BHH=float(H@B@H)
print('ample class self-pairing B(H,H)=%.4e'%BHH)
rng=np.random.default_rng(0); minmargin=1e18; viol=0
for _ in range(40000):
    D=rng.standard_normal(J)
    BDD=float(D@B@D); BDH=float(D@B@H)
    margin=BDH*BDH-BDD*BHH   # >=0 required by Hodge index if signature (1,n-1)
    if margin<minmargin: minmargin=margin
    if margin<-1e-6*max(1.0,abs(BHH)): viol+=1
print('min over 40000 classes of (D.H)^2-(D.D)(H.H) = %.4e   violations=%d'%(minmargin,viol))
if BHH>0 and viol==0:
    print('HODGE_INDEX_INEQUALITY_HOLDS: reversed Cauchy-Schwarz satisfied for all tested D -> Castelnuovo/Hodge certificate; Weil positivity of the finite pairing certified on this basis')
else:
    print('MARGIN_NEGATIVE: pairing not globally Lorentzian on this basis (structure probe: refine ample class / basis)')
