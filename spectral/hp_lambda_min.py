"""High-precision (mpmath dps=40) build of the Suzuki/Weil form Q at moderate (a,K) to resolve:
is λ_min a TRUE small-positive saturation, or a float64 artifact (collapses to ~0 at high precision)?
Closed forms for pole/core/prime; tail via mp.quad of the Lerch-weighted trig."""
import mpmath as mp
mp.mp.dps=40

def build_Q_hp(a, K, NGtail=60):
    a=mp.mpf(a); L=2*a
    om=[ (k*mp.pi)/(2*a) for k in range(1,K+1)]
    def Ep(p,k):
        w=om[k-1]; return w*(mp.e**(-p*a) - ((-1)**k)*mp.e**(p*a))/(p*p+w*w)
    Ehp=[Ep(mp.mpf(1)/2,k) for k in range(1,K+1)]
    Ehm=[Ep(-mp.mpf(1)/2,k) for k in range(1,K+1)]
    # pole
    Qp=[[Ehp[i]*Ehm[j]+Ehm[i]*Ehp[j] for j in range(K)] for i in range(K)]
    # core
    Ccore=-(mp.digamma(mp.mpf(1)/4)-mp.log(mp.pi))/2
    # prime: S_kl(s) closed form
    def Smat(s,i,j):
        wi,wj=om[i],om[j]
        if (i+1+j+1)%2==1: return mp.mpf(0)
        if i==j:
            return (L-s)*mp.cos(wi*s)+mp.sin(wi*s)/wi
        return (mp.sin(wj*s)-mp.sin(wi*s))/(wi-wj)+(mp.sin(wi*s)+mp.sin(wj*s))/(wi+wj)
    # von Mangoldt up to e^L
    Nmax=int(mp.floor(mp.e**L))
    lam={}
    isc=[False]*(Nmax+1)
    for p in range(2,Nmax+1):
        if not isc[p]:
            for m in range(2*p,Nmax+1,p): isc[m]=True
            pk=p
            while pk<=Nmax:
                lam[pk]=mp.log(p); pk*=p
    # tail: Psi_kl(t) and integrand
    def Psi(t,i,j):
        wi,wj=om[i],om[j]
        if (i+1+j+1)%2==1: return mp.mpf(0)
        if i==j: return wi*(L-t)*mp.sin(wi*t)
        return 2*wi*wj*(mp.cos(wi*t)-mp.cos(wj*t))/(wi*wi-wj*wj)
    def lerch1(t):
        return mp.lerchphi(mp.e**(-2*t),1,mp.mpf(1)/4)
    Q=mp.zeros(K,K)
    for i in range(K):
        for j in range(i,K):
            val=Qp[i][j]
            if i==j: val+= -2*Ccore*a
            # prime
            for n in range(2,Nmax+1):
                if n in lam and mp.log(n)<L:
                    val+= -(lam[n]/mp.sqrt(n))*Smat(mp.log(n),i,j)
            # tail: (1/2)∫_0^L e^{-t/2} Lerch(t) Psi_ij(t) dt   (log-sing at 0)
            f=lambda t,i=i,j=j: mp.e**(-t/2)*lerch1(t)*Psi(t,i,j)
            tail=mp.mpf(0)
            try: tail=mp.quad(f,[0, L/4, L/2, L])/2
            except Exception: tail=mp.quad(f,[mp.mpf('1e-6'),L])/2
            val+=tail
            Q[i,j]=val; Q[j,i]=val
    return Q

for (a,K) in [(1.0,20),(1.0,28),(1.5,24)]:
    Q=build_Q_hp(a,K)
    E,_=mp.eigsy(Q)
    ev=sorted([mp.re(x) for x in E])
    lam_min=ev[0]/a
    print(f"a={a} K={K} dps=40:  λ_min = {mp.nstr(lam_min,8)}   next: {mp.nstr(ev[1]/a,6)}, {mp.nstr(ev[2]/a,6)}")
print("Read: if λ_min stays small-POSITIVE and stable across K (not ~1e-40), it's TRUE saturation, not artifact.")
