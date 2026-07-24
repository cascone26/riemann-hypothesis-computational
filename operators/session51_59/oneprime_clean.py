# ONE-PRIME (p=2) DECOMPOSITION: does prime p=2 CANCEL the archimedean negative direction?
# tau = AB - P,  AB = arch(A)+pole(B),  P = prime term (only k=2 active when 2<=lam^2<3).
# Compute: min eig of AB, min eig of full tau, and Rayleigh quotient of tau, AB, P on the
# arch-min eigenvector u0. If <u0|P|u0> <= <u0|AB|u0> (P even more negative there) => -P rescues.
import mpmath as mp
mp.mp.dps=35
def build(LSQ):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    N=10; idx=list(range(-N,N+1)); D=len(idx)
    Rq=TWO_PI*N/L+50; G=1400
    grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
    wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
    Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
    def ehat(r,n):
        s=TWO_PI*n/L; d=s-r
        if abs(d)<mp.mpf('1e-28'): return mp.sqrt(L)
        return (1/mp.sqrt(L))*mp.e**(0.5j*d*L)*2*mp.sin(0.5*d*L)/d
    E=[[ehat(r,n) for r in grid] for n in idx]
    def arch(i,j): return mp.re(mp.fsum(E[i][g]*mp.conj(E[j][g])*wq[g]*Om[g] for g in range(G+1)))/(2*mp.pi)
    sh=mp.sinh(L/4)**2
    def Bnm(n,m):
        num=L**2-16*mp.pi**2*m*n; den=(L**2+16*mp.pi**2*m**2)*(L**2+16*mp.pi**2*n**2); return 32*L*sh*num/den
    def vm(k):
        n=k;p=0;f=2;ok=True
        while f*f<=n:
            if n%f==0:
                if p==0:p=f
                elif p!=f:ok=False;break
                while n%f==0:n//=f
            f+=1
        if n>1:
            if p==0:p=n
            elif p!=n:ok=False
        return mp.log(p) if(ok and p>0) else mp.mpf(0)
    VM=[(k,vm(k)) for k in range(2,int(mp.floor(LSQ))+1) if vm(k)>0]
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    def Pnm(n,m): return mp.fsum(lp*mp.mpf(k)**mp.mpf('-0.5')*qk(n,m,mp.log(k)) for k,lp in VM)
    AB=mp.matrix(D,D); P=mp.matrix(D,D); T=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j in range(i,D):
            m=idx[j]; ab=arch(i,j)+Bnm(n,m); pp=Pnm(n,m)
            AB[i,j]=ab;AB[j,i]=ab; P[i,j]=pp;P[j,i]=pp; T[i,j]=ab-pp;T[j,i]=ab-pp
    return AB,P,T,D,VM
def rq(Mx,v,D):  # Rayleigh quotient v^T M v / v^T v
    num=mp.fsum(v[i]*Mx[i,j]*v[j] for i in range(D) for j in range(D))
    den=mp.fsum(v[i]*v[i] for i in range(D)); return num/den
print(f"{'lam^2':>6} {'min(AB)':>12} {'min(tau)':>12} | on arch-min-vec u0: {'<u0|AB|u0>':>11} {'<u0|P|u0>':>11} {'<u0|tau|u0>':>12}",flush=True)
for LSQ in [mp.mpf('2.0'),mp.mpf('2.2'),mp.mpf('2.5'),mp.mpf('2.72')]:
    AB,P,T,D,VM=build(LSQ)
    Eab,Qab=mp.eigsy(AB); oab=sorted(range(D),key=lambda k:Eab[k])
    Et,Qt=mp.eigsy(T); ot=sorted(range(D),key=lambda k:Et[k])
    u0=[Qab[i,oab[0]] for i in range(D)]  # arch-min eigenvector
    print(f"{float(LSQ):6.2f} {mp.nstr(Eab[oab[0]],5):>12} {mp.nstr(Et[ot[0]],5):>12} |            "
          f"{mp.nstr(rq(AB,u0,D),4):>11} {mp.nstr(rq(P,u0,D),4):>11} {mp.nstr(rq(T,u0,D),4):>12}  (nprime={len(VM)})",flush=True)
print("\nIf min(tau)>min(AB) and <u0|P|u0> <= <u0|AB|u0>: prime p=2 CANCELS arch-negativity => rescues positivity.",flush=True)
