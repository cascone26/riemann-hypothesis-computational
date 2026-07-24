# NEAR-THRESHOLD EXPANSION (lam^2 -> 2+): how does arch-min descend, prime-rescue grow, min(tau) behave?
# Goal: leading-order law. If min(tau) >= 0 provably via prime cancellation >= arch descent.
import mpmath as mp
mp.mp.dps=35
def build(LSQ):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    N=10; idx=list(range(-N,N+1)); D=len(idx)
    Rq=TWO_PI*N/L+50; G=1200
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
    # one prime p=2 only:
    lg2=mp.log(2); w2=lg2/mp.sqrt(2)  # Lambda(2)/sqrt(2) = log2/sqrt2
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    def Pnm(n,m): return w2*qk(n,m,lg2)  # single term k=2
    AB=mp.matrix(D,D); T=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j in range(i,D):
            m=idx[j]; ab=arch(i,j)+Bnm(n,m); pp=Pnm(n,m) if LSQ>=2 else mp.mpf(0)
            AB[i,j]=ab;AB[j,i]=ab; T[i,j]=ab-pp;T[j,i]=ab-pp
    return AB,T,D,L,w2
print(f"{'lam^2':>7} {'eps=lam^2-2':>11} {'min(AB)':>13} {'min(tau)':>13} {'window(1-log2/L)':>16}",flush=True)
for eps in [mp.mpf('0.001'),mp.mpf('0.01'),mp.mpf('0.05'),mp.mpf('0.1'),mp.mpf('0.2'),mp.mpf('0.4'),mp.mpf('0.7')]:
    LSQ=2+eps
    AB,T,D,L,w2=build(LSQ)
    mab=min(mp.eigsy(AB,eigvals_only=True))
    mt=min(mp.eigsy(T,eigvals_only=True))
    window=1-mp.log(2)/L
    print(f"{float(LSQ):7.3f} {float(eps):11.4f} {mp.nstr(mab,5):>13} {mp.nstr(mt,5):>13} {mp.nstr(window,5):>16}",flush=True)
print("\nRead: does min(tau) stay >=0 while min(AB)<0? Is min(tau) ~ O(eps^2) (tangent) or bounded below >0?",flush=True)
