# Is the one-prime term P2(n,m)=w2*q(n,m,log2) LOW RANK? If rank<=2 (or 4), one-prime positivity
# = exact low-rank perturbation interlacing of (A+B). Test singular values of P2.
import mpmath as mp
mp.mp.dps=35
for LSQ in [mp.mpf('2.2'),mp.mpf('2.5'),mp.mpf('2.72')]:
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi
    N=10; idx=list(range(-N,N+1)); D=len(idx)
    lg2=mp.log(2); w2=lg2/mp.sqrt(2)
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    P=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j,m in enumerate(idx):
            P[i,j]=w2*qk(n,m,lg2)
    # symmetric => eigenvalues; count how many are non-negligible
    ev=sorted([abs(e) for e in mp.eigsy(P,eigvals_only=True)],reverse=True)
    top=[float(e) for e in ev[:6]]
    print(f"lam^2={float(LSQ):.2f}: top |eig(P2)| = {['%.3e'%e for e in top]}  (rank = # non-tiny)",flush=True)
    # signed eigenvalues of P2 (is it rank-2 with one +, one -?)
    evs=sorted([float(e) for e in mp.eigsy(P,eigvals_only=True)],key=abs,reverse=True)
    print(f"          signed top eigs: {['%.3e'%e for e in evs[:4]]}",flush=True)
