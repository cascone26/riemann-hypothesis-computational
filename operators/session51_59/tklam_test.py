# SHARP CRITERION test: is k_lambda (k_hat = 1/4 Xi) an exact 0-EIGENVECTOR of tau (tau k=0), not just
# Rayleigh-0? Build k in the ehat basis (c_n = 1/4 Xi(2*pi*n/L) = 1/4 xi(1/2 + i 2pi n/L)), compute
# ||tau k|| / ||k|| and compare to (a) random even vectors, (b) the numeric ground state. If k is a
# near-eigenvector while random v are O(1), tau k=0 is substantiated (=> RH-consistent, sharp form).
import mpmath as mp
mp.mp.dps=30
def run(LSQ,N=12):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
    Rq=TWO_PI*N/L+50; G=1000
    grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
    wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
    Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
    def ehat(r,n):
        s=TWO_PI*n/L; d=s-r
        if abs(d)<mp.mpf('1e-24'): return mp.sqrt(L)
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
    T=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j,m in enumerate(idx): T[i,j]=arch(i,j)+Bnm(n,m)-Pnm(n,m)
    # k_lambda coefficients: c_n = 1/4 xi(1/2 + i*2*pi*n/L). xi via mpmath.
    def xi(s): return mp.mpf('0.5')*s*(s-1)*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    c=mp.matrix(D,1)
    for i,n in enumerate(idx):
        r=TWO_PI*n/L; c[i]=mp.re(mp.mpf('0.25')*xi(mp.mpf('0.5')+1j*r))
    # normalize
    nk=mp.sqrt(mp.fsum(c[i]**2 for i in range(D)))
    for i in range(D): c[i]/=nk
    Tk=T*c
    resid=mp.sqrt(mp.fsum(Tk[i]**2 for i in range(D)))  # ||tau k|| (k normalized)
    rayleigh=mp.fsum(c[i]*Tk[i] for i in range(D))       # <k|tau|k>
    # random even vectors: build a few, measure ||tau v||/||v||
    import hashlib
    def randeven(seed):
        v=mp.matrix(D,1)
        for i,n in enumerate(idx):
            h=int(hashlib.md5(f"{seed}_{abs(n)}".encode()).hexdigest(),16)%1000
            val=mp.mpf(h)/1000-mp.mpf('0.5')
            v[i]=val  # even: depends on |n|
        nv=mp.sqrt(mp.fsum(v[i]**2 for i in range(D)))
        for i in range(D): v[i]/=nv
        Tv=T*v; return float(mp.sqrt(mp.fsum(Tv[i]**2 for i in range(D))))
    rands=[randeven(s) for s in range(5)]
    gs=min(mp.eigsy(T,eigvals_only=True))
    return float(resid),float(rayleigh),rands,float(gs)
for LSQ in [mp.mpf('2.5'),mp.mpf('2.72'),mp.mpf('5.0')]:
    resid,ray,rands,gs=run(LSQ)
    print(f"lam^2={float(LSQ):.2f}: ||tau k_lam||/||k||={resid:.4e}  <k|tau|k>={ray:+.4e}  |  random-even ||tau v||: {['%.3f'%r for r in rands]}",flush=True)
print("\nIf ||tau k_lam|| << random ||tau v||: k_lam IS a near-eigenvector (tau k=0), sharp RH criterion substantiated.",flush=True)
