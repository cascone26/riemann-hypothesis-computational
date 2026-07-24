# Is the §56 criticality FORCED? Test: does each prime-p form near-ANNIHILATE the even ground state v_a?
# Build even ground state v_a of tau; for each active prime p compute <v_a|P_p|v_a> (p-th prime Rayleigh)
# and d(min-eig)/dw_p = -<v_a|P_p|v_a>. Compare to random even vectors. Small => criticality explained.
import mpmath as mp
mp.mp.dps=28
def run(LSQ,N=12):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
    Rq=TWO_PI*N/L+55; G=900
    grid=[(-Rq+2*Rq*mp.mpf(i)/G) for i in range(G+1)]; hh=2*Rq/G
    wq=[hh]*(G+1); wq[0]=hh/2; wq[-1]=hh/2
    Om=[mp.re(mp.digamma(mp.mpf('0.25')+0.5j*r))-LOGPI for r in grid]
    def ehat(r,n):
        s=TWO_PI*n/L; d=s-r
        if abs(d)<mp.mpf('1e-22'): return mp.sqrt(L)
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
        return (mp.log(p),p) if(ok and p>0) else (mp.mpf(0),0)
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    # full prime term and per-base-prime pieces
    bases=sorted(set(vm(k)[1] for k in range(2,int(mp.floor(LSQ))+1) if vm(k)[0]>0))
    def Pmat(restrict_base=None):
        Pm=mp.matrix(D,D)
        for k in range(2,int(mp.floor(LSQ))+1):
            lam,base=vm(k)
            if lam>0 and (restrict_base is None or base==restrict_base):
                for i,n in enumerate(idx):
                    for j,m in enumerate(idx):
                        Pm[i,j]+=lam*mp.mpf(k)**mp.mpf('-0.5')*qk(n,m,mp.log(k))
        return Pm
    A=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j,m in enumerate(idx): A[i,j]=arch(i,j)+Bnm(n,m)
    Pfull=Pmat(); T=A-Pfull
    # even ground state
    def col_e(n):
        v=mp.matrix(D,1)
        if n==0: v[idx.index(0)]=1
        else: v[idx.index(n)]=1/mp.sqrt(2); v[idx.index(-n)]=1/mp.sqrt(2)
        return v
    Ec=[col_e(0)]+[col_e(n) for n in range(1,N+1)]; kk=len(Ec)
    BE=mp.matrix(kk,kk)
    for a in range(kk):
        for b in range(kk): BE[a,b]=mp.fsum(Ec[a][i]*T[i,j]*Ec[b][j] for i in range(D) for j in range(D))
    ev,Vec=mp.eigsy(BE); o=sorted(range(kk),key=lambda t:ev[t])
    # ground state in full-D coords
    vg=mp.matrix(D,1)
    for a in range(kk):
        ca=Vec[a,o[0]]
        for i in range(D): vg[i]+=ca*Ec[a][i]
    ng=mp.sqrt(mp.fsum(vg[i]**2 for i in range(D)))
    for i in range(D): vg[i]/=ng
    def ray(M,v): return float(mp.fsum(v[i]*(M*v)[i] for i in range(D)))
    # random even vector for comparison
    import hashlib
    vr=mp.matrix(D,1)
    for i,n in enumerate(idx):
        h=int(hashlib.md5(f"{abs(n)}".encode()).hexdigest(),16)%1000; vr[i]=mp.mpf(h)/1000-mp.mpf('0.5')
    nr=mp.sqrt(mp.fsum(vr[i]**2 for i in range(D)))
    for i in range(D): vr[i]/=nr
    print(f"lam^2={float(LSQ):.1f} primes={bases}: min-eig(even)={float(ev[o[0]]):+.3e}",flush=True)
    print(f"  <v_a|P_full|v_a>={ray(Pfull,vg):+.4f}  | per prime p: <v_a|P_p|v_a> (=-d(mineig)/dw_p):",flush=True)
    for p in bases:
        Pp=Pmat(p); print(f"     p={p}: <v_a|P_p|v_a>={ray(Pp,vg):+.5f}   (random-even: {ray(Pp,vr):+.4f})",flush=True)
run(mp.mpf('4.0'))
run(mp.mpf('6.0'))
