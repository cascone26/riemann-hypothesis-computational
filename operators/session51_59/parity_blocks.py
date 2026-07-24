# PARITY-BLOCK decomposition of one-prime positivity. gamma: n<->-n (functional-eq reflection).
# tau=(A+B)-P2 is gamma-symmetric => splits into EVEN + ODD blocks. Positivity <=> both blocks >=0.
# Compute min eig of tau, A+B in EACH parity block, across the one-prime regime. Is the arch-neg
# direction (and its prime rescue) confined to ONE parity sector? That would make the sign FORCED.
import mpmath as mp
mp.mp.dps=35
def build_blocks(LSQ,N=10):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
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
    lg2=mp.log(2); w2=lg2/mp.sqrt(2)
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    AB=mp.matrix(D,D); P=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j,m in enumerate(idx):
            AB[i,j]=arch(i,j)+Bnm(n,m); P[i,j]=w2*qk(n,m,lg2)
    # parity basis: even = e0, (e_n+e_-n)/sqrt2 (n>0); odd = (e_n-e_-n)/sqrt2 (n>0)
    pos=[k for k in range(D) if idx[k]>0]; zero=[k for k in range(D) if idx[k]==0]
    def col_e(n):  # even combo
        v=mp.matrix(D,1); 
        if n==0: v[idx.index(0)]=1
        else: v[idx.index(n)]=1/mp.sqrt(2); v[idx.index(-n)]=1/mp.sqrt(2)
        return v
    def col_o(n):
        v=mp.matrix(D,1); v[idx.index(n)]=1/mp.sqrt(2); v[idx.index(-n)]=-1/mp.sqrt(2); return v
    Ecols=[col_e(0)]+[col_e(n) for n in range(1,N+1)]
    Ocols=[col_o(n) for n in range(1,N+1)]
    def blk(Mx,cols):
        k=len(cols); B=mp.matrix(k,k)
        for a in range(k):
            for b in range(k):
                B[a,b]=mp.fsum(cols[a][i]*Mx[i,j]*cols[b][j] for i in range(D) for j in range(D))
        return B
    return AB,P,Ecols,Ocols,blk
print(f"{'lam^2':>6} | {'min tau EVEN':>13} {'min tau ODD':>12} | {'min AB EVEN':>12} {'min AB ODD':>11}",flush=True)
for LSQ in [mp.mpf('2.2'),mp.mpf('2.5'),mp.mpf('2.72')]:
    AB,P,Ec,Oc,blk=build_blocks(LSQ)
    T=AB-P
    tE=min(mp.eigsy(blk(T,Ec),eigvals_only=True)); tO=min(mp.eigsy(blk(T,Oc),eigvals_only=True))
    aE=min(mp.eigsy(blk(AB,Ec),eigvals_only=True)); aO=min(mp.eigsy(blk(AB,Oc),eigvals_only=True))
    print(f"{float(LSQ):6.2f} | {mp.nstr(tE,5):>13} {mp.nstr(tO,5):>12} | {mp.nstr(aE,5):>12} {mp.nstr(aO,5):>11}",flush=True)
print("\nWhich parity sector holds the arch-negativity? Does the prime rescue live in the SAME sector?",flush=True)
