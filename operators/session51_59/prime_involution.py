# Is P2 = w2 * R2 with R2 an EXACT reflection (R2^2 = I)?  Test P2^2 vs w2^2 * I.
# And: is the sign of <u0|R2|u0> (u0=arch-min-vec) FORCED by parity (gamma: n<->-n reflection)?
import mpmath as mp
mp.mp.dps=35
def pieces(LSQ,N=10):
    L=2*mp.log(mp.sqrt(LSQ)); TWO_PI=2*mp.pi; LOGPI=mp.log(mp.pi)
    idx=list(range(-N,N+1)); D=len(idx)
    lg2=mp.log(2); w2=lg2/mp.sqrt(2)
    def qk(n,m,y):
        if n==m: return 2*(1-abs(y)/L)*mp.cos(TWO_PI*n*y/L)
        return (mp.sin(TWO_PI*m*y/L)-mp.sin(TWO_PI*n*y/L))/(mp.pi*(n-m))
    P=mp.matrix(D,D)
    for i,n in enumerate(idx):
        for j,m in enumerate(idx): P[i,j]=w2*qk(n,m,lg2)
    return P,w2,D,idx,L
for LSQ in [mp.mpf('2.5'),mp.mpf('2.72')]:
    for N in [10,16]:
        P,w2,D,idx,L=pieces(LSQ,N)
        P2sq=P*P
        # ||P^2 - w2^2 I|| / ||w2^2 I||
        diff=mp.mpf(0); nrm=mp.mpf(0)
        for i in range(D):
            for j in range(D):
                tgt=w2*w2 if i==j else mp.mpf(0)
                diff+=abs(P2sq[i,j]-tgt)**2; nrm+=abs(tgt)**2
        rel=mp.sqrt(diff/nrm)
        print(f"lam^2={float(LSQ):.2f} N={N}: ||P2^2 - w2^2 I||/||w2^2 I|| = {mp.nstr(rel,4)}  (0 => exact reflection; w2={mp.nstr(w2,5)})",flush=True)
