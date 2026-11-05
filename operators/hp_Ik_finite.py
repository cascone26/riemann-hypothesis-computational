"""S2877: exact finite-sum form of I_k via Plancherel: I_k = pi * sum_n P_k(nu_n)^2 (nu_n^2-g_k^2)/Q'(nu_n)^2."""
import mpmath as mp
mp.mp.dps=40
Z=[mp.zetazero(k).imag for k in range(1,61)]
for lam2,N in [(14,8),(14,16),(4,20),(4,30),(2,8)]:
    L=mp.log(lam2); h=2*mp.pi/L; nus=[h*n for n in range(-N,N+1)]; G=Z[:N]
    Qp=lambda x: mp.fprod(x-v for v in nus if v!=x)
    out=[]
    for k in range(1,4):
        gk=G[k-1]
        Pk=lambda t: mp.fprod(t*t-g*g for j,g in enumerate(G) if j!=k-1)
        terms=[Pk(v)**2*(v*v-gk*gk)/Qp(v)**2 for v in nus]
        I=mp.pi*mp.fsum(terms); neg=mp.fsum(x for x in terms if x<0); pos=mp.fsum(x for x in terms if x>0)
        out.append(f"I_{k}={mp.nstr(I,4)} (pos/|neg|={mp.nstr(pos/abs(neg),3)})")
    print(f"lam2={lam2} N={N}: "+"; ".join(out))
