"""S2876: sign of I_k = int_0^inf w P^2/(Q^2 (t^2-g_k^2)) dt across (lam2, N); also count j<=N with gamma_j<=nu_j."""
import sys, mpmath as mp
mp.mp.dps=30
Z=[mp.zetazero(k).imag for k in range(1,61)]
for lam2,Ns in [(14,[8,16,24,40]),(4,[6,10,14,20,30]),(2,[4,8,12])]:
    L=mp.log(lam2)
    for N in Ns:
        nu=[2*mp.pi*n/L for n in range(1,N+1)]; G=Z[:N]
        viol=sum(1 for g,v in zip(G,nu) if g<=v)
        R2=lambda t: mp.fprod((t*t-g*g)/(t*t-v*v) for g,v in zip(G,nu))**2
        w=lambda t: 4*mp.sin(t*L/2)**2/L
        signs=[]
        for k in range(1,min(N,5)+1):
            gk=G[k-1]
            f=lambda t: (w(t)/t**2 if t!=0 else L)*R2(t)/(t*t-gk*gk)
            pts=sorted(set([mp.mpf(0)]+G+nu+[G[-1]*2,G[-1]*8]))+[mp.inf]
            I=mp.quad(f,pts,maxdegree=8)
            signs.append(mp.nstr(I,3))
        print(f"lam2={lam2} N={N}: #(gamma_j<=nu_j)={viol}  I_k (k=1..5) = {signs}",flush=True)
