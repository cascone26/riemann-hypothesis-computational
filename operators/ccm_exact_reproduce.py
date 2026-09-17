#!/usr/bin/env python3
"""
CCM exact spectral-triple operator reproduction (arXiv:2511.22755), cycle 317.
QW_lambda(V_n,V_m) = W_R + W_p + W_02 in the periodic scaling basis U_n=L^{-1/2}e^{2pi i n x/L}.
  omega(x)=q(U_n,U_m)(x):  n!=m -> [sin(2pi m x/L)-sin(2pi n x/L)]/(pi(n-m)); n=m -> 2(L-x)cos(2pi n x/L)/L.
  W_R(n,m)= (omega(0)/2)(gamma+log(4pi (e^L-1)/(e^L+1))) + INT_0^L [e^{x/2}omega(x)-omega(0)]/(e^x-e^{-x})dx
  W_p(n,m)= sum_{1<k<=e^L} Lambda(k) k^{-1/2} omega(log k)
  W_02(n,m)= 32 L sinh^2(L/4)(L^2-16pi^2 m n)/((L^2+16pi^2 m^2)(L^2+16pi^2 n^2))
Then: smallest eigenvalue eps_N (simple,even), eigenvector xi; xihat(z)=2 L^{-1/2} sin(zL/2) sum_j xi_j/(z-2pi j/L);
zeros of xihat = spectrum ~ zeta zeros. Compare first zero to gamma_1=14.1347251417.
"""
import mpmath as mp
mp.mp.dps = 40
PI = mp.pi
GAMMA = mp.euler

def build(lam, N, Kprime=None):
    L = 2*mp.log(lam)
    def rho(x): return mp.e**(x/2)/(mp.e**x - mp.e**(-x))
    # single-frequency arch moments (n>=0); J_sin(n)=INT rho*sin(2pi n x/L) (sin kills 0-singularity)
    def Jsin(n):
        if n==0: return mp.mpf(0)
        return mp.quad(lambda x: rho(x)*mp.sin(2*PI*n*x/L), [0, L/8, L/2, L])
    # diagonal regularized arch integral D(n)=INT [e^{x/2}*2(L-x)cos(2pi n x/L)/L - 2]/(e^x-e^{-x}) dx
    def Ddiag(n):
        def integ(x):
            num = mp.e**(x/2)*2*(L-x)*mp.cos(2*PI*n*x/L)/L - 2
            den = mp.e**x - mp.e**(-x)
            if x < mp.mpf(10)**(-18):   # regularized limit at 0
                # num ~ x*(...); den ~ 2x -> ratio finite; use series-safe small value
                return integ(mp.mpf(10)**(-15))
            return num/den
        return mp.quad(integ, [mp.mpf(10)**(-15), L/8, L/2, L])
    idx = list(range(-N, N+1))
    Js = {}
    for n in range(0, N+1): Js[n] = Jsin(n)
    def Jsig(n): return Js[abs(n)]*(1 if n>=0 else -1)
    Dd = {n: Ddiag(n) for n in range(-N, N+1)}
    archconst = GAMMA + mp.log(4*PI*(mp.e**L-1)/(mp.e**L+1))
    # von Mangoldt up to e^L
    Nmax = int(mp.floor(mp.e**L))
    if Kprime: Nmax = min(Nmax, Kprime)
    Lam = [mp.mpf(0)]*(Nmax+1); comp=[False]*(Nmax+1)
    for p in range(2, Nmax+1):
        if not comp[p]:
            for mmm in range(2*p, Nmax+1, p): comp[mmm]=True
            lp=mp.log(p); pk=p
            while pk<=Nmax: Lam[pk]=lp; pk*=p
    primes = [(k, Lam[k]) for k in range(2, Nmax+1) if Lam[k]>0]
    def omega(n, m, y):   # q(U_n,U_m)(y), y in [0,L]
        if n==m:
            return 2*(L-y)*mp.cos(2*PI*n*y/L)/L
        return (mp.sin(2*PI*m*y/L)-mp.sin(2*PI*n*y/L))/(PI*(n-m))
    def W02(n,m):
        return 32*L*mp.sinh(L/4)**2*(L**2-16*PI**2*m*n)/((L**2+16*PI**2*m**2)*(L**2+16*PI**2*n**2))
    n2=2*N+1
    Q = mp.matrix(n2, n2)
    for ii,n in enumerate(idx):
        for jj in range(ii, n2):
            m = idx[jj]
            # arch
            if n==m:
                WR = archconst + Dd[n]
            else:
                WR = (Jsig(m)-Jsig(n))/(PI*(n-m))
            # prime
            Wp = mp.mpf(0)
            for (k,lk) in primes:
                Wp += lk*mp.power(k, mp.mpf('-0.5'))*omega(n,m,mp.log(k))
            v = WR + Wp + W02(n,m)
            Q[ii,jj]=v; Q[jj,ii]=v
    return Q, L, idx

def first_zero(lam, N, Kprime=None):
    Q, L, idx = build(lam, N, Kprime)
    E, V = mp.eigsy(Q)
    # smallest eigenvalue index
    mn=min(range(len(idx)), key=lambda i: E[i])
    xi = [V[r, mn] for r in range(len(idx))]
    eps = E[mn]
    def xihat(z):
        s = mp.mpf(0)
        for r,j in enumerate(idx):
            s += xi[r]/(z - 2*PI*j/L)
        return 2*L**mp.mpf('-0.5')*mp.sin(z*L/2)*s
    # find zero near 14.1347 by sign change + refine
    z0 = mp.findroot(xihat, mp.mpf('14.13'))
    return eps, z0

if __name__=="__main__":
    G1 = mp.mpf('14.134725141734693')
    print("CCM exact operator: first xihat zero vs gamma_1=14.134725141734693")
    for lam2, N in [(12,30),(13,30),(13,50)]:
        lam = mp.sqrt(lam2)
        eps, z0 = first_zero(lam, N)
        print(f"  lambda=sqrt({lam2}), N={N}: eps_N={mp.nstr(eps,5)}  first zero={mp.nstr(z0,18)}  |err|={mp.nstr(abs(z0-G1),4)}", flush=True)
