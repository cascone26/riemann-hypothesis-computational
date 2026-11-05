"""CCM Weil matrix split into W02 / W_R / W_p parts (same formulas as ccm_exact_reproduce.build)."""
import mpmath as mp, numpy as np
PI = mp.pi; GAMMA = mp.euler
def parts(lam2, N, dps=30, as_mp=False):
    mp.mp.dps = dps
    L = mp.log(lam2)
    rho = lambda x: mp.e**(x/2)/(mp.e**x - mp.e**(-x))
    def Jsin(n):
        if n == 0: return mp.mpf(0)
        return mp.quad(lambda x: rho(x)*mp.sin(2*PI*n*x/L), [0, L/8, L/2, L])
    def Ddiag(n):
        f = lambda x: (mp.e**(x/2)*2*(L-x)*mp.cos(2*PI*n*x/L)/L - 2)/(mp.e**x - mp.e**(-x))
        return mp.quad(f, [0, L/8, L/2, L])
    idx = list(range(-N, N+1)); M = len(idx)
    Js = {n: Jsin(n) for n in range(N+1)}
    Dd = {n: Ddiag(n) for n in range(-N, N+1)}
    Jsig = lambda n: Js[abs(n)]*(1 if n >= 0 else -1)
    archconst = GAMMA + mp.log(4*PI*(mp.e**L-1)/(mp.e**L+1))
    Nmax = int(mp.floor(mp.e**L + mp.mpf('1e-20')))
    pw = []
    for k in range(2, Nmax+1):
        # von Mangoldt
        p = None
        for q in range(2, k+1):
            if k % q == 0: p = q; break
        kk = k
        while kk % p == 0: kk //= p
        if kk == 1 and k < mp.e**L: pw.append((k, mp.log(p)))
    def omega(n, m, y):
        if n == m: return 2*(L-y)*mp.cos(2*PI*n*y/L)/L
        return (mp.sin(2*PI*m*y/L)-mp.sin(2*PI*n*y/L))/(PI*(n-m))
    W02 = mp.matrix(M, M); WR = mp.matrix(M, M); WP = mp.matrix(M, M)
    for i, n in enumerate(idx):
        for j, m in enumerate(idx):
            W02[i, j] = (32*L*mp.sinh(L/4)**2*(L**2-16*PI**2*m*n)/((L**2+16*PI**2*m**2)*(L**2+16*PI**2*n**2)))
            WR[i, j] = (archconst + Dd[n]) if n == m else ((Jsig(m)-Jsig(n))/(PI*(n-m)))
            WP[i, j] = (sum(lk*mp.power(k, -0.5)*omega(n, m, mp.log(k)) for k, lk in pw))
    if as_mp: return W02, WR, WP, L, pw
    f = lambda A: np.array([[float(A[i, j]) for j in range(M)] for i in range(M)])
    return f(W02), f(WR), f(WP), float(L), pw
