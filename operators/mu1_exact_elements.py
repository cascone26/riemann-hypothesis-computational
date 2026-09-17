import mpmath as mp
mp.mp.dps=30
def ker(th,ph):
    t=abs(mp.cos(th)-mp.cos(ph))
    return t*mp.log(t) if t>0 else mp.mpf(0)
def Kel(m,n):
    def inner(th):
        f=lambda ph: ker(th,ph)*mp.cos(n*ph)
        return mp.quad(f,[0,th,mp.pi])
    g=lambda th: inner(th)*mp.cos(m*th)
    return mp.quad(g,[0,mp.pi])
def Mel(m,n):  # closed-form-ish mass: (1/mn) INT_0^pi sin(m p)sin(n p)sin p dp
    f=lambda p: mp.sin(m*p)*mp.sin(n*p)*mp.sin(p)
    return mp.quad(f,[0,mp.pi])/(m*n)
def lam_min(K):
    Kmat=mp.matrix(K,K); M=mp.matrix(K,K)
    for i in range(1,K+1):
        for j in range(i,K+1):
            v=Kel(i,j); Kmat[i-1,j-1]=v; Kmat[j-1,i-1]=v
            w=Mel(i,j); M[i-1,j-1]=w; M[j-1,i-1]=w
    # generalized eig via Cholesky of M: M=L L^T, solve L^-1 K L^-T
    L=mp.cholesky(M); Li=L**-1
    A=Li*Kmat*Li.T; A=(A+A.T)/2
    ev=mp.eigsy(A,eigvals_only=True)
    return min(ev)
print("EXACT-element spectral mu1 (mpmath dps=30):",flush=True)
prev=None
for K in [6,8,10,12,14]:
    lm=lam_min(K); mu1=1+lm/2
    print(f"  K={K:3d}: lambda={mp.nstr(lm,14)}  mu1={mp.nstr(mu1,14)}",flush=True)
    prev=(lm,mu1)
lm,mu1=prev
for basis in [['euler','log(2)','pi'],['euler','log(2)','log(pi)','pi','catalan'],['euler','log(2)','pi','zeta(3)','catalan','pi^2']]:
    print("  identify(mu1):",mp.identify(mu1,basis),flush=True)
