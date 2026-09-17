import mpmath as mp
mp.mp.dps = 35
# High-precision Weil-form min eigenvalue mu_lambda, to see BELOW the float64 ~1e-12 floor.
# tau = arch + boundary - prime, small N, mpmath eigsy. Measure how mu_lambda decays with lambda.

def vm(kmax):
    tab=[]
    for k in range(2,kmax+1):
        n,p,f,ok=k,0,2,True
        while f*f<=n:
            if n%f==0:
                if p==0:p=f
                elif p!=f:ok=False;break
                while n%f==0:n//=f
            f+=1
        if n>1:
            if p==0:p=n
            elif p!=n:ok=False
        if ok and p>0: tab.append((k, mp.log(p)))
    return tab

def omega(r):  # Re psi(1/4 + i r/2) - log pi
    return mp.re(mp.digamma(mp.mpf(1)/4 + mp.mpf(1)/2*1j*r)) - mp.log(mp.pi)

def ehat(r,n,L):
    s=2*mp.pi*n/L; d=s-r
    if abs(d)<mp.mpf(10)**-20: return mp.sqrt(L)
    return (1/mp.sqrt(L))*mp.e**(mp.mpf(1)/2*1j*d*L)*2*mp.sin(mp.mpf(1)/2*d*L)/d

def qk(n,m,y,L):
    if n==m: return 2*(1-abs(y)/L)*mp.cos(2*mp.pi*n*y/L)
    return (mp.sin(2*mp.pi*m*y/L)-mp.sin(2*mp.pi*n*y/L))/(mp.pi*(n-m))

def w02(n,m,L):
    sh=mp.sinh(L/4)**2
    num=(L**2-16*mp.pi**2*m*n); den=(L**2+16*mp.pi**2*m**2)*(L**2+16*mp.pi**2*n**2)
    return 32*L*sh*num/den

def mu_lambda(LSQ, N, G=600):
    L=2*mp.log(mp.sqrt(LSQ)); VM=vm(LSQ)
    idx=list(range(-N,N+1)); D=len(idx)
    R=2*mp.pi*N/L+40
    grid=[(-R + 2*R*i/(G)) for i in range(G+1)]; h=grid[1]-grid[0]
    w=[h]*(G+1); w[0]=w[-1]=h/2
    Om=[omega(r) for r in grid]
    # arch[a,b] = (1/2pi) sum_g w_g Om_g ehat_a(g) conj(ehat_b(g))
    E=[[ehat(grid[g],n,L) for g in range(G+1)] for n in idx]
    tau=mp.zeros(D,D)
    for a in range(D):
        for b in range(a,D):
            arch=mp.mpf(0)
            for g in range(G+1):
                arch+= w[g]*Om[g]*E[a][g]*mp.conj(E[b][g])
            arch=mp.re(arch)/(2*mp.pi)
            pr=mp.mpf(0)
            for k,Lam in VM:
                pr+= Lam*k**mp.mpf('-0.5')*qk(idx[a],idx[b],mp.log(k),L)
            val=arch + w02(idx[a],idx[b],L) - pr
            tau[a,b]=val; tau[b,a]=val
    ev=mp.eigsy(tau, eigvals_only=True)
    return min(ev)

print("HIGH-PRECISION (dps=35) mu_lambda below float64 floor:")
print(f"{'LSQ':>5} {'lam':>6} {'N':>3} {'mu_min (mpmath)':>28}")
for LSQ,N in [(14,7),(30,8),(60,9),(120,10)]:
    mu=mu_lambda(LSQ,N)
    print(f"{LSQ:>5} {float(mp.sqrt(LSQ)):>6.2f} {N:>3} {mp.nstr(mu,12):>28}", flush=True)
