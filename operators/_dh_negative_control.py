import numpy as np, math
from scipy.special import digamma
# NEGATIVE CONTROL: Davenport-Heilbronn D(s) violates RH (has off-line zeros), NO Euler product.
# If mu_lambda detects RH, mu(D) should be ROBUSTLY NEGATIVE vs zeta's razor-zero.
# D(s)=Σ a_n n^-s, a=[1,κ,-κ,-1,0] period 5 (χ mod 5, χ(2)=i). Log-deriv coeffs Λ_D via recursion.

kappa = (math.sqrt(10-2*math.sqrt(5)) - 2)/(math.sqrt(5)-1)
def a_coef(n):  # Dirichlet coeff of D
    r = n % 5
    return {1:1.0, 2:kappa, 3:-kappa, 4:-1.0, 0:0.0}[r]

def lambda_D(N):
    # -D'/D = Σ c_n n^-s ; recursion: -a_n log n = Σ_{d|n} c_d a_{n/d}; Λ_D(n) = -c_n
    a = np.array([0.0]+[a_coef(n) for n in range(1,N+1)])  # 1-indexed
    c = np.zeros(N+1)
    for n in range(1,N+1):
        s = 0.0
        for d in range(1,n):
            if n % d == 0:
                s += c[d]*a[n//d]
        # a[1]=1 so c[n] = (-a_n log n - s)/a[1]
        c[n] = (-a[n]*math.log(n) - s) if n>1 else 0.0
    return -c   # Λ_D(n)

def omega_D(r):  # D archimedean: Γ((s+1)/2), conductor 5 -> Re ψ(3/4+ir/2) - log(π/5)
    return np.real(digamma(0.75 + 0.5j*r)) - np.log(np.pi/5)
def omega_Z(r):  # zeta: Re ψ(1/4+ir/2) - log π
    return np.real(digamma(0.25 + 0.5j*r)) - np.log(np.pi)

def build_mu(omega_fn, coef_arr, Ncut, N=40, fourW=6.0, a=14, b=44):
    L = fourW/4
    RG = np.arange(-(b+70), (b+70), 0.004); OM = omega_fn(RG)
    CEN = np.linspace(a, b, N)
    def fej(x):
        y=L*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN])
    arch = (P*OM[None,:])@P.T*0.004/(2*np.pi)
    ln = np.array([math.log(n) for n in range(2,Ncut+1)])
    lw = np.array([coef_arr[n]/math.sqrt(n) for n in range(2,Ncut+1)])
    cu = np.cos(np.outer(ln, RG)); prime=np.zeros((N,N))
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*0.004/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+prime.T)
    tau = arch - prime
    ev = np.linalg.eigvalsh(tau)
    return ev[0], np.max(np.abs(ev))

# zeta coeffs = von Mangoldt
def vm(N):
    lam=np.zeros(N+1)
    for p in range(2,N+1):
        pk=p
        while pk<=N:
            lam[pk]=math.log(p); pk*=p
        # only prime powers; sieve properly
    lam=np.zeros(N+1)
    for n in range(2,N+1):
        m,p=n,0; f=2; ok=True
        while f*f<=m:
            if m%f==0:
                if p==0:p=f
                elif p!=f:ok=False;break
                while m%f==0:m//=f
            f+=1
        if m>1:
            if p==0:p=m
            elif p!=m:ok=False
        if ok and p>0: lam[n]=math.log(p)
    return lam

Ncut=120
LD = lambda_D(Ncut)
VM = vm(Ncut)
print(f"kappa = {kappa:.6f}")
print(f"Λ_D(n) n=2..10: {[round(LD[n],3) for n in range(2,11)]}  (NOT prime-supported: composites+negatives)")
muZ, nZ = build_mu(omega_Z, VM, Ncut)
muD, nD = build_mu(omega_D, LD, Ncut)
print(f"\nZETA (RH true):              mu_min = {muZ:+.4e}   |mu|/||tau|| = {abs(muZ)/nZ:.3e}")
print(f"DAVENPORT-HEILBRONN (RH false): mu_min = {muD:+.4e}   |mu|/||tau|| = {abs(muD)/nD:.3e}")
print(f"\nratio muD/muZ = {muD/muZ:.1f}x   {'✓ DH robustly more negative = mu DETECTS RH-violation' if muD < muZ*5 else 'inconclusive'}")
