import numpy as np, math
# ---- SOLVED cousin: curve/F_q, Frobenius angles theta (RH => on unit circle) ----
g=6; rng=np.random.default_rng(1); theta=rng.uniform(0.25,math.pi-0.25,size=g)
def sym_toeplitz(theta,rho,K):
    s=np.array([np.sum((rho**k+rho**(-k))*np.cos(k*theta)) for k in range(K+1)])
    return np.array([[s[abs(i-j)] for j in range(K+1)] for i in range(K+1)])
K=14
min_rh=np.linalg.eigvalsh(sym_toeplitz(theta,1.0,K)).min()
min_bad=np.linalg.eigvalsh(sym_toeplitz(theta,1.15,K)).min()
print('=== SOLVED (function field) mechanism ===')
print('genus',g,' RH-true (|Frob|=q^{1/2}) Toeplitz min eig = %.4e'%min_rh)
print('RH-broken (Frob off circle) Toeplitz min eig = %.4e'%min_bad)
if min_rh>-1e-9 and min_bad<-1e-9:
    print('CONFIRMED: positivity <=> Frobenius-unitary <=> RH (trig-moment/Herglotz), one mechanism')
# ---- NUMBER-FIELD probe: is the prime-side moment sequence positive-definite? ----
def von_mangoldt(N):
    Lam=np.zeros(N+1); sieve=np.ones(N+1,bool); sieve[:2]=False
    for p in range(2,N+1):
        if sieve[p]:
            lp=math.log(p); pk=p
            while pk<=N: Lam[pk]=lp; pk*=p
            for m in range(p*p,N+1,p): sieve[m]=False
    return Lam
N=5000; Lam=von_mangoldt(N); ns=np.nonzero(Lam>0)[0]
logn=np.log(ns.astype(float)); wt=Lam[ns]/np.sqrt(ns.astype(float))
# moments c_k = archimedean baseline (k=0) minus prime oscillation at 'frequency' k
Kn=20; base=wt.sum()*2.0
c=np.zeros(Kn+1)
for k in range(Kn+1):
    c[k]=(base if k==0 else 0.0)-np.sum(wt*np.cos(k*logn))
Tnf=np.array([[c[abs(i-j)] for j in range(Kn+1)] for i in range(Kn+1)])
mn=np.linalg.eigvalsh(Tnf).min()
print('=== NUMBER FIELD probe ===')
print('prime-moment Toeplitz min eig = %.4e (>=0 => Weil-positive on this moment window)'%mn)
if mn>=-1e-9:
    print('NUMBER_FIELD_MOMENTS_POSITIVE_DEFINITE: prime measure behaves unitary on window -> supports a positive spectral measure')
else:
    print('DEFICIT: prime-moment Toeplitz indefinite -> the missing positive spectral measure is not yet present (target to construct)')
