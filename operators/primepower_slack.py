import numpy as np
from scipy.special import digamma
from scipy.linalg import eigh
def vonmangoldt(N):
    lam=np.zeros(N+1);isc=np.zeros(N+1,dtype=bool)
    for p in range(2,N+1):
        if not isc[p]:
            isc[2*p::p]=True;lp=np.log(p);pk=p
            while pk<=N: lam[pk]=lp;pk*=p
    return lam
def is_prime_power_ge2(n):  # True if n=p^k, k>=2
    m=n;p=2;f=0;base=0
    while p*p<=m:
        if m%p==0:
            base=p
            while m%p==0: m//=p;f+=1
            return (m==1 and f>=2)
        p+=1
    return False
NMAX=3000;LAM=vonmangoldt(NMAX);NS=np.arange(2,NMAX+1);LOGN=np.log(NS);LAMW0=LAM[2:]/np.sqrt(NS)
a,b,W=14.0,44.0,1.0;fourW=4*W;dr=0.004
RGRID=np.arange(-(b+120),(b+120),dr);OMEGA=np.real(digamma(0.25+0.5j*RGRID))-np.log(np.pi)
def fej(x):
    y=W*x;out=np.ones_like(y);nz=np.abs(y)>1e-12;out[nz]=(np.sin(y[nz])/y[nz])**2;return out
def fc(c):
    y=W*(0.5j-c);y2=W*(0.5j+c);return (np.sin(y)/y)**2+(np.sin(y2)/y2)**2
mask=LOGN<fourW;ln=LOGN[mask];nn=NS[mask];w0=LAMW0[mask];cosu=np.cos(np.outer(ln,RGRID))
print("Does dropping prime-POWERS (k>=2) flip arch>=prime to arch<prime? (=> powers provide the RH slack)")
for N in [13,20]:
    CEN=np.linspace(a,b,N)
    P=np.array([fej(RGRID-c)+fej(RGRID+c) for c in CEN]);Gm=(P@P.T)*dr
    arch=(P*OMEGA[None,:])@P.T*dr/(2*np.pi)
    fp=np.array([fc(c) for c in CEN]);pole=2.0*np.real(np.outer(fp,fp))
    Mn=[]
    for jj in range(len(ln)):
        col=cosu[jj];Mk=np.zeros((N,N))
        for j in range(N): Mk[j,:]=2.0*(((P[j][None,:]*P)@col)*dr/(2*np.pi))
        Mn.append(0.5*(Mk+Mk.T))
    prime=np.tensordot(w0,np.array(Mn),axes=1)
    try: w,V=eigh(pole+arch-prime,Gm)
    except Exception: print(f"  N={N}: singular");continue
    v=V[:,0]
    av=float(v@arch@v);pv_full=float(v@prime@v);polv=float(v@pole@v)
    pv_primes=sum(w0[jj]*float(v@Mn[jj]@v) for jj in range(len(ln)) if not is_prime_power_ge2(int(nn[jj])))
    pv_powers=sum(w0[jj]*float(v@Mn[jj]@v) for jj in range(len(ln)) if is_prime_power_ge2(int(nn[jj])))
    print(f"  N={N}: arch={av:.5f}  pole={polv:+.2e}  prime_full={pv_full:.5f}  prime_primesONLY={pv_primes:.5f}  powers(k>=2)={pv_powers:+.5f}")
    print(f"        arch+pole - prime_full   = {av+polv-pv_full:+.3e}  (true margin, >0 = RH-ok)")
    print(f"        arch+pole - prime_primesONLY = {av+polv-pv_primes:+.3e}  (if <0: prime powers ARE the slack keeping it positive)")
