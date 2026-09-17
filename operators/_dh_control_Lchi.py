import numpy as np, math
from scipy.special import digamma
import operators._dh_negative_control as DH

# CONTROL: single Dirichlet L(s,χ) mod5 (entire, RH-TRUE, SAME arch as DH). Complex Hermitian Weil form.
# coeff term = χ(n)Λ(n) (prime-supported, complex). If κ(Lχ)=0 where κ(DH)>0 => genuine off-line detection.
chi = {1:1+0j, 2:1j, 3:-1j, 4:-1+0j, 0:0j}
def vmL(N):  # χ(n)Λ(n): von Mangoldt on prime powers, times χ
    out=np.zeros(N+1, complex)
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
        if ok and p>0: out[n]=math.log(p)*chi[n%5]
    return out

def build_tau_cplx(coef_arr, Ncut, N, fourW, a, b, real=False):
    L = fourW/4
    RG = np.arange(-(b+80), (b+80), 0.004); OM = DH.omega_D(RG)
    CEN = np.linspace(a, b, N)
    def fej(x):
        y=L*x; o=np.ones_like(y); nz=np.abs(y)>1e-12; o[nz]=(np.sin(y[nz])/y[nz])**2; return o
    P = np.array([fej(RG-c)+fej(RG+c) for c in CEN])
    arch = (P*OM[None,:])@P.T*0.004/(2*np.pi)
    ln = np.array([math.log(n) for n in range(2,Ncut+1)])
    lw = np.array([coef_arr[n]/math.sqrt(n) for n in range(2,Ncut+1)])
    cu = np.cos(np.outer(ln, RG))
    prime = np.zeros((N,N), complex)
    for j in range(N):
        gp=((P[j][None,:]*P)@cu.T)*0.004/(2*np.pi); prime[j,:]=2.0*(gp*lw[None,:]).sum(axis=1)
    prime=0.5*(prime+np.conj(prime.T))
    tau = arch.astype(complex) - prime
    tau=0.5*(tau+np.conj(tau.T))
    ev=np.linalg.eigvalsh(tau)
    return int((ev<-1e-6).sum()), ev[0]

Ncut=3000
LC = vmL(Ncut); LD = DH.lambda_D(Ncut)
print("CONTROL: L(s,χ) mod5 (RH-TRUE entire) vs DH (RH-FALSE), SAME arch, kappa=#(eig<-1e-6)")
print(f"{'window':>10} {'N':>4} {'kappa(Lχ,RH-true)':>18} {'kappa(DH,RH-false)':>18}")
for (a,b,fW,N) in [(14,44,5.0,60),(78,92,11.0,110)]:
    kL,muL = build_tau_cplx(LC, Ncut, N, fW, a, b)
    kD,muD = build_tau_cplx(LD, Ncut, N, fW, a, b)
    print(f"[{a:>3},{b:>3}] {N:>4} {kL:>18} {kD:>18}", flush=True)
print("\nDH resolution-stability at [78,92] (does kappa converge or grow=artifact?):")
for N in [90,110,140]:
    kD,muD = build_tau_cplx(LD, Ncut, N, 11.0, 78, 92)
    print(f"  N={N}: kappa(DH)={kD}  mu={muD:.3e}", flush=True)
