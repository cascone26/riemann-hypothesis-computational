import numpy as np, math
from scipy.special import digamma
import operators._dh_negative_control as DH   # reuse lambda_D, omega_D, vm, a_coef

def build_tau(omega_fn, coef_arr, Ncut, N, fourW, a, b):
    L = fourW/4
    RG = np.arange(-(b+80), (b+80), 0.004); OM = omega_fn(RG)
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
    return arch - prime   # DH is ENTIRE: no pole/boundary term

def kappa(tau, thr=1e-6):
    ev = np.linalg.eigvalsh(tau)
    return int((ev < -thr).sum()), ev[0], np.max(np.abs(ev))

Ncut=3000
LD = DH.lambda_D(Ncut)
print(f"DH entire Weil form: kappa = #(eigenvalues < -1e-6) [genuine negative squares, not machine noise]")
print(f"{'window':>12} {'4W':>4} {'N':>4} {'kappa(DH)':>10} {'mu_min':>12}")
# window WITHOUT off-line zeros (all on-line) vs WITH (off-line at t~85.7)
for (a,b,fW,N) in [(14,44,5.0,60),(50,80,8.0,90),(78,92,11.0,110)]:
    tau = build_tau(DH.omega_D, LD, Ncut, N, fW, a, b)
    k,mu,nrm = kappa(tau)
    print(f"[{a:>3},{b:>3}]   {fW:>4} {N:>4} {k:>10} {mu:>12.4e}", flush=True)
print("\nPREDICTION: kappa=0 where all zeros on-line; kappa>=1 at [78,92] where DH's off-line zeros live.")
print("If kappa jumps 0->>=1 at the off-line window, mu/kappa DETECTS the RH-violation.")
