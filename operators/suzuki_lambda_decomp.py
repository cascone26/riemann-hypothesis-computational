"""
Decompose lambda_a on its ACTUAL minimizer w* into archimedean vs prime parts (eq 4.5).
Q_full = int int g(x-y) w'(x)w'(y);  Q_arch uses g_no_primes (arch+pole);  Q_prime = Q_full - Q_arch.
lambda_a = min eig of Q_full/a with eigenvector c* -> w*. Evaluate arch-Rayleigh and prime-Rayleigh on w*.
Shows how the archimedean surplus vs prime deficit balances at the knife-edge as a grows.
Uses the VALIDATED screw_g / screw_g_no_primes from suzuki_screw_lambda_a.
"""
import numpy as np, time
from suzuki_screw_lambda_a import screw_g, screw_g_no_primes, vonmangoldt_upto

def build(a, M, K):
    Nmax=int(np.floor(np.exp(2*a)))+1
    Lam=vonmangoldt_upto(min(Nmax,30_000_000))
    x=np.linspace(-a,a,M+1); h=x[1]-x[0]
    diffs=np.arange(-M,M+1)*h
    gf=np.array([screw_g(d,Lam,Nmax) for d in diffs])
    ga=np.array([screw_g_no_primes(d,Lam,Nmax) for d in diffs])
    idx=np.arange(M+1)
    shift=(idx[:,None]-idx[None,:])+M
    Gf=gf[shift]; Ga=ga[shift]
    kk=np.arange(1,K+1)
    arg=np.pi*np.outer((x+a),kk)/(2*a)
    Vp=(kk*np.pi/(2*a))*np.cos(arg)
    w=np.full(M+1,h); w[0]=w[-1]=0.5*h
    Wp=Vp*w[:,None]
    Qf=Wp.T@Gf@Wp; Qf=(Qf+Qf.T)/2
    Qa=Wp.T@Ga@Wp; Qa=(Qa+Qa.T)/2
    return Qf,Qa

def decomp(a,M=4000,K=80):
    Qf,Qa=build(a,M,K)
    ev,V=np.linalg.eigh(Qf/a)     # generalized w/ mass a I -> divide by a; mass matrix = a*I for sine basis? see note
    c=V[:,0]                      # minimizer coeffs (lowest)
    mass=a*np.dot(c,c)            # ||w||^2 = a * sum c_k^2
    lam=ev[0]
    arch=(c@Qa@c)/mass
    prime=(c@(Qf-Qa)@c)/mass
    return lam, arch, prime

if __name__=="__main__":
    print("lambda_a decomposed on its minimizer w*: lambda = arch_part + prime_part  (eq 4.5 arch vs prime)")
    print(f"{'a':>6} {'lambda_a':>13} {'arch part':>13} {'prime part':>13} {'check sum':>13}")
    t0=time.time()
    for a in [0.2,0.35,0.5,0.75,1.0,1.5,2.0,3.0]:
        M=min(8000,max(3000,int(600/a)))
        lam,arch,prime=decomp(a,M=M,K=80)
        print(f"{a:6.2f} {lam:13.6e} {arch:13.6e} {prime:13.6e} {arch+prime:13.6e}",flush=True)
    print(f"done {time.time()-t0:.1f}s")
    print("\nNote: 'arch part' includes -log a-(2A+1)+L(w*)/||w*||^2 (it IS the g_no_primes Rayleigh quotient).")
