#!/usr/bin/env python3
"""
mpmath_weil_hiprec.py — LEDGER §104/§108: ARBITRARY-PRECISION Weil operator (the untried unblock).
Every numerical floor in this campaign (§74 |v̂|, §98 eig0, §100 QW(k), §103 cluster, §104 arch↔prime) is a FLOAT64
catastrophic-cancellation floor at ~1e-7..1e-14: Q_arch (O(1)) + Q_prime (O(1)) → the true knife-edge ~1e-14, invisible
in double. We have NEVER run the linear algebra in high precision. This computes, at dps digits:
  - λ_a (lowest gen-eigenvalue of the FULL screw Weil form) and λ_1 (gap) — is λ_a definitively ≥0 below the float64 floor?
  - the ARCH-ONLY vs FULL split (Q_arch, Q_prime separately at full precision, so their cancellation resolves).
  - the ground-state PARITY (even? = Connes even-simple), definitively, not floored.
Sine-Galerkin basis v_k(x)=sin(kπ(x+a)/2a) (Dirichlet, mean-zero v'), same as suzuki_screw_lambda_a but all-mpf.
"""
import sys, time
import mpmath as mp

def screw_g_mp(t, primes_logs, a):
    """Suzuki screw g(t) at high precision. primes_logs = list of (Lambda(n), log n) with n<=e^{2a} (finite)."""
    at = abs(t)
    # boundary/pole  -4(e^{t/2}+e^{-t/2}-2)
    boundary = -4*(mp.e**(t/2)+mp.e**(-t/2)-2)
    # arch core  -(|t|/2)(psi(1/4)-log pi)  = +2.686|t|
    om0 = mp.digamma(mp.mpf(1)/4) - mp.log(mp.pi)   # = -5.3722...
    arch_core = -(at/2)*om0
    # arch tail  -(1/4)[Phi(1,2,1/4) - e^{-|t|/2} Phi(e^{-2|t|},2,1/4)]
    phi0 = mp.zeta(2, mp.mpf(1)/4)                   # Phi(1,2,1/4)=sum (n+1/4)^-2
    if at < mp.mpf(10)**(-20):
        phit = phi0
    else:
        z = mp.e**(-2*at)
        phit = mp.nsum(lambda n: z**n/(n+mp.mpf(1)/4)**2, [0, mp.inf])
    arch_tail = -(mp.mpf(1)/4)*(phi0 - mp.e**(-at/2)*phit)
    # primes  sum_{log n <= |t|} Lambda(n)/sqrt(n) (|t|-log n)
    prime = mp.mpf(0)
    for (Lam, logn, sqrtn) in primes_logs:
        if logn <= at:
            prime += Lam/sqrtn*(at-logn)
    return boundary + arch_core + arch_tail + prime, boundary + arch_core + arch_tail  # (full, arch-only)

def prime_data(a):
    """(Lambda(n), log n, sqrt n) for prime powers n<=e^{2a}."""
    N = int(mp.e**(2*a))+1
    out=[]
    is_c=[False]*(N+1)
    for p in range(2,N+1):
        if not is_c[p]:
            for m in range(2*p,N+1,p): is_c[m]=True
            lp=mp.log(p); pk=p
            while pk<=N:
                out.append((lp, mp.log(pk), mp.sqrt(mp.mpf(pk)))); pk*=p
    return out

def weil_matrices(a, K, M):
    """Full and arch-only Weil matrices Q_{kl}=∫∫ g(x-y) v_k'(x) v_l'(y), plus Mass, at high precision."""
    pl = prime_data(a)
    xs = [(-a + 2*a*i/M) for i in range(M+1)]
    h = 2*a/M
    # g on the difference grid (Toeplitz): d = (i-j)*h, i-j in [-M,M]
    gfull=[None]*(2*M+1); garch=[None]*(2*M+1)
    for d in range(0, M+1):
        gf,ga = screw_g_mp(d*h, pl, a)
        gfull[M+d]=gf; gfull[M-d]=gf; garch[M+d]=ga; garch[M-d]=ga   # even
    # basis derivatives v_k'(x_i) = (kπ/2a) cos(kπ(x_i+a)/2a); values v_k(x_i)=sin(...)
    import math
    def deriv(k,i): return (k*mp.pi/(2*a))*mp.cos(k*mp.pi*(xs[i]+a)/(2*a))
    def val(k,i):   return mp.sin(k*mp.pi*(xs[i]+a)/(2*a))
    w=[h]*(M+1); w[0]=h/2; w[M]=h/2
    # Wp[k][i] = v_k'(x_i)*w_i ; then Q = Wp G Wp^T
    Wp=[[deriv(k,i)*w[i] for i in range(M+1)] for k in range(1,K+1)]
    Vv=[[val(k,i) for i in range(M+1)] for k in range(1,K+1)]
    def assemble(g):
        # (G Wp^T)[i][l] = sum_j g[i-j] Wp[l][j]
        GW=[[mp.fsum(g[M+i-j]*Wp[l][j] for j in range(M+1)) for l in range(K)] for i in range(M+1)]
        Q=[[mp.fsum(Wp[k][i]*GW[i][l] for i in range(M+1)) for l in range(K)] for k in range(K)]
        return Q
    Qf=assemble(gfull); Qa=assemble(garch)
    Mass=[[mp.fsum(Vv[k][i]*w[i]*Vv[l][i] for i in range(M+1)) for l in range(K)] for k in range(K)]
    return mp.matrix(Qf), mp.matrix(Qa), mp.matrix(Mass), Vv, xs, w

def gen_eig_sym(Q, Mass):
    """lowest eigenpairs of Q c = λ Mass c via Cholesky, all mpf."""
    L = mp.cholesky(Mass)
    Linv = L**-1
    A = Linv * Q * Linv.T
    A = (A+A.T)/2
    E, V = mp.eigsy(A)                 # ascending eigenvalues, orthonormal V
    # eigenvectors in original: c = Linv.T * V[:,j]
    return E, (Linv.T * V)

def main():
    dps = int(sys.argv[1]) if len(sys.argv)>1 else 30
    K   = int(sys.argv[2]) if len(sys.argv)>2 else 24
    M   = int(sys.argv[3]) if len(sys.argv)>3 else 600
    alist = [mp.mpf(x) for x in (sys.argv[4].split(',') if len(sys.argv)>4 else ['0.6','0.8','1.0','1.2'])]
    mp.mp.dps = dps
    print(f"# mpmath high-precision Weil operator | dps={dps} K={K} M={M}", flush=True)
    print(f"# {'a':>5} {'lambda_a(FULL)':>26} {'lambda_1(gap)':>16} {'parity(ground)':>14} {'arch_only lam0':>26}", flush=True)
    for a in alist:
        t0=time.time()
        Qf,Qa,Mass,Vv,xs,w = weil_matrices(a, K, M)
        Ef, Vf = gen_eig_sym(Qf, Mass)
        Ea, Va = gen_eig_sym(Qa, Mass)
        # ground-state parity: reconstruct v(x)=Σ c_k sin_k, check symmetry
        c0 = [Vf[k,0] for k in range(K)]
        vx = [mp.fsum(c0[k]*Vv[k][i] for k in range(K)) for i in range(M+1)]
        num = mp.fsum(vx[i]*vx[M-i]*w[i] for i in range(M+1))
        den = mp.fsum(vx[i]*vx[i]*w[i] for i in range(M+1))
        par = num/den
        parity = 'EVEN' if par>0.9 else ('ODD' if par<-0.9 else f'mix({mp.nstr(par,3)})')
        print(f"  {mp.nstr(a,3):>5} {mp.nstr(Ef[0],18):>26} {mp.nstr(Ef[1],8):>16} {parity:>14} {mp.nstr(Ea[0],18):>26}  ({time.time()-t0:.0f}s)", flush=True)

if __name__=='__main__':
    main()
