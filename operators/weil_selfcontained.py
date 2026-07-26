#!/usr/bin/env python3
"""§272: self-contained high-precision Riemann-Weil quadratic form (mpmath), foundation verified in §271.
Weil form on even test functions phi (work in frequency domain: F(r)=phi_hat(r), real even, r in R):
  Q(phi) = |F(i/2)|^2 + |F(-i/2)|^2            [pole terms, s=0,1]
         + (1/2pi) int F(r)^2 Omega(r) dr       [archimedean, Omega=Re psi(1/4+ir/2)-log pi, VERIFIED §271]
         - 2 Sum_n Lambda(n)/sqrt(n) * g(log n)  [primes, g(u)=(1/2pi) int F(r)^2 e^{-iur} dr = FT of F^2]
RH <=> Q>=0 for all phi. As a quadratic form Q=sum c_j c_k B[j,k] on a basis F_j; compute signature of B.
Basis: F_j(r) = r^{2j} exp(-r^2/(2 sig^2)) (even, well-conditioned, real).  Guard: arch+pole must be pos-def;
full (with primes) is RH => should stay pos-def for a reasonable basis (any NEG => bug or genuine near-failure)."""
import mpmath as mp
mp.mp.dps = 18
def Omega(r): return mp.re(mp.digamma(mp.mpf(1)/4 + 1j*r/2)) - mp.log(mp.pi)

M = 5                      # even Hermite functions, NORMALIZED (orthonormal -> well-conditioned)
# F_j(r) = H_{2j}(r) e^{-r^2/2} / sqrt(2^{2j}(2j)! sqrt(pi))  (norm of H_n e^{-r^2/2} in L^2(dr) is sqrt(2^n n! sqrt pi))
def hnorm(n): return mp.sqrt(2**n * mp.factorial(n) * mp.sqrt(mp.pi))
def F(j, r): return mp.hermite(2*j, r) * mp.e**(-r**2/2) / hnorm(2*j)
def Fhat_at_imag(j, a):   # F_j(i a) = H_{2j}(i a) e^{a^2/2} / norm
    return mp.hermite(2*j, 1j*a) * mp.e**(a**2/2) / hnorm(2*j)

def g_jk(j, k, u):        # (1/2pi) int F_j F_k e^{-iur} dr = FT of product at u (real even)
    return (1/(2*mp.pi)) * 2*mp.quad(lambda r: F(j,r)*F(k,r)*mp.cos(u*r), [0, mp.inf])

def vonmangoldt(N):
    Lam = [mp.mpf(0)]*(N+1)
    for p in range(2, N+1):
        if Lam[p]==0:
            lp = mp.log(p); pk=p
            while pk<=N:
                Lam[pk]=lp; pk*=p
    return Lam

Nmax = 60
Lam = vonmangoldt(Nmax)
print(f"§272 self-contained Weil form, M={M} even-Hermite basis, primes to {Nmax}, dps={mp.mp.dps}", flush=True)

# build the three matrices
Bpole = mp.zeros(M); Barch = mp.zeros(M); Bprime = mp.zeros(M)
for j in range(M):
    for k in range(M):
        Bpole[j,k] = mp.re(Fhat_at_imag(j,mp.mpf(1)/2))*mp.re(Fhat_at_imag(k,mp.mpf(1)/2))*2  # i/2 and -i/2 same (even)
        Barch[j,k] = (1/(2*mp.pi))*2*mp.quad(lambda r: F(j,r)*F(k,r)*Omega(r), [0, mp.inf])
        s = mp.mpf(0)
        for n in range(2, Nmax+1):
            if Lam[n]!=0:
                s += Lam[n]/mp.sqrt(n)*g_jk(j,k,mp.log(n))
        Bprime[j,k] = -2*s
    print(f"  row {j} done", flush=True)

def sig_of(Bmat, label):
    ev = mp.eigsy(Bmat, eigvals_only=True)
    npos=sum(1 for e in ev if e>1e-20); nneg=sum(1 for e in ev if e<-1e-20)
    print(f"  {label:>18}: eigs={[mp.nstr(e,3) for e in ev]}  ->  {npos}+ / {nneg}-", flush=True)
    return ev
print("SIGNATURES:", flush=True)
sig_of(Barch+Bpole, "arch+pole (guard)")
sig_of(Bprime, "prime only")
Bfull = Barch+Bpole+Bprime
sig_of(Bfull, "FULL Weil (=RH)")
print("DONE. arch+pole pos-def = guard OK. FULL pos-def => Weil positivity holds on this basis (RH consistent).", flush=True)
print("Any NEG in FULL => debug construction OR a genuine near-failure direction (inspect).", flush=True)
