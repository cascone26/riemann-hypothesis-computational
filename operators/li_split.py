#!/usr/bin/env python3
"""§306: exact arch/prime split of the Li zero-power-sums S_j. From xi=1/2 s(s-1)pi^{-s/2}Gamma(s/2)zeta(s):
  xi'/xi(s) = B + sum_rho[1/(s-rho)+1/rho] = [1/s+1/(s-1)-1/2 ln pi + 1/2 psi(s/2)] + zeta'/zeta(s).
Expanding sum_rho 1/(s-rho) = -sum_{j>=0} s^j S_{j+1} gives  -S_{j+1} = [s^j] G(s),
  G(s) = h(s) + 1/(s-1) - 1/2 ln pi + zeta'/zeta(s),  h(s)=1/2 psi(s/2)+1/s (pole-cancelled arch).
PRIME content of S_{j+1} = -[s^j] zeta'/zeta(s)|_0. Question: is that a convergent {log n} moment, or a
continuation constant? Compute the split; compare arch vs prime magnitude; verify -S_{j+1}=[s^j]G vs zeros. REPORT."""
import mpmath as mp
mp.mp.dps=30
gammas=[mp.mpf(x) for x in open("operators/_zeros150.txt").read().split()]
def Szero(j):
    s=mp.mpf(0)
    for g in gammas:
        rho=mp.mpf(1)/2+1j*g
        s+=rho**(-j)+mp.conj(rho)**(-j)
    return mp.re(s)
# Taylor coefficients at s=0
Kc=7
def zlogder(s): return mp.zeta(s,derivative=1)/mp.zeta(s)
# arch Taylor coeffs in closed form: h(s)=1/2 psi(s/2)+1/s = -gamma/2 + 1/2 sum_{k>=1}(-1)^{k+1} zeta(k+1)(s/2)^k
# arch(s)=h(s)+1/(s-1)-1/2 ln pi ; [s^k]1/(s-1) = -1 all k>=0
ta=[]
for k in range(Kc):
    if k==0: hk=-mp.euler/2
    else: hk=mp.mpf(1)/2*(-1)**(k+1)*mp.zeta(k+1)/2**k
    ak=hk-1 + (-mp.log(mp.pi)/2 if k==0 else 0)
    ta.append(ak)
tz=mp.taylor(zlogder,0,Kc)       # zeta'/zeta Taylor coeffs (PRIME content, via continuation)
print("j :  S_{j+1}(zeros)   -[s^j]G(check)   arch-part -[s^j]arch   prime-part -[s^j](z'/z)",flush=True)
for j in range(1,Kc):
    G_j=ta[j]+tz[j]
    Sz=Szero(j+1)
    print(f"{j}:  {mp.nstr(Sz,8):>12}   {mp.nstr(-G_j,8):>12}   arch={mp.nstr(-ta[j],8):>11}   prime={mp.nstr(-tz[j],8):>11}",flush=True)
print("",flush=True)
print("So S_{j+1} = arch + prime, with prime = -[s^j] zeta'/zeta(0). Now: is prime a CONVERGENT {log n} sum?",flush=True)
# naive convergent von Mangoldt moment (would-be if series converged): sum_n Lambda(n)/n^0 * (log n)^j — DIVERGES.
# Show it diverges (partial sums grow) => prime content is NOT a convergent {log n} moment; needs continuation:
def mang(N):
    L={}
    for p in range(2,N):
        if all(p%d for d in range(2,int(p**.5)+1)):
            q=p
            while q<N: L[q]=mp.log(p); q*=p
    return L
L=mang(2000)
for j in [1,2]:
    partial=[sum(L[n]*mp.log(n)**j for n in sorted(L) if n<X) for X in (100,500,2000)]
    print(f"  would-be prime moment sum_n Lambda(n)(log n)^{j} partial @X=100,500,2000: {[mp.nstr(x,5) for x in partial]}  (DIVERGES => not convergent; continuation needed)",flush=True)
print("DONE",flush=True)
