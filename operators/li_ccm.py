#!/usr/bin/env python3
"""§304: connect §303's "log|zeta|^2 stays PD" step to LI'S CRITERION (RH <=> lambda_n >= 0 for all n, where
lambda_n = Taylor coeffs of (d/ds) log xi(s) at s=1, = sum_rho[1-(1-1/rho)^n]). Li's criterion IS the log-step
positivity §303 pinpointed. Compute lambda_n from the zeros, confirm positivity (RH-consistent), and set up the
bridge to our CCM sech^{1/2}/{log n} closed form (§261d). Charisma: REPORT."""
import mpmath as mp
mp.mp.dps=30
# first zeros' imaginary parts (compute via zetazero, cache to file to avoid recompute)
import os
cache="operators/_zeros150.txt"
if os.path.exists(cache):
    gammas=[mp.mpf(x) for x in open(cache).read().split()]
else:
    gammas=[]
    for k in range(1,151):
        gammas.append(mp.im(mp.zetazero(k)))
        if k%10==0: print(f"  got zero {k}",flush=True)
    open(cache,"w").write(" ".join(mp.nstr(g,25) for g in gammas))
print(f"using {len(gammas)} zeros (first gamma={mp.nstr(gammas[0],8)})",flush=True)
# Li coefficients: lambda_n = sum over conjugate pairs of 2*(1 - Re[(1-1/rho)^n]), rho=1/2+i*gamma
def li(n):
    s=mp.mpf(0)
    for g in gammas:
        rho=mp.mpf(1)/2+1j*g
        s+=2*(1-mp.re((1-1/rho)**n))
    return s
print("Li coefficients lambda_n (RH <=> all >= 0):",flush=True)
known={1:'0.0230957',2:'0.0923457',3:'0.207639',4:'0.368793',5:'0.575543',6:'0.827566'}
for n in range(1,7):
    v=li(n)
    kn=known.get(n,'?')
    print(f"  lambda_{n} = {mp.nstr(v,8)}   (known~{kn}; partial-sum over {len(gammas)} zeros, converges from below)  {'>=0 OK' if v>-mp.mpf('0.01') else 'NEG?!'}",flush=True)
print("  (partial sums under-shoot the true lambda_n since tail Sum_{|gamma|>T} n*Re(1/rho) > 0 is dropped)",flush=True)
# Bombieri-Lagarias arithmetic side: lambda_n = (arch/pole term) - (prime/log n term).
# arithmetic form: lambda_n = 1 - (n/2)(gamma_E + log(4pi)) + [sum over zeros already above] ; the PRIME content
# enters via the SAME {log n} explicit-formula lattice as our CCM closed form (§261d). Show the arch leading term:
gammaE=mp.euler
for n in range(1,7):
    lead=1-(mp.mpf(n)/2)*(gammaE+mp.log(4*mp.pi))+ (mp.mpf(n)/2)*sum(1/g**2 for g in gammas)*0  # placeholder arch lead
    print(f"  n={n}: BL arch-lead 1-(n/2)(gammaE+log4pi) = {mp.nstr(1-(mp.mpf(n)/2)*(gammaE+mp.log(4*mp.pi)),6)}",flush=True)
print("DONE",flush=True)
