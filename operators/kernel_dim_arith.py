#!/usr/bin/env python3
"""ATTACK A11: is dim(near-kernel of tau) a clean ARITHMETIC count of the primes <= c?
If dim = f(#prime-powers, pi(c), ...) exactly, that's a new identity tying the operator's
RH-content directly to the primes -- a structural result, not a reshuffle."""
import numpy as np, sys
sys.path.insert(0,"operators")
from critical_eigvec_id import build_tau_np, prime_powers_up_to

def pi_count(c):
    n=0;x=2
    while x<=c:
        if all(x%p for p in range(2,int(x**0.5)+1)): n+=1
        x+=1
    return n

print(f"{'c':>4} {'ker_dim':>8} {'ppw<=c':>7} {'pi(c)':>6} {'2*ppw+1':>8} {'band_M':>7}")
rows=[]
for c in range(2,31):
    tau,ns=build_tau_np(96,c)
    ev=np.linalg.eigvalsh(tau)
    kd=int(np.sum(np.abs(ev)<1e-6))
    ppw=len(prime_powers_up_to(c)); pic=pi_count(c)
    M=(kd-1)//2
    print(f"{c:>4} {kd:>8} {ppw:>7} {pic:>6} {2*ppw+1:>8} {M:>7}")
    rows.append((c,kd,ppw,pic))
# hunt for a clean relation kd = a*ppw + b or kd = 2*pi(c)+... etc
import numpy as np
c=np.array([r[0] for r in rows]);kd=np.array([r[1] for r in rows])
ppw=np.array([r[2] for r in rows]);pic=np.array([r[3] for r in rows])
print("\nfit kd ~ a*ppw+b:", np.polyfit(ppw,kd,1))
print("fit kd ~ a*pi(c)+b:", np.polyfit(pic,kd,1))
print("fit kd ~ a*c+b:", np.polyfit(c,kd,1))
print("kd - 2*ppw:", (kd-2*ppw).tolist())
print("kd - 2*pi(c):", (kd-2*pic).tolist())
