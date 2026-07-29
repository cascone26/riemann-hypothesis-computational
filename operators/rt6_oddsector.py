#!/usr/bin/env python3
"""§395 (angle 1): odd/complex-h sector — is the non-even part of the explicit formula identically
degenerate for zeta (self-dual)? Build the ANTISYMMETRIC probe A_c(r)=phi(r-c)-phi(r+c) (odd), form
its zero-side Gram and analytic side. If self-duality forces this sector trivial, closes in one run."""
import numpy as np
import l2_edge_sweep as L
w = 3.0
def Aodd(cs, r):  # antisymmetric probe row at r
    return np.exp(-((r[None,:]-cs[:,None])**2)/(2*w*w)) - np.exp(-((r[None,:]+cs[:,None])**2)/(2*w*w))
print("ODD-SECTOR: antisymmetric probes A_c=phi(r-c)-phi(r+c). Zero-side Gram 2*sum A(g)A(g):")
print(f"{'s':>5} {'n':>3} {'lmin(oddZero)':>13} {'||oddZero||':>11} {'lmin(oddOvl)':>12}")
for s in (5.0, 4.0, 3.0, 2.5):
    cs = np.arange(14.0, 44.0+1e-9, s)
    Az = Aodd(cs, L.ZEROS)
    oddZ = 2.0 * Az @ Az.T
    P = Aodd(cs, L.RGRID); dr = L.RGRID[1]-L.RGRID[0]
    oddOvl = P @ P.T * dr
    print(f"{s:5.1f} {len(cs):3d} {float(np.linalg.eigvalsh(oddZ)[0]):+13.3e} "
          f"{np.linalg.norm(oddZ,2):11.3e} {float(np.linalg.eigvalsh(oddOvl)[0]):+12.3e}", flush=True)
print("\nreading: oddZero identically ~0 (all eigenvalues tiny) => odd sector trivial for self-dual zeta,")
print("sector closes. oddZero O(1) with its own collapse => a NEW independent channel (odd explicit formula).")
