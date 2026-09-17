#!/usr/bin/env python3
"""§956 UHP-ZERO COUNT of the de Branges characteristic function E(z)=W(0,z) — the correct N₁ test.

§955 ruled out the naive Weyl ratio. The right object: the de Branges structure function E(z). Its negative-
square index κ = number of zeros in the OPEN upper half-plane (Hermite-Biehler ⟺ 0 UHP zeros ⟺ genuine N₀;
each UHP zero is one negative square). In the z-variable (s=½+iz): nontrivial ζ-zeros are REAL under RH; the
completed ξ's trivial factor s(s−1) puts zeros at z=±i/2, so s=0 gives ONE UHP zero at z=+i/2. PREDICTION under
RH: E(z)=W(0,z) (exact Suzuki characteristic function, real-axis zeros already verified = ζ-zeros in §242 main)
has EXACTLY ONE UHP zero, near z=i/2 — the single negative square = ζ's pole shadow. Count via the argument
principle (1/2πi)∮ E'/E dz around a UHP rectangle (avoiding the real axis). If =1 near i/2, the ζ de Branges
function is N₁ as the Krein-Langer route requires. HANDLE (own code). Honest re finite-K spurious zeros: report
the location of every UHP zero found, not just the count."""
import numpy as np, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from suzuki_exact_Q import W_family

def count_zeros_rect(f, x0, x1, y0, y1, n=1600):
    """(1/2πi)∮ f'/f via f'/f ~ d/dz log f; integrate log-derivative numerically around rectangle."""
    xs = np.linspace(x0, x1, n); ys = np.linspace(y0, y1, n)
    # build closed contour (counterclockwise): bottom, right, top, left
    pts = np.concatenate([xs + 1j*y0, x1 + 1j*ys, xs[::-1] + 1j*y1, x0 + 1j*ys[::-1]])
    vals = np.array([f(z) for z in pts])
    # winding number of f(contour) around 0 = (1/2π) total change of arg
    ang = np.unwrap(np.angle(vals))
    return (ang[-1] - ang[0] + (np.angle(vals[0]) - np.angle(vals[-1])))  # close the loop
def winding(f, x0, x1, y0, y1, n=2000):
    xs = np.linspace(x0, x1, n); ys = np.linspace(y0, y1, n)
    pts = np.concatenate([xs + 1j*y0, x1 + 1j*ys, xs[::-1] + 1j*y1, x0 + 1j*ys[::-1], [x0+1j*y0]])
    vals = np.array([f(z) for z in pts])
    ang = np.unwrap(np.angle(vals))
    return int(round((ang[-1] - ang[0])/(2*np.pi)))

def locate_uhp_zeros(f, x0, x1, y0, y1, gx=140, gy=60):
    """coarse |f| minima scan in the UHP box to locate zeros (for honesty about where they are)."""
    xs = np.linspace(x0, x1, gx); ys = np.linspace(y0, y1, gy)
    found = []
    A = np.array([[abs(f(x+1j*y)) for x in xs] for y in ys])
    med = np.median(A)
    for i in range(1, gy-1):
        for j in range(1, gx-1):
            if A[i,j] < A[i-1,j] and A[i,j] < A[i+1,j] and A[i,j] < A[i,j-1] and A[i,j] < A[i,j+1] and A[i,j] < 0.05*med:
                found.append((round(xs[j],2), round(ys[i],3), float(A[i,j])))
    return found

print("="*94)
print("§956 UHP-zero count of the de Branges characteristic function E(z)=W(0,z). RH ⟹ expect ONE (at z≈i/2).")
print("="*94)
print(f"{'a':>5} {'K':>4} | {'winding # (UHP zeros in box)':>28} | UHP zero locations (x, y, |E|)")
for a in [1.0, 1.5, 2.0]:
    for K in [80, 110]:
        W, lam_a = W_family(a, K, NG=400)
        E = lambda z: W(0.0, z)
        # box: x across the spectral window (avoid the far tails), y from just above axis to 1.5
        xlo, xhi = 3.0, 8.0*a + 6.0
        wn = winding(E, xlo, xhi, 0.06, 1.5, n=2400)
        zeros = locate_uhp_zeros(E, xlo, xhi, 0.08, 1.5)
        near_half = [z for z in zeros if abs(z[1]-0.5) < 0.25]
        loc = "; ".join(f"({x},{y})" for x,y,_ in zeros[:6]) if zeros else "none found in scan"
        print(f"{a:>5.1f} {K:>4} | {wn:>28} | {loc}")
print("="*94)
print("READ: winding # = number of E-zeros in the UHP box. RH prediction = 1 (the z=i/2 shadow of ζ's s=0/pole).")
print("If ≈1 and located near y≈0.5, ζ's de Branges function is N₁ — the negative square IS the pole, confirming")
print("the Krein-Langer object directly. >1 or scattered ⟹ finite-K spurious zeros (report honestly) or the box")
print("caught real-axis zeros — refine. This is the direct de Branges realization of the §946-955 Π₁ finding.")
print("="*94)
