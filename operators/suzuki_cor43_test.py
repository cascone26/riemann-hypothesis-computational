"""
Numerical test of Suzuki 2301.00421 Theorem 4.2 / Corollary 4.3 (eq 4.4/4.7), OCR-free.

Under RH, the screw-line identity gives
    (1/2pi) * || S_t ||^2_{L2(R)}  =  -g(t)  =  (1/2) G_g(t,t)  =  (1/2) sum_{gamma in Gamma} (2 - 2 cos(gamma t)) / gamma^2
with Gamma = {+/- gamma_n} (eq 4.6, m_gamma=1).  Both sides are clean (no garbled (1.6)/(4.3) OCR):
  LHS: our validated S_t evaluator (suzuki_screwline_St), norm computed on the real line.
  RHS: zero-sum over zeta ordinates.

For real z: E#_xi(z) = conj(E_xi(z)), so
    S_t(z) = (i/2)(1 + E_xi/conj(E_xi)) P_t   ==>   |S_t(z)|^2 = Re(E_xi(z))^2 / |E_xi(z)|^2 * |P_t(z)|^2.
One xi, xi' evaluation per z.  This tests the ENTIRE screw-line construction (Thm 4.2) end to end.
"""
import mpmath as mp
from suzuki_screwline_St import E_xi, P_t, gammas

mp.mp.dps = 20

def abs2_St(z, t, N):
    z = mp.mpf(z)
    E = E_xi(z)                      # complex
    reE = mp.re(E)
    absE2 = mp.re(E)**2 + mp.im(E)**2
    P = P_t(z, t, N)
    absP2 = mp.re(P)**2 + mp.im(P)**2
    return (reE*reE/absE2) * absP2

def lhs_norm(t, N=100, Z=80.0, npts=1600):
    """(1/2pi) integral_{-Z}^{Z} |S_t(z)|^2 dz  (Simpson).  |S_t|^2 ~ 1/z^2 tail, add analytic-ish tail est."""
    t = mp.mpf(t)
    h = mp.mpf(2*Z)/npts
    s = mp.mpf(0)
    prev = abs2_St(-Z, t, N)
    # Simpson on uniform grid
    acc = prev
    vals = [prev]
    for k in range(1, npts+1):
        z = -Z + k*h
        v = abs2_St(z, t, N)
        vals.append(v)
    total = vals[0] + vals[-1]
    for k in range(1, npts):
        total += (4 if k % 2 == 1 else 2) * vals[k]
    integ = total * h / 3
    return integ / (2*mp.pi), vals

def rhs_zerosum(t, N=100):
    t = mp.mpf(t)
    acc = mp.mpf(0)
    for g in gammas(N):
        # pair +g and -g: each contributes (2-2cos(g t))/g^2 ; -g gives same
        acc += 2 * (2 - 2*mp.cos(g*t)) / (g*g)
    return acc / 2   # -g(t) = (1/2) G_g(t,t)

if __name__ == "__main__":
    print("Suzuki Cor 4.3 test:  (1/2pi)||S_t||^2  vs  -g(t)=(1/2)sum(2-2cos gt)/g^2")
    print(f"{'t':>6} {'LHS (norm/2pi)':>20} {'RHS (zero-sum)':>20} {'ratio':>12} {'tail?':>8}")
    N = 120
    for t in [1.0, 2.0, 3.0, 4.0]:
        lhs, vals = lhs_norm(t, N=N, Z=90.0, npts=1800)
        rhs = rhs_zerosum(t, N=N)
        ratio = lhs/rhs if rhs != 0 else mp.inf
        tailmag = f"{mp.nstr(vals[-1],2)}"
        print(f"{t:6.2f} {mp.nstr(lhs,12):>20} {mp.nstr(rhs,12):>20} {mp.nstr(ratio,8):>12} {tailmag:>8}")
