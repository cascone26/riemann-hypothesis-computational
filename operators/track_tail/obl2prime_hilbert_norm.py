#!/usr/bin/env python3
"""
obl2prime_hilbert_norm.py  -- Rie, obligation 2' (off-diagonal tail bound).

Result chain, all from the deposited Q51 matrix (results/Q51_a0549_nl4000_float.npy):

1. Extract the scalar sequence behind the Loewner collapse Q[k,l]=(kl/(k^2-l^2))(H(l)-H(k))
   via Hawk's cocycle: s(m) = Q[k0,m](k0^2-m^2)/(k0 m) = H(m)-H(k0), and fit
        m*s(m) = D*m + C + B1 cos(w m) + B2 sin(w m),  w = b*log2.
   -> u(m) := m*H(m) = u_c + B1 cos + B2 sin, with u_c = C (Cauchy/gamma amplitude),
      osc amp = hypot(B1,B2) (modulated/beta amplitude at frequency b*log2).
   Two-term law confirmed: u_c ~ -0.907, osc amp ~ 0.624, rms ~2e-3, stable across k0.
   A-free prediction beta/gamma = 4 w_p/pi = 0.624 (w_p = log2/sqrt2) verified to ~9%
   (measured 0.683 -- a stable next-order/coefficient residual, not noise).

2. Clean operator form:  Q[k,l] = [k u(l) - l u(k)] / ((k-l)(k+l))
   = u_c/(k+l) + (modulated channel).  The CONSTANT part is EXACTLY the Cauchy kernel.

3. Extrapolate the verified law to the infinite tail (modes>=25, odd) and measure:
   - ROW SUM diverges (~log M): Gershgorin/row-sum route is DEAD on the infinite tail.
     This kills the v1-draft framing AND every in-window "near+far" assembly.
   - OPERATOR NORM converges (log-slow, from below) to |u_c|*pi/2 = 1.425, the
     ODD-MODE discrete Hilbert norm (pi/2, not pi -- spacing 2 halves it).
   - Margin: lambda_min(T_tail) >= d0 - ||Off|| = 2.057 - 1.425 = 0.632 > 0.
     Closes -- but ~half the in-window number, and only via Hilbert, never row sums.

Remaining obligation 2': pin |u_c| analytically (= log-singularity coeff of H' times
normalization; the A-free ratio already fixes beta once u_c is fixed) and bound the
modulated remainder (numerically o(1) vs the Hilbert norm; ratio->1 below). Then the
theorem margin is d0 - |u_c|*(pi/2), gated on Isa's d0>=2.057 and this u_c.
"""
import numpy as np

Q = np.load('results/Q51_a0549_nl4000_float.npy')
n = Q.shape[0]; modes = np.array([1 + 2*i for i in range(n)])
a = 0.549; b = np.pi/(2*a); w = b*np.log(2.0)
w_p = np.log(2.0)/np.sqrt(2.0)

# --- 1. extract + fit the scalar sequence ---
k0i = 7; k0 = modes[k0i]; m = []; s = []
for j in range(n):
    if j == k0i or modes[j] < 21:
        continue
    s.append(Q[k0i, j]*(k0**2 - modes[j]**2)/(k0*modes[j])); m.append(modes[j])
m = np.array(m, float); s = np.array(s)
D, C, B1, B2 = np.linalg.lstsq(
    np.column_stack([m, np.ones_like(m), np.cos(w*m), np.sin(w*m)]), m*s, rcond=None)[0]
u_c = C
print(f"u_c={u_c:.4f}  osc_amp={np.hypot(B1,B2):.4f}  ratio={np.hypot(B1,B2)/abs(u_c):.3f}"
      f"  pred 4w_p/pi={4*w_p/np.pi:.3f}")
print(f"|u_c|*pi/2 (odd-mode Hilbert asymptote) = {abs(u_c)*np.pi/2:.4f}")

# --- 2/3. build the extrapolated off-diagonal, measure row sum vs op norm ---
def u(mm): return u_c + B1*np.cos(w*mm) + B2*np.sin(w*mm)
def Off(mm):
    K, L = np.meshgrid(mm.astype(float), mm.astype(float), indexing='ij')
    den = (K - L)*(K + L)
    return np.where(K == L, 0.0, (K*u(L) - L*u(K))/np.where(den == 0, 1, den))

print(f"\n{'M':>6} {'max row-sum':>12} {'op-norm':>9}  d0-opnorm")
for Mmax in [51, 201, 801, 3201, 12801]:
    mm = np.arange(25, Mmax + 1, 2); Am = Off(mm)
    rs = np.max(np.sum(np.abs(Am), axis=1))
    op = np.max(np.abs(np.linalg.eigvalsh(Am)))
    print(f"{Mmax:6d} {rs:12.3f} {op:9.4f}  {2.057-op:8.4f}")
print("\nrow sum -> diverges (Gershgorin dead); op norm -> |u_c|*pi/2=1.425; margin -> 0.632")
