"""Does the Weil explicit-formula quadratic form share the sign theorem's diagonal-dominance?
Weil: sum_rho h(gamma) = [archimedean + pole] - (1/2pi) sum_p sum_k (log p) p^{-k/2} 2 g(k log p).
For a test function, the quadratic form Q(g) = sum_rho |ghat(gamma)|^2 (>=0 iff RH).
Decompose: DIAGONAL (archimedean/main, from the density) vs OFF-DIAGONAL (prime terms).
Does off-diagonal/diagonal ~ 1/log(scale)? (the sign-theorem mechanism, tested on the Weil form)"""
import numpy as np
from sympy import primerange
# Use a family of test functions g_a(u) = indicator-ish (Gaussian width sigma centered 0), vary sigma.
# Weil functional W(sigma) = archimedean(sigma) - prime_sum(sigma); RH <=> W>=0 for all admissible g.
# archimedean (diagonal-like): ~ integral involving log; prime_sum (off-diagonal): explicit.
def weil_pieces(sigma):
    # g(u) = e^{-u^2/(2 sigma^2)}, ghat(t) = sigma sqrt(2pi) e^{-sigma^2 t^2/2}
    # archimedean main term ~ (1/2pi) int ghat(t) [Re psi(1/4+it/2) - log pi] dt  (approx: log-density)
    # simplest proxy: A(sigma) = g(0) x <log gamma density weight> -- use the leading (1/2pi) int ghat Omega
    # For scaling purposes compute the PRIME SUM (off-diagonal) and the ARCH via zero-side.
    from scipy.special import loggamma as slg
    ts = np.linspace(-60/sigma, 60/sigma, 4000)
    ghat = sigma*np.sqrt(2*np.pi)*np.exp(-sigma**2*ts**2/2)
    dt = ts[1]-ts[0]
    # archimedean density Omega(t) = Re[psi(1/4+it/2)]/... ~ (1/2)log(t^2/4pi+..) ; use (1/2)log((|t|+2)/(2pi))
    Omega = 0.5*np.log((np.abs(ts)+2)/(2*np.pi))
    A = np.trapezoid(ghat*Omega, ts)/(2*np.pi)*2
    # prime sum
    P = 0.0
    for p in primerange(2, int(np.exp(8*sigma))+3):
        lp=np.log(p); k=1
        while k*lp < 8*sigma+2:
            gk = np.exp(-(k*lp)**2/(2*sigma**2))  # g(k log p)
            P += 2*lp*p**(-k/2)*gk/(2*np.pi)
            k+=1
    return A, P
print(" sigma   arch(diag)   primes(offdiag)   off/diag   W=diag-off")
for sigma in (0.3, 0.5, 1.0, 2.0, 3.0, 4.0):
    A,P = weil_pieces(sigma)
    print(f" {sigma:.1f}    {A:.4f}      {P:.4f}       {P/A:.3f}     {A-P:+.4f}", flush=True)
print("DONE")
