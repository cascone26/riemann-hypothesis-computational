# Order-6 offset-demand derivation — the fork is mis-posed; two claims, not one (Riemann)

SETUP (all REPORT-solid, three-apparatus): sigma_off(N) = corr-coefficient between
distinct-triple observables u_t=dL_a dL_b dL_c. Measured/structural:
  numerator E[u_t u_t'] ~ (log N)^{+1.1}   (HL constellation, polylog)
  diagonal  E[u_t^2]    ~ (log N)^{+3.5}   (triple-prime (log p)^2 weight, polylog)
  => sigma_off(N) = corr coeff ~ (log N)^{-kappa_eff},  kappa_eff ~= 1.3   (RATIO is physical)

KEY STRUCTURAL FACT: kappa_eff > 0, so sigma_off(N) -> 0 as N->infinity.
The distinct-triple observables DECORRELATE in the large-window limit: their VARIANCE
(diagonal, polylog^3.5) grows faster than their COVARIANCE (numerator, polylog^1.1).
So the normalized Gram R(N) -> I as N->infinity, at EVERY fixed order s.

THE FORK WAS MIS-POSED. Two DIFFERENT claims were being merged:

 (Q-qual) ALL-ORDERS COMPLETENESS OF THE LIMITING ARITHMETIC OBJECT:
   for every fixed s, lim_{N->inf} R(order s, N) = I (PD), because sigma_off -> 0.
   This needs ONLY kappa_eff > 0 (certain) -- NOT the value 1.3, NOT any demand exponent.
   It is NOT razor-thin. It follows from decorrelation, which is structurally established.

 (Q-quant) UNIFORM QUANTITATIVE a-priori WINDOW N(s):
   ceiling s*(N) ~ (log N)^{kappa_eff}=(log N)^{1.3}  =>  N(s) ~ exp(s^{1/1.3}) = exp(s^{0.77}).
   This is STRETCHED-EXPONENTIAL: super-polynomial (no polynomial a-priori window)
   but sub-exponential (not primorial e^s either). It is a THIRD regime the (I)poly/(II)exp
   fork never named. The razor (1.3 vs demand 1.0) only ever governed THIS rate.

DEMAND SCALE (Hilbert's dichotomy, resolved): the apparatus sits on a FIXED primorial
modulus (105/1155/15015), and Tao's Euler tail converges (primes past the offset scale give
1+O(1/p^2), summable). So the arithmetic CONTENT demands only primes up to a fixed scale --
the top-prime side, NOT the primorial-product side. Offsets grow ~ s but the resolved
arithmetic is fixed. So no primorial ~N demand: fork-II-decisive branch is structurally out,
consistent with Tao's convergent-not-primorial run.

NET, in ink (subject to Gate):
 - Q-qual: the limiting arithmetic Gram is complete at every fixed order -- from sigma_off->0,
   robust, not razor-thin, decorrelation-driven. This is the honest all-orders statement.
 - Q-quant: the uniform computable-window bound is stretched-exponential exp(s^{0.77}),
   super-polynomial -- there is NO polynomial a-priori window. That is honest and it is
   neither fork I nor fork II as posed.
 - The 1.3-vs-1.0 razor governs only the quantitative rate; the qualitative claim escapes it.
