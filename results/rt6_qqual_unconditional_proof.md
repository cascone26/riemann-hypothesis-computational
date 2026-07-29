# Q-qual is UNCONDITIONAL — proof via one-sided sieve bounds (Riemann)
# CATCH: the "kappa_cov = O(1) Hardy-Littlewood constant" framing is HL-conditional
# (exact 6-pt von Mangoldt correlation = HL k-tuple asymptotic; k=2 is twin primes).
# FIX: replace exact constants with ONE-SIDED bounds; both are unconditional.

CLAIM (Q-qual). Fix an order s (finite admissible-triple set, |cell| = n(s) fixed).
Let G[t,t'] = E_{n<=N}[ u_t u_t' ], u_t = prod_i dL(n+t_i), dL = Lambda - 1,
R = D^{-1/2} G D^{-1/2} (unit diagonal). Then R(N) -> I and lambda_min(R(N)) -> 1
as N -> infinity. Unconditionally. Hence the limiting order-s Gram is PD (complete).

INPUTS (all PNT / Chebyshev / Selberg-sieve grade — NO Hardy-Littlewood):
 (P1) PNT: #{p <= N} ~ N/log N ; psi(N) ~ N ; prime powers p^{k>=2} number O(sqrt N).
 (P2) Selberg upper sieve (Halberstam-Richert, Thm 2.x): for a fixed admissible
      k-tuple h, #{n<=N : n+h_1,...,n+h_k all prime} <= C_k * S(h) * N/(log N)^k,
      C_k absolute, unconditional. (Only the matching LOWER bound needs HL.)

LEMMA 1 (diagonal lower bound). D_t = E[ prod_{i=1}^3 (Lambda(n+t_i)-1)^2 ] >= (1/4)log N (1+o(1)).
 Proof. Every summand is a product of squares >= 0. Restrict the average to
 A = { n<=N : n+t_1 prime, n+t_2 and n+t_3 not prime powers }.
 On A: summand = (log(n+t_1)-1)^2 * (0-1)^2 * (0-1)^2 = (log(n+t_1)-1)^2.
 |A| >= #{n+t_1 prime} - #{n+t_1,n+t_2 both prime} - #{n+t_1,n+t_3 both prime}
       - O(sqrt N)   [prime-power exclusions]
     = N/log N (1+o(1)) - o(N/log N)      [ (P1) then (P2) k=2 ]
     ~ N/log N.
 Restrict further to n+t_1 >= sqrt N (drops O(sqrt N) = o(N/log N)); there
 (log(n+t_1)-1)^2 >= ((1/2)log N - 1)^2 ~ (1/4)(log N)^2.
 So D_t >= (1/N)(N/log N)(1/4)(log N)^2 (1+o(1)) = (1/4) log N (1+o(1)).  QED

LEMMA 2 (off-diagonal upper bound). For distinct 6-offset (t,t'), |C_{t,t'}| = O(1),
 bounded independent of N.
 Proof. C = E[ prod_{i=1}^6 dL(n+h_i) ] = sum_{S subset of 6} (-1)^{6-|S|} E[ prod_{i in S} Lambda(n+h_i) ].
 For each S: prod_{i in S} Lambda(n+h_i) != 0 only if all {n+h_i}_{i in S} are prime powers;
 dropping p^{k>=2} (O(sqrt N), negligible) and bounding each Lambda <= log N,
   E[prod_S Lambda] <= (log N)^{|S|} * #{n: all S prime}/N <= (log N)^{|S|} * C_{|S|} S(S)/(log N)^{|S|}
                      = C_{|S|} S(S) = O(1)          [ (P2) ].
 (S nonadmissible => count 0.) Sum of 2^6 terms each O(1) => |C| = O(1).  QED

THEOREM (Q-qual, unconditional). |R[t,t']| = |C_{t,t'}| / sqrt(D_t D_{t'})
 <= O(1) / ((1/4) log N) = O(1/log N) -> 0  for every fixed distinct pair.
 R(N) is a fixed n(s)-dim symmetric matrix, unit diagonal, off-diagonals -> 0,
 so R(N) -> I entrywise; eigenvalues continuous => lambda_min(R(N)) -> 1 > 0.
 Hence for each fixed s there is N_0(s) with R PD for all N > N_0(s).  QED

SCOPE (honest, matches Q-quant): N_0(s) is INEFFECTIVE (sieve constants + o(1)'s) and
grows fast in s -- consistent with the Q-quant stretched-exponential window exp(s^0.77).
The QUALITATIVE limit is unconditional; the EFFECTIVE window is not polynomial. No conflict.

WHAT IS NOT CLAIMED: the exact decay exponent kappa_eff ~ 1.3 (that IS HL-flavored, the
matching sieve LOWER bound). Q-qual needs only R->0, i.e. the UPPER bound O(1/log N),
which the one-sided sieve gives unconditionally. kappa_eff stays REPORT, not theorem.
