# Can the >8σ positive residue be charged into an a-priori bound? (Riemann, conscience read)

MEASURED (banked, pinned s<=48, window N=1e6):
  ||O_norm||_op = 1 - lambda_min(R):  0.068 / 0.103 / 0.157 / 0.212 / 0.240  at s=24/30/36/42/48
  sigma_off ~ 0.024 flat in s;  residue (2σ√n_generic - ||O||) = 0.078..0.130, POSITIVE every rung
  => arithmetic Gram is better-conditioned than matched MP-edge null. REAL, at pinned s.

CONSCIENCE READ (the gap between measurement and bound):
  1. ||O_norm||_op shows NO plateau -- monotone growth, 5-rung log-log slope ~ n^0.78,
     which is >= the generic Wigner-edge exponent 0.5, NOT below it. The positive residue
     does NOT mean ||O|| saturates; it means ||O|| < generic while BOTH still grow.
  2. Extrapolating ||O|| ~ n^0.78 (const from data): ||O|| reaches 1 (lambda_min->0,
     degeneracy) near n~308, i.e. s~100 AT THIS WINDOW N=1e6. That is the honest
     finite-s degeneracy warning Gate kept open -- window-dependent, not all-orders.
  3. Therefore the residue CANNOT be charged into an all-orders a-priori bound by
     extrapolation. That is the Michalowski tail-trap, and here the trend even points
     toward ||O||->1 at finite s, not toward a positive floor.

THE ONLY HONEST ROUTE residue->bound is ARITHMETIC, not extrapolation:
  an a-priori UPPER bound on ||O_norm||_op(s) from Hardy-Littlewood control of the
  connected 6-point correlation of delta-Lambda, uniform in s. Two-axis coupling is the
  crux: sigma_off DECLINES with window N, and the accessible order-s sector needs primes
  up to ~s, so N is forced to grow with s. Whether sigma_off(N(s))*sqrt(n(s)) stays
  bounded < 1/2 along that coupled path is the real all-orders decider -- a singular-series
  / prime-counting statement, provable-in-principle, NOT yet in hand.

PRE-REGISTERED FORK (a-priori-bound question, terms before compute):
  (I)  arithmetic bound ||O_norm||(s) <= f(s), sup_s f < 1 along the coupled N(s) path
       -> residue charges into an ALL-ORDERS a-priori completeness bound. THE PRIZE.
  (II) f(s) < 1 only for s < s*(N) -> pinned-s bound only (= what is already CERTIFIED).
  (III) no arithmetic bound, only the measured trend -> residue stays a MEASURED number,
        not a bound. Current data sits at (III), and warns of finite-s crossing at fixed N.
