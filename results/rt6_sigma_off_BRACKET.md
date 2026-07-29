# rt6 sigma_off ladder — HONEST BRACKET, not a pinned decider (Riemann, 2026-07-28)

TASK (Dyson's all-orders decider): is sigma_off(s)*sqrt(n) bounded < ~0.5 (=> lambda_min(R)>0
all-orders) or does it grow like the generic bulk-edge?

RESULT: my singular-series kernel O CANNOT supply this number self-consistently.
Two normalizations of R = O / sqrt(D_t D_t') bracket the truth, neither valid:
  - D = m2^3 * Wself  (m2=13.2): sigma_off*sqrt(n) = 0.095..0.244 (bounded) BUT lambda_min
    crosses 0 (-0.42 @ s=54) -> IMPOSSIBLE for a true Gram (Hilbert: mean of a square >=0),
    so this normalization under-scales the off-diagonal (extra ~2300 factor). ARTIFACT.
  - D = Wself (m2-free): |R_offdiag| ~ 55 >> 1 -> violates Cauchy-Schwarz for any Gram. INVALID.

ROOT CAUSE: O (dimensionless singular-series density ratio) has the SIGN + LOCATION of the
connected correlation (3-apparatus confirmed) but NOT the field-moment SCALE of the connected
6-point moment E[u_t u_t']. That scale is neither m2^3 (overshoots) nor 1 (undershoots).

CONSEQUENCE / ROUTE:
  - The valid sigma_off(s) must be read from Hilbert's DIRECT Gram G[t,t']=E[u_t u_t'] computed
    in physical delta-Lambda units (Cauchy-Schwarz-respecting by construction), NOT from my proxy.
  - What my kernel DOES certify (banked): sign structure + location + shape-mismatch pattern of
    the off-diagonal. The WHERE, not the magnitude.
  - Dyson: take sigma_off(s) from Hilbert's G, not from me. Then sigma_off*sqrt(n) vs bulk-edge
    is a physically-normalized decider.
