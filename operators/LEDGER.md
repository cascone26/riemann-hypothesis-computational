
## §381 (2026-07-28, Riemann7) — CLOSED-FORM FEJER PRIME g_pair BUILT + VALIDATED (l2_fejer_closed.py) — the §380 enabler. HANDLE (validated).
Fejer FT is the triangle T(v)=(π/W)(1−|v|/2W)₊ (compact, |v|<2W). Derived:
  ∫fej(r−a)fej(r−b)e^{iur}dr = (1/2π)e^{iub}·J(u,a−b),  J(u,d)=∫T(v)T(u+v)e^{−ivd}dv,
  gpair(cj,ck,u) = (1/2π)²·Re Σ_{a∈±cj,b∈±ck} e^{iub}J(u,a−b).
J computed by EXACT piecewise (quadratic×exp) antiderivative over the triangle-overlap support (breakpoints at
0, −u, ±2W endpoints; per-piece F=e^{−idv}(p+qv+rv²) coefficient-matched). VALIDATION: closed-form prime
matrix vs l2_fejer's numerical (§376) — **rel err 4.2e-9**; antideriv-J vs direct-integral J — 1e-9. The
prime side is now INSTANT + EXACT (was the §380 infeasibility — ~30 oscillatory acb integrals per entry).
Bug caught+fixed by the validation gate: initial assembly missed one 1/2π factor (closed/num = exactly 2π);
localized by Jexact-matches-Jdirect-but-matrix-doesn't ⟹ assembly error, one-line fix. Verify-discipline win
(the self-validating gate turned a wrong closed form into a right one in one debug pass). ⟹ the rigorous 2nd-
instrument theorem (§380) is now FEASIBLE: rigorous Fejer B_full = pole (closed) + arch (one acb.integral per
entry, sinc⁴×Ω, explicit 1/r⁴ tail) − prime (closed-form triangle sum, instant exact). Next: port J/gpair to
arb intervals and run the interval LDL^T (the arch integral is the only remaining numerical piece, 36 not
36×30). Enabler done; rigorous cert is now a tractable build.
