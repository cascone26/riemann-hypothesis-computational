
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

## §463 (BUILD attempt 2 — correct-extension shooting; does NOT reproduce zeros either; W_sa singular-point matching needs the full CCM paper) — shooting on |x|>1 with the derived BCs (regular Frobenius start at x=1: ξ(1)=1,ξ'(1)=(μ-4π²)/2; pure-sin B=0 at ∞) gives E=√|μ| ≈ 0.33-0.41×γ_k with ~3× excess density. Verify-caught. HANDLE (honest, no false success).
rt6_prolate_shoot*.py. Derived the W_sa boundary conditions correctly from CM arXiv:2112.05500: L_{β±}=0 (regular/
flux→0 at ±1, giving the Frobenius start ξ'(1)=(μ-4π²)/2) and L_{β̂}=0 ⟺ B=0 (pure sin(2πx)/x at ∞, verified from
eq17: the functional → -2πB). Shot on |x|>1 (where negative eigenvalues live) for B(μ)=0. RESULT: E-values 4.6,7.2,
9.2,11.0,12.7,14.2,... vs zeros 14.1,21.0,25.0,30.4,... — ratio 0.33→0.41 (NOT 1), and ~3× too many (20 roots vs 4
zeros in [3,32]). So this extension attempt ALSO fails (verify-caught, like §462's naive one). The persistent 3×
excess density + ~0.4 ratio ACROSS both methods = a SYSTEMATIC issue: the correct W_sa requires the singular-point
MATCHING between [-1,1] and |x|>1 (the even-function connection through the regular-singular points ±1, with the
discontinuity-allowed extension), which couples the two sides — my free-start x>1 shooting doesn't enforce it. The
exact ζ-connection lives in the paper's LATER sections (Borel summation of the finite-solution Taylor series
U_n(μ), Lemma 3.1 → the FT of φ_μ → ζ) which I don't have from the first 8 pages. HONEST STATE: I obtained the EXACT
operator, derived the correct BCs, and made two careful implementation attempts — BOTH verify-caught as not
reproducing the zeros. Building it correctly = a genuine substantial task needing the full paper (the singular-point
matching + the Borel-summation ζ-connection). Did NOT claim false success (Charisma: no Handle for "realized zeros").
