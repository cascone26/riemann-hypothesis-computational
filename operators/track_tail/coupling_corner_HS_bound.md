# The coupling corner bound: ‖B_Cauchy‖ ≤ ½√(log 2) — provable, cut-independent, both parities

**Status (2026-07-24, roundtable):** closes the open derivation Ein named as the M* optimizer's
"real content" (the provable rectangle coupling bound). Verified numerically; the bound is within
3–7% of the measured norm (near-sharp because the corner operator is almost rank-one).
Script: `coupling_corner_HS_bound_check.py` (= tmp corner_hs_check), output in
`results/corner_hs_check_out.txt`.

## Statement

Let B be the Cauchy corner on a parity grid (spacing 2): B_{kl} = 1/(k+l), k ∈ [25, M] (grid
modes), l > M (grid modes). Then for every cut M and both parity grids:

    ‖B‖ ≤ ‖B‖_HS ≤ ½·√(log 2) = 0.41627…

## Proof (three lines)

1. ‖B‖ ≤ ‖B‖_HS always, and ‖B‖²_HS = Σ_{k≤M<l} 1/(k+l)².
2. The grid points (k,l) (both ≡ 1 mod 2, or both ≡ 0 mod 2) are midpoints of disjoint 2×2
   cells covering ⊂ [24, M+1] × [M+1, ∞). The function (x+y)^{−2} is jointly convex, so the
   2D midpoint value is ≤ the cell average:  Σ 1/(k+l)² ≤ ¼ ∫_{24}^{M+1}∫_{M+1}^{∞} (x+y)^{−2} dy dx.
3. The integral is ∫_{24}^{M+1} dx/(x+M+1) = log((2M+2)/(M+25)) < log 2.  ∎

Continuum origin (why the constant is universal): substituting y ↦ 1/s maps the corner kernel
1/(x+y) on (0,1)×(1,∞) isometrically to the kernel 1/(1+ts) on (0,1)² — a bounded, nearly-rank-one
kernel with HS² = ∫∫ (1+ts)^{−2} = log 2 exactly. Scale invariance of 1/(x+y) makes the bound
cut-independent; the discrete grid contributes the ¼.

## Measured (col window to 120001)

| M | Cauchy corner norm | HS | u-law ‖B‖ | ‖B‖² |
|---|---|---|---|---|
| 259 | 0.3877 | 0.3882 | 0.4740 | 0.2247 |
| 801 | 0.4047 | 0.4053 | 0.4992 | 0.2492 |

Norm ≈ HS to 3 digits (near-rank-one confirmed); both crawl toward the 0.4163 limit from below;
u-law ‖B‖² matches Ein's independent 0.19→0.24 measurement.

## Assembled provable coupling budget

    ‖B_ulaw‖ ≤ |u_c|·0.4163  +  (π/2)·β  +  ‖rank-1 corner‖
             ≤ 0.4163 + 0.9803 + 0.032  =  1.4286        ⇒  ‖B‖² ≤ 2.041

(β-corner priced crude at π/2·β — measured contribution is ~0.1, so this line is the remaining
optimization target; a rigorous improvement route exists via Abel summation with the geometric
bound |Σ_{l grid} e^{iθl}| ≤ 1/|sin θ| ≈ 1.09 at θ = b·log 2, but it is NOT needed to close.)

## Consequence for M* (three-term test, both sectors)

λ_A·(d_C − ‖Off_C‖) ≥ ‖B‖² with λ_A ≥ 2.0 (to be read from the block's own Arb cert; measured
2.03–2.07 both sectors), ‖Off_C‖ ≤ π/2 + (π/2)β + 0.03 = 2.581, d_C(k) ≥ log(kb/2π) − 0.3614 −
0.9/k (closed-form floor, both parities):

    need d_C ≥ 2.581 + 2.041/2.0 = 3.60   ⇒   log(M*b/2π) ≥ 3.965   ⇒   **M* ≈ 117**

[25,117] on one parity grid ≈ 47 modes ≈ 1.1k Arb entries ≈ half an hour laptop — down from
M* ≈ 259 (which assumed an unproven 3.7-class coupling bound) and from ~1160 crude. Every
constant in the chain is now proof-backed: HS integral (above), Hilbert π/2 (classical, odd/even
spacing), β crude (Hilbert twice), rank-1 (Σ 1/k²), d₀ floor (closed forms, C_d = 0.9 pending
only the ε_k bookkeeping).

Caveat registered: the ε_k term in C_d = 0.9 is the one not-yet-proven constant feeding this
arithmetic; measured model-residuals ≤ 0.014 at k ≥ 25 (envelope 0.812 vs working 0.9). If the
final proven C_d lands anywhere ≤ 1.5, M* moves by < 2 modes — insensitive.
