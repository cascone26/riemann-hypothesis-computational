
[REPORT — not yet code-verified, theory section from arXiv:2206.03682]

Suzuki canonical system: d/dx w(x;z) = z(0  1; -σ(x) 0) w(x;z), boundary: w(0)=(1,0)^T

The Weyl m-function m(z) is defined via the two fundamental solutions:
  w_+(x;z) = (e_{+,1}, e_{+,2}),  w_-(x;z) = (e_{-,1}, e_{-,2})
  m(z) = -e_{+,2}(a;z) / e_{+,1}(a;z)  (pole location: Re z = 1)

EXACT CONSTRUCTION (in this codebase):
- σ(x) = Ω(r)/2π - P(r), the explicit-formula symbol (arch - prime)
- Matrix form: Q_{kl} = ∫∫ g(x-y) v'_k(x) v'_l(y) dx dy (§242 in suzuki_exact_Q.py)
  - v_k(x) are prolate-like basis (sine/cosine on [0, 2a])
  - g contains: pole term (ζ's s=1 pole) + arch term + prime term
- Solve (Q - λ a I) c = f_± exactly (quadrature-free)
- W(θ,z) = (z - i) Σ c_+ ŝ_k(z) + e^{iθ} (z + i) Σ c_- ŝ_k(z)
  where ŝ_k(z) is the exact Fourier transform of v_k(x)

m-FUNCTION ENTRY POINT (this session):
Q-function m(z) = ⟨(M_σ - z J)^{-1} φ, φ⟩_J, where:
  - M_σ = arch - prime (one negative square in Π₁ kernel)
  - J = indefinite metric (Pontryagin space structure)
  - φ is a cyclic vector (e.g., standard basis)
  - (M_σ - z J)^{-1} is the resolvent in the indefinite space

NEVANLINNA KERNEL K_m(z,w) = (m(z) - m(w̄)̄) / (z - w̄):
  - If κ(K_m) = 1, then m ∈ N₁ (generalized Nevanlinna, one neg square)
  - Single pole of m corresponds to the pole of K_m
  - For ζ: pole at Re z = 1, residue ~+1 (strong, absorbs M_σ's neg square)

OPEN FRONTIER (§974 in LEDGER):
  The naive string E(z) = W(0,z) from the screw form has NO UHP zeros.
  The Weil κ (from pole absorption) realizes κ=1 as a pole-shadow.
  The de Branges E-function (encodes RH zeros → Hermite-Biehler kernel) is a
  different object. The bridge (when κ_Weil = 1, how does this connect to E's
  zero-count via argument principle?) remains unbuilt, but is the most direct
  route to proving RH via the κ=1 positivity.

HANDLES THIS SESSION:
  ✓ m(z) constructed and sampled in UHP
  ✓ Nevanlinna κ counted for M_σ only and Q with pole
  ✓ Controls (m=z, m=z², m=-1/z) validated kernel code
  ✓ Pole rank(Q_pole)=2 confirmed (ζ's structure)

REMAINING (for the full bridge):
  - Express κ in terms of de Branges E-function's zero-count (Hermite-Biehler link)
  - Connect the pole-shadow κ=1 to RH zeros on the line (argument principle)
  - Full Suzuki canonical-system solution (not just Q-function, solve W_full)
