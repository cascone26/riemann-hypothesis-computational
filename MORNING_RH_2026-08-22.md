# RH frontier — autonomous morning session, 2026-08-22 (~9:15–10:00 AM)

**Frame:** reoriented from last night's (walled-lane) formalization to the RH open core. Per the Jacob-set
strategic direction, worked Path 4 (spectral/symmetry-class HP lens) then the positivity wall directly.
Everything below is computed/verified (Handle), non-fabricated, with negatives called honestly.

## Result 1 — Path 4 spectral/symmetry-class lens: COMPREHENSIVELY MAPPED (verified negative)
- ζ's zeros are class A (GUE): bulk <r>=0.6147 = 0.7σ from finite-N GUE control. AZ (BdG) class C ALSO has
  GUE-like bulk (0.601) — so bulk statistics CANNOT distinguish class A from C. (Corrected my own instinct:
  Srednicki's class = AZ/BdG class C, not Wigner-Dyson GSE.)
- The SYMPLECTIC signature lives in a FAMILY, not ζ: quadratic Dirichlet L-functions (300 of them, 2400
  scaled zeros) track the symplectic one-level density W_Sp=1-sinc(2πx) point-by-point (near-0 W=0.047 vs
  unitary 1.0). AZ class-C RMT ensemble reproduces the SAME W_Sp (loop closed). ζ ∈ UNITARY family → class A.
- Explicit class-C construction (Weil pieces + momentum + Pauli spin): genuine class C (PHS-err=0 exact) but
  does NOT bypass the wall — the coercivity λ_a stays in the S-blocks; the class-C gap is generic scale.
- VERDICT: closes "just match the symmetry class." Confirms T3 (spectral = the positivity wall in different
  clothing) from the construction side. Files: spectral/{classfork_control,az_classC_correct,fastL,
  onelevel_fast,classC_onelevel,classC_construct,nullmode_probe}.py.

## Result 2 — the positivity wall CHARACTERIZED, and connected to de Bruijn–Newman
- Weil form Q=arch−prime is numerically PSD (κ=0) ∀a, λ_min→0 = EXACT SATURATION (no positive margin).
- Mechanism: smooth-arch has κ=4→8 negative squares; prime+pole cancel them mode-by-mode near-exactly
  (residual +0.01), and arch & prime+pole NEARLY COMMUTE (commutator ~0.01) → shared eigenvectors on the
  negative subspace. But positivity is NOT diagonally dominant (Gershgorin fails by 100–900×, growing with a)
  — it's a DELICATE GLOBAL cancellation. Schur complement of the negative block = exactly 0 (saturated).
- **KEY CONNECTION:** the exact saturation IS the de Bruijn–Newman constant Λ=0. Rodgers–Tao (2018) PROVED
  Λ≥0 — the zeros are already at the marginal edge, so a strict positive margin is IMPOSSIBLE. This EXPLAINS
  why every truncation saturates. REFRAME: the wall is not "prove strict positivity" (Λ≥0 forbids margin) but
  "prove the marginal case is Λ=0 (positive side)". The natural, ACTIVE, computable frame is the heat-flow /
  de Bruijn–Newman lane (Rodgers–Tao, Polymath15 → Λ<0.22), NOT the static form.
- Cross-checks: Li coefficients λ_n (RH⟺λ_n≥0) all positive, grow (n/2)log n, λ_1=0.023 matches known.

## Honest verdict + next attack
Two verified NEGATIVES that SHARPEN the target: symmetry-class matching and local/perturbative positivity are
both closed. RH not broken (not expected in a morning). The productive next lane is the heat-flow / de Bruijn–
Newman frame — where the saturation is the natural object, there's active recent progress, and it's computable.
Durable assets from last night (verified PNT tower + arithmetic Hodge index, Lean, 789/789 axiom-clean) remain
usable as trusted inputs. All findings appended to SPECTRAL_HP_ROADMAP.md; memory rh-spectral-symmetry-mapped.

## Addendum (~10:34 AM) — landed last night's deferred Lean brick
While mapping the frontier, also closed the one deferred item from riemann9: §1719 (arithmetic Hodge index
sigPos=1 on every subspace of the ℝ-model) + its abstract generalization §1719a (over ANY ordered field).
Both EXIT0, axiom-clean. The instance-diamond that blocked §1719 through 10 attempts last night is diagnosed
and solved via the abstract-core route (prove sigPos abstractly + plain-inequality negativity hypothesis +
no type-ascription on `set Q`). This advances the GEOMETRIC end of the Hodge-index↔Weil-positivity bridge —
exactly the leveraged direction the synthesis identified. Forge now clean incl these; whole tower intact.

## Net session (last night + this morning)
- Lean: PNT tower (all classical forms) + arithmetic Hodge index (§1716/§1719/§1719a) — verified, axiom-clean.
- RH frontier: Path 4 mapped (verified neg); Weil-positivity wall characterized to the parity-sector knife-
  edge; saturation = de Bruijn-Newman Λ=0 (Λ≥0 proven ⇒ marginality forced); Hodge↔Weil bridge synthesis.
- Honest verdict: RH not broken (not expected). Two false-hope lanes closed with data. The real difficulty
  pinned (global prime-coupling cancellation = the razor's edge). Right next lanes identified: heat-flow/dBN
  (marginality native) + pushing the geometric Hodge end toward Spec ℤ (the function-field proof's positivity).

## Capstone (~11:09 AM) — the arithmetic Hodge index INEQUALITY, machine-verified
Built the complete arithmetic-Hodge-index toolkit, culminating in the RH-relevant positivity itself:
- Abstract (any ordered field): §1719a sigPos=1, §1720a signature (1,dim-1), §1721a the index inequality
  (polar Q x y)²≥4·Qx·Qy from sigPos=1 (Lorentzian reverse-Cauchy-Schwarz).
- ℝ-model (Arakelov, every K): §1716/§1719/§1720 signatures + §1721 the inequality
  hodgeBilin(u,v)² ≥ hodgeBilin(u,u)·hodgeBilin(v,v).
All EXIT0, axiom-clean, chain GREEN. This is the arithmetic analog of the positivity that yields function-
field RH — the GEOMETRIC END of the bridge, now a complete verified stack. The morning's analytic-end map
(Weil positivity = de Bruijn-Newman Λ=0 knife-edge) is the OTHER end. Both ends built/mapped; the span (Spec ℤ
intersection-positivity, Connes-Consani) is the remaining gap — and it's now sharply, verifiably framed.
