# §199.3 PLAN — Rigorous certification of p=2-window Weil positivity (the finite theorem)

**Target theorem (computer-assisted):** For all support half-lengths a ∈ [½log2, ½log3], the Weil
quadratic form Q_a on H¹₀(−a,a) is ⪰ 0, unconditionally. This extends the Yoshida/Bombieri
unconditional support threshold past the first prime — the first such extension since 1992/2000.
Then the machine re-runs for p=3, 5, … (the §189 relay), each a new theorem.

Status of inputs: §184–198 (structure fully mapped; 4×4 even block numerically PD, margin ~1.6e-5,
minors down to 3.5e-8); §199.1 Arb screw kernel VALIDATED; acb.integral validated; sympy in flint_venv.
Failed first pass (§199.2 timeout) taught the key lesson: do NOT numerically integrate the composite
kernel — decompose exactly first.

---

## Phase A — EXACT decomposition of the form (the design centerpiece)

Integrate by parts twice (v(±a)=0 kills boundary terms): Q(v)=∫∫g(x−y)v'(x)