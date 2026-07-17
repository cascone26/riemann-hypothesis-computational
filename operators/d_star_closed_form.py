#!/usr/bin/env python3
"""
Plan #3: d* = 2.46310724598894... — closed-form hypothesis engine.

Test ~2000 candidate expressions built from known constants.
A closed form would reveal WHY the RG flow has this fixed-point dimension
and potentially give a road to proving IR stability analytically.
"""
import math, json, itertools

D_STAR  = 2.46310724598894   # from rg_fixedpoint.json (most precise)
D_STAR2 = 2.4630441410268573  # from heat_kernel_analysis.json
ALPHA   = 1.23155362299447    # heat-kernel decay exponent (d* = 2*alpha)
TOL = 1e-6   # match threshold (6 decimal places)
TIGHT = 1e-9 # tight match

GAMMA1  = 14.134725141734693   # first Riemann zero imaginary part
GAMMA2  = 21.022039638771555
GAMMA3  = 25.010857580145688
EULER_GAMMA = 0.5772156649015329  # Euler-Mascheroni γ
PI    = math.pi
E_    = math.e
PHI   = (1 + math.sqrt(5)) / 2   # golden ratio
LOG2  = math.log(2)
LOG3  = math.log(3)
LOGPI = math.log(PI)
LOG2PI = math.log(2*PI)
CATALAN = 0.9159655941772190   # Catalan's constant
APERY   = 1.2020569031595943   # Apéry's constant ζ(3)
ZETA2   = PI**2 / 6            # ζ(2)
ZETA4   = PI**4 / 90           # ζ(4)

# ══════════════════════════════════════════════════════════════
# Level 0: single constants and simple transforms
# ══════════════════════════════════════════════════════════════

base_consts = {
    "pi":          PI,
    "e":           E_,
    "phi":         PHI,
    "gamma_EM":    EULER_GAMMA,
    "log2":        LOG2,
    "log3":        LOG3,
    "logpi":       LOGPI,
    "log2pi":      LOG2PI,
    "catalan":     CATALAN,
    "apery":       APERY,
    "zeta2":       ZETA2,
    "zeta4":       ZETA4,
    "gamma1":      GAMMA1,
    "gamma2":      GAMMA2,
    "sqrt2":       math.sqrt(2),
    "sqrt3":       math.sqrt(3),
    "sqrt5":       math.sqrt(5),
    "pi2":         PI**2,
    "pi/2":        PI/2,
    "pi/3":        PI/3,
    "pi/4":        PI/4,
    "2pi":         2*PI,
    "4pi":         4*PI,
    "e2":          E_**2,
    "1/e":         1/E_,
    "1/pi":        1/PI,
    "loglog2":     math.log(LOG2) if LOG2 > 0 else None,
}

matches = []

def check(val, expr):
    if val is None: return
    try:
        v = float(val)
        if not math.isfinite(v): return
        for target, tname in [(D_STAR, "d*"), (ALPHA, "alpha=d*/2")]:
            if abs(v - target) < TOL:
                quality = "TIGHT" if abs(v - target) < TIGHT else "close"
                matches.append({
                    "expr": expr, "val": float(v),
                    "target": tname, "target_val": target,
                    "error": float(abs(v - target)),
                    "quality": quality
                })
                tag = "***TIGHT***" if quality == "TIGHT" else "   match   "
                print(f"  {tag}  {expr} = {v:.10f}  (err={abs(v-target):.2e}  vs {tname}={target})")
    except: pass

print(f"d* = {D_STAR}")
print(f"alpha = {ALPHA}  (alpha = d*/2)")
print(f"Searching for closed form... (tolerance={TOL})")
print()

# ── Level 0: single constants ──
print("── Level 0: single constants ──")
for name, val in base_consts.items():
    if val: check(val, name)

# ── Level 1: simple arithmetic combos ──
print("── Level 1: simple combos (a+b, a*b, a/b, a^b for small a,b) ──")
nums = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12]
for n in nums:
    for name, c in base_consts.items():
        if c is None: continue
        check(c + n,   f"{name}+{n}")
        check(c - n,   f"{name}-{n}")
        check(n - c,   f"{n}-{name}")
        check(c * n,   f"{n}*{name}")
        check(c / n,   f"{name}/{n}")
        check(n / c,   f"{n}/{name}")
        try: check(c**n,    f"{name}^{n}")
        except: pass
        try: check(n**c,    f"{n}^{name}")
        except: pass
        try: check(math.exp(c/n), f"exp({name}/{n})")
        except: pass
        try: check(math.exp(n/c), f"exp({n}/{name})")
        except: pass
        try: check(math.log(c*n) if c*n>0 else None, f"log({n}*{name})")
        except: pass
        try: check(math.sqrt(c*n) if c*n>=0 else None, f"sqrt({n}*{name})")
        except: pass

# ── Level 2: ratios of pairs ──
print("── Level 2: ratios and products of constant pairs ──")
clist = list(base_consts.items())
for (na, va), (nb, vb) in itertools.combinations(clist, 2):
    if va is None or vb is None: continue
    check(va + vb, f"{na}+{nb}")
    check(va - vb, f"{na}-{nb}")
    check(vb - va, f"{nb}-{na}")
    check(va * vb, f"{na}*{nb}")
    if abs(vb) > 1e-12: check(va / vb, f"{na}/{nb}")
    if abs(va) > 1e-12: check(vb / va, f"{nb}/{na}")
    try: check(va**vb, f"{na}^{nb}")
    except: pass
    try: check(vb**va, f"{nb}^{na}")
    except: pass
    try: check(math.sqrt(va*vb) if va*vb>=0 else None, f"sqrt({na}*{nb})")
    except: pass

# ── Level 3: gamma1-specific ──
print("── Level 3: gamma1-specific combos ──")
g1 = GAMMA1
for n in [1,2,3,4,5,6]:
    for m in [1,2,3,4,5,6,7,8,9,10,12,14,16]:
        check(n*PI/g1 + m,   f"{n}*pi/{m}+... " if False else f"{n}*pi/gamma1+{m}")
        check(n*PI/g1,       f"{n}*pi/gamma1")
        check(g1/(n*PI) + m, f"gamma1/({n}*pi)+{m}")
        check(n/g1,          f"{n}/gamma1")
        check(g1/n,          f"gamma1/{n}")
        check(math.exp(n/g1) if n < 30 else None, f"exp({n}/gamma1)")
        check(math.log(g1/n) if g1/n > 0 else None, f"log(gamma1/{n})")
        check(math.sqrt(g1/n), f"sqrt(gamma1/{n})")

# ── Level 4: zeta-value combos ──
print("── Level 4: zeta-value combos ──")
# ζ(s) for small s
for s in [2,3,4,5,6]:
    zs = math.pi**(2*s) / (2**(2*s-1) * abs(math.factorial(s) if s < 20 else 1)) if s%2==0 else None
for combo in [
    ("2+pi/gamma1",   2 + PI/g1),
    ("2+2/gamma1",    2 + 2/g1),
    ("pi-gamma1/10",  PI - g1/10),
    ("2*gamma_EM+pi/2", 2*EULER_GAMMA + PI/2),
    ("log(gamma1)-1", math.log(g1)-1),
    ("log(2*gamma1)", math.log(2*g1)),
    ("e^(1/gamma_EM)", E_**(1/EULER_GAMMA)),
    ("pi^gamma_EM",   PI**EULER_GAMMA),
    ("e^gamma_EM",    E_**EULER_GAMMA),
    ("pi/e+1",        PI/E_+1),
    ("(pi+e)/pi",     (PI+E_)/PI),
    ("pi^2/4",        PI**2/4),
    ("log(pi^2)",     2*LOGPI),
    ("sqrt(pi+e)",    math.sqrt(PI+E_)),
    ("gamma1/2pi*e",  g1/(2*PI)*E_),
    ("2*log(2*pi)",   2*LOG2PI),
    ("3/2*log(pi)",   1.5*LOGPI),
    ("pi/log(pi)",    PI/LOGPI),
    ("e/log(e)",      E_/1.0),   # log(e)=1
    ("2+2*gamma_EM",  2+2*EULER_GAMMA),
    ("pi-e/pi",       PI - E_/PI),
    ("pi^2/6-1",      ZETA2-1),
    ("1+pi/2/log(2)", 1+PI/2/LOG2),
    ("gamma1/(2*pi)+2", g1/(2*PI)+2),
]:
    check(combo[1], combo[0])

# ── Level 5: targeted —— alpha = d*/2 = 1.2315... ──
print("── Level 5: alpha=1.2315... specific ──")
# 1/alpha = 0.8119...?
inv_alpha = 1/ALPHA
print(f"  1/alpha = {inv_alpha:.10f}")
check(inv_alpha, "1/alpha")
# Is alpha = pi/e^gamma?
check(PI/math.exp(EULER_GAMMA), "pi/exp(gamma_EM)")
check(PI/(2*E_**EULER_GAMMA), "pi/(2*exp(gamma_EM))")
check(math.sqrt(PI/E_), "sqrt(pi/e)")
check(math.log(PI), "log(pi)")
check(math.log(3.42), "log(3.42)")
check(PI**0.75/E_**0.25, "pi^(3/4)/e^(1/4)")

print()
print(f"Total matches found: {len(matches)}")
print()
print("=== ALL MATCHES (sorted by error) ===")
matches.sort(key=lambda x: x["error"])
for m in matches[:30]:
    print(f"  [{m['quality']:5s}] {m['expr']:40s} = {m['val']:.12f}  err={m['error']:.3e}  vs {m['target']}")

with open("results/d_star_closed_form.json","w") as fh:
    json.dump({"d_star": D_STAR, "alpha": ALPHA, "matches": matches}, fh, indent=2)
print("\nSaved → results/d_star_closed_form.json")
