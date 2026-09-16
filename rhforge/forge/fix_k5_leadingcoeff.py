import re, sys

path = sys.argv[1]
with open(path) as f:
    src = f.read()

# Extract each `set quad{i} := BODY with hquad{i}def` block (BODY may span many lines).
bodies = {}
for i in range(1, 6):
    pat = re.compile(
        rf"  set quad{i} :=\n(.*?) with hquad{i}def\n",
        re.DOTALL,
    )
    m = pat.search(src)
    if not m:
        print(f"FAILED to find quad{i} block", file=sys.stderr)
        sys.exit(1)
    bodies[i] = m.group(1)

# Cut out the whole region from the first `set quad1 :=` through the `have hpoly' := ... linear_combination ...` line
start_pat = re.compile(r"  set quad1 :=\n")
start = start_pat.search(src)
end_pat = re.compile(
    r"  have hpoly' : cubic1 \* quad1 \+ cubic2 \* quad2 \+ cubic3 \* quad3 \+ cubic4 \* quad4 \+ cubic5 \* quad5 = 0 := by\n"
    r"    linear_combination hpoly - cubic1 \* hqpe1 - cubic2 \* hqpe2 - cubic3 \* hqpe3 - cubic4 \* hqpe4 - cubic5 \* hqpe5\n"
)
end = end_pat.search(src)
if not (start and end):
    print("FAILED to find region bounds", file=sys.stderr)
    sys.exit(1)

replacement_parts = []
replacement_parts.append(
    "  -- S2623: quad_i via `set` gave `ring1` a transparent let-body it kept unfolding inside the\n"
    "  -- hpoly' linear_combination, timing out even at maxHeartbeats=40M (5x raise, S2622). Fix:\n"
    "  -- CrossK4Key's proven generalize-based opacity trick (revert the hypothesis that mentions the\n"
    "  -- huge body, generalize it to a fresh opaque atom, re-intro) -- true opacity for ring1, unlike `set`.\n"
)
for i in range(1, 6):
    replacement_parts.append(
        f"  revert hqpe{i}\n"
        f"  generalize hQ{i} : ({bodies[i].strip()}) = Q{i}\n"
        f"  intro hqpe{i}\n"
    )
replacement_parts.append(
    "  have hpoly' : cubic1 * Q1 + cubic2 * Q2 + cubic3 * Q3 + cubic4 * Q4 + cubic5 * Q5 = 0 := by\n"
    "    linear_combination hpoly - cubic1 * hqpe1 - cubic2 * hqpe2 - cubic3 * hqpe3 - cubic4 * hqpe4 - cubic5 * hqpe5\n"
)
replacement = "".join(replacement_parts)

new_src = src[: start.start()] + replacement + src[end.end() :]

# The coefficient-extraction rw currently unfolds hquad{i}def; swap to `<- hQ{i}` (Q_i = BODY_i, so
# `<- hQi` rewrites Qi back to its raw polynomial body for the coeff/simp pipeline).
old_rw = "rw [hcubic1def, hcubic2def, hcubic3def, hcubic4def, hcubic5def, hquad1def, hquad2def, hquad3def, hquad4def, hquad5def] at hc19"
new_rw = "rw [hcubic1def, hcubic2def, hcubic3def, hcubic4def, hcubic5def, ← hQ1, ← hQ2, ← hQ3, ← hQ4, ← hQ5] at hc19"
if old_rw not in new_src:
    print("FAILED to find coeff-extraction rw line", file=sys.stderr)
    sys.exit(1)
new_src = new_src.replace(old_rw, new_rw)

with open(path, "w") as f:
    f.write(new_src)
print("OK, rewrote", path)
