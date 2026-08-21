#!/usr/bin/env bash
# forge-brick <localfile.lean> <ModuleName> ["<ledger one-liner>"]
# ONE command for the whole brick ritual: ship → compile+olean → AXIOM-AUDIT → mirror → bank.
# Refuses to bank unless the module is axiom-clean {propext, Classical.choice, Quot.sound}.
# Kills the ~7-step-per-brick overhead AND makes a dirty brick un-bankable.
set -uo pipefail
LF="${1:?usage: forge-brick <file.lean> <Module> [desc]}"
MOD="${2:?need module name}"
DESC="${3:-}"
[ -f "$LF" ] || { echo "no such file: $LF"; exit 1; }
FORGE='rhforge/forge'
LOCALMIRROR="$HOME/projects/riemann/rhforge/forge"
LEDGER="$HOME/projects/riemann/LEDGER.md"

echo "→ shipping $MOD.lean"
scp -q "$LF" "olympus:$FORGE/$MOD.lean" || { echo "scp failed"; exit 1; }

echo "→ compile + olean + axiom-audit (remote)"
RESULT=$(ssh olympus "export PATH=\$HOME/.elan/bin:\$PATH; cd \$HOME/rhforge/mathlib4 && LP=\$(lake env printenv LEAN_PATH); cd \$HOME/rhforge/forge && export LEAN_PATH=\"\$LP:\$HOME/rhforge/forge\"; \
  if lean -o $MOD.olean $MOD.lean 2>/tmp/fb_$MOD.log; then \
    lean --run axiom_audit.lean $MOD 2>&1 | tail -3; \
  else echo COMPILE_FAIL; grep -iE 'error' /tmp/fb_$MOD.log | head -3; fi")
echo "$RESULT"

if echo "$RESULT" | grep -q 'bad=0'; then
  cp "$LF" "$LOCALMIRROR/$MOD.lean"
  [ -n "$DESC" ] && printf -- '- %s: %s\n' "$MOD" "$DESC" >> "$LEDGER"
  echo "✓ BANKED $MOD — axiom-clean, mirrored, ledgered."
else
  echo "✗ NOT BANKED ($MOD) — not axiom-clean or compile failed. Fix and retry."
  exit 1
fi
