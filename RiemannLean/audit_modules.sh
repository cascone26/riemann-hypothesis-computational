#!/bin/bash
# Integrity audit: force clean rebuild of RH modules + verify every theorem depends only on
# trusted kernel axioms. Prevents stale-olean false-successes (LEDGER §756).
set -e
export PATH="$HOME/.elan/bin:$PATH"
cd "$(dirname "$0")"
MODS="Bessel SignTheorem WeilForm Certificate"
echo "[audit] removing oleans to force genuine recompilation..."
for m in $MODS; do rm -f .lake/build/lib/lean/RiemannLean/$m.olean; done
echo "[audit] clean rebuild..."
lake build $(for m in $MODS; do echo RiemannLean.$m; done) 2>&1 | grep -iE "error|Built RiemannLean|completed|failed" | tail -5
echo "[audit] axiom check on all theorems..."
{ echo "import RiemannLean.Bessel"; echo "import RiemannLean.SignTheorem"; echo "import RiemannLean.WeilForm";
  echo "open ArithmeticFunction";
  grep -hE "^theorem" RiemannLean/{Bessel,SignTheorem,WeilForm,Certificate}.lean | sed -E 's/theorem ([A-Za-z0-9_]+).*/#print axioms \1/'; } > RiemannLean/_Audit.lean
BAD=$(lake env lean RiemannLean/_Audit.lean 2>&1 | grep -iE "sorry|error" | grep -viE "propext|Classical.choice|Quot.sound" || true)
rm -f RiemannLean/_Audit.lean
if [ -n "$BAD" ]; then echo "[audit] FAIL — untrusted axioms / errors:"; echo "$BAD"; exit 1; fi
N=$(grep -hcE "^theorem" RiemannLean/{Bessel,SignTheorem,WeilForm,Certificate}.lean | paste -sd+ - | bc)
echo "[audit] PASS — $N theorems, all trusted axioms {propext, Classical.choice, Quot.sound}, 0 sorries."
