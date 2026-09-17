#!/usr/bin/env bash
# RH Forge AXIOM GATE — the un-foolable correctness check.
# Compiling successfully is NOT enough: a file can compile yet use `sorryAx` or a rogue axiom.
# This loads every forge module and flags any declaration whose transitive axiom set escapes
# the standard trio {propext, Classical.choice, Quot.sound}.  Exit 0 iff the whole forge is clean.
#
# Usage:  ./axiom_gate.sh            # audit every *.lean module in the forge
#         ./axiom_gate.sh Mod1 Mod2  # audit only the named modules
set -u
export PATH="$HOME/.elan/bin:$PATH"
cd "$HOME/rhforge/mathlib4"
LP=$(lake env printenv LEAN_PATH)
cd "$HOME/rhforge/forge"
export LEAN_PATH="$LP:$HOME/rhforge/forge"

if [ "$#" -gt 0 ]; then
  MODS="$*"
else
  # all modules except the auditor itself; basenames without .lean
  MODS=$(for f in *.lean; do b="${f%.lean}"; [ "$b" = axiom_audit ] && continue; [ -f "$b.olean" ] && echo "$b"; done | tr '\n' ' ')
fi

echo "=== RH Forge axiom gate: auditing $(echo $MODS | wc -w) modules ==="
lean --run axiom_audit.lean $MODS
rc=$?
if [ "$rc" -eq 0 ]; then
  echo "GREEN BOARD — every audited declaration uses only {propext, Classical.choice, Quot.sound}."
else
  echo "RED — declarations above escape the standard axiom set (sorryAx / rogue axiom). Fix or quarantine."
fi
exit $rc
