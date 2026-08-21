import Lean
open Lean

/-- RH Forge axiom-audit gate.  For each declaration DEFINED in the given modules, compute its
    transitive axiom set and flag any that escapes the standard trio {propext, Classical.choice,
    Quot.sound} (or that uses `sorryAx`).  Usage:
      lake env lean --run axiom_audit.lean Module1 Module2 ...
    Exit code 0 iff every checked declaration is clean. -/

def allowedAxioms : List Name := [`propext, `Classical.choice, `Quot.sound]

unsafe def main (args : List String) : IO UInt32 := do
  initSearchPath (← findSysroot)
  let modNames := args.toArray.map String.toName
  let env ← importModules (modNames.map (fun m => { module := m })) Options.empty
  let ctx : Core.Context := { fileName := "<audit>", fileMap := default }
  let coreState : Core.State := { env := env }
  let mut checked := 0
  let mut bad := 0
  for modName in modNames do
    let some modIdx := env.getModuleIdx? modName
      | do IO.println s!"WARN  module not found: {modName}"; continue
    for (declName, _) in env.constants.toList do
      if env.getModuleIdxFor? declName != some modIdx then continue
      if declName.isInternal || declName.isImplementationDetail then continue
      let eio := (Lean.collectAxioms declName : CoreM (Array Name)).run ctx coreState
      let (axs, _) ← eio.toIO (fun _ => IO.userError "core error")
      checked := checked + 1
      let extra := axs.filter (fun a => !allowedAxioms.contains a)
      unless extra.isEmpty do
        IO.println s!"BAD   {declName}  uses  {extra.toList}"
        bad := bad + 1
  IO.println s!"=== AXIOM AUDIT: checked={checked} clean={checked - bad} bad={bad} ==="
  return (if bad == 0 then 0 else 1)
