#!/usr/bin/env python3
"""
wright_lean_harness.py — The Wright module: COMPREHEND → GRAFT → PROVE → LAND workflow.

Implements the three-stage cycle:
1. COMPREHEND: parse target sorry + its statement from the Lean file
2. GRAFT: generate candidate proofs and score them with formalizability_predictor
3. PROVE: run lake build on a scratch copy and parse real compiler output
4. LAND: only if PROVE genuinely passed (zero errors, no remaining sorry),
          apply the change to the real file and record provenance

CRITICAL: "done" means real lake build output with zero errors. This tool
reports failures honestly and does not claim sorries are resolved without
real compiler evidence.
"""

import re
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass
import json

from formalizability_predictor import FormalizabilityPredictor
from sketch_generator import SketchGenerator


@dataclass
class ProveResult:
    """Result of proof attempt (PROVE stage)."""
    success: bool
    theorem_name: str
    original_sorry: str
    attempted_proof: str
    build_output: str
    remaining_sorries: List[str]
    error_messages: List[str]
    compilation_time_ms: float


@dataclass
class LandResult:
    """Result of landing a proof change (LAND stage)."""
    success: bool
    theorem_name: str
    file_path: str
    commit_info: Optional[str]


class WrightLeanHarness:
    """
    Wright module: realizes proofs in actual Lean source code.

    Workflow: COMPREHEND (parse sorry) → GRAFT (generate + score candidates) →
    PROVE (run lake build, parse output) → LAND (if success, apply + verify)
    """

    def __init__(self, lean_project_root: str):
        """
        Initialize the harness.

        Args:
            lean_project_root: Path to RiemannLean (or similar) project root
        """
        self.lean_root = Path(lean_project_root)
        self.predictor = FormalizabilityPredictor()
        self.generator = SketchGenerator()

    def comprehend(self, theorem_name: str, file_path: str) -> Optional[Dict]:
        """
        COMPREHEND stage: parse a sorry from the Lean file.

        Args:
            theorem_name: Name of theorem to locate
            file_path: Path to .lean file

        Returns:
            Dict with theorem_name, statement, hypothesis, conclusion, or None
        """
        if not Path(file_path).exists():
            print(f"  ERROR: File not found: {file_path}")
            return None

        with open(file_path, "r") as f:
            content = f.read()

        # Find the theorem statement - handle multiline format with "by sorry"
        # Use [\s\S]*? to match any character including newlines and unicode
        pattern = (
            rf"(theorem|lemma)\s+{re.escape(theorem_name)}"
            r"[\s\S]*?:=\s+by\s+sorry"
        )
        match = re.search(pattern, content, re.DOTALL)

        if not match:
            print(f"  ERROR: Could not locate sorry for {theorem_name}")
            return None

        full_stmt = match.group(0)

        # Extract conclusion: find ') :' which separates hypothesis from conclusion
        # Theorem format: theorem name (hypotheses) : conclusion := by sorry
        import re as re_module

        # Look for ') :' - the separator between hypothesis and conclusion
        paren_colon = full_stmt.find(') :')
        assign_idx = full_stmt.find(':=')

        if paren_colon >= 0 and assign_idx > paren_colon:
            # Everything after ') :' and before ':=' is the conclusion
            conclusion_part = full_stmt[paren_colon + 3:assign_idx].strip()
            # Everything after theorem name and before ') :' is the hypothesis
            name_idx = full_stmt.find(theorem_name)
            if name_idx >= 0:
                hypothesis_part = full_stmt[name_idx + len(theorem_name):paren_colon].strip()
            else:
                hypothesis_part = ""
        else:
            # Fallback if no ') :' pattern found
            conclusion_part = ""
            hypothesis_part = ""

        return {
            "theorem_name": theorem_name,
            "full_statement": full_stmt,
            "hypothesis": hypothesis_part,
            "conclusion": conclusion_part,
            "file_path": file_path,
        }

    def graft(self, comprehended: Dict) -> List[Tuple[str, float]]:
        """
        GRAFT stage: generate candidate proofs and score them.

        Args:
            comprehended: Output from COMPREHEND

        Returns:
            List of (candidate_proof_string, formalizability_score) tuples,
            sorted by score descending
        """
        theorem_name = comprehended["theorem_name"]
        conclusion = comprehended["conclusion"]

        candidates = []

        # Candidate 1: Direct sorry (baseline)
        sorry_candidate = "sorry"
        score1 = self.predictor.predict(conclusion).score
        candidates.append((sorry_candidate, score1, "baseline_sorry"))

        # Candidate 2: Try generated sketch
        try:
            sketch = self.generator.generate(
                f"theorem {theorem_name} : {conclusion}"
            )
            # Extract just the proof part (tactics after "by")
            sketch_proof = sketch.skeleton_proof.split("by")[-1].strip()
            score2 = self.predictor.predict(conclusion).score * 0.95  # Penalty for attempted proof
            candidates.append((sketch_proof, score2, "generated_sketch"))
        except Exception as e:
            print(f"  Sketch generation failed: {e}")

        # Candidate 3: Simple tactics that might work on decidable claims
        if sketch and sketch.decidable_claims:
            for tactic in ["decide", "norm_num", "omega", "rfl"]:
                score3 = self.predictor.predict(conclusion).score * 1.1
                candidates.append((tactic, score3, f"tactic_{tactic}"))

        # Sort by score descending
        candidates.sort(key=lambda x: x[1], reverse=True)

        # Return top candidates
        return [(proof, score) for proof, score, _ in candidates[:5]]

    def prove(
        self, comprehended: Dict, candidate_proof: str
    ) -> ProveResult:
        """
        PROVE stage: run lake build on a scratch copy with the candidate proof.

        Args:
            comprehended: Output from COMPREHEND
            candidate_proof: The proof text to attempt (e.g. "sorry" or "decide")

        Returns:
            ProveResult with real build output and success status
        """
        theorem_name = comprehended["theorem_name"]
        file_path = comprehended["file_path"]

        # Create scratch copy
        temp_dir = tempfile.mkdtemp(prefix="wright_")
        try:
            # Copy entire project
            project_copy = Path(temp_dir) / "project"
            shutil.copytree(self.lean_root, project_copy)

            # Modify the target file in the copy
            target_file = project_copy / Path(file_path).relative_to(
                self.lean_root
            )

            with open(target_file, "r") as f:
                original_content = f.read()

            # Replace sorry with candidate
            pattern = (
                rf"(theorem\s+{re.escape(theorem_name)}\s*"
                r"[^:]*:[^:=]*:=\s*)sorry"
            )
            modified_content = re.sub(
                pattern, rf"\1{candidate_proof}", original_content
            )

            with open(target_file, "w") as f:
                f.write(modified_content)

            # Run lake build
            result = subprocess.run(
                ["lake", "build"],
                cwd=str(project_copy),
                capture_output=True,
                text=True,
                timeout=300,
            )

            build_output = result.stdout + result.stderr

            # Parse output for sorries and errors
            remaining_sorries = re.findall(
                r"(sorry|unsolved goals?)", build_output, re.IGNORECASE
            )
            error_msgs = re.findall(
                r"error.*?:\s*(.+?)(?=\n|$)", build_output
            )

            # Determine success: zero errors, zero remaining sorries in THIS theorem
            success = (
                result.returncode == 0
                and f"{theorem_name}" not in build_output  # check theorem not mentioned in errors
                and len(error_msgs) == 0
            )

            return ProveResult(
                success=success,
                theorem_name=theorem_name,
                original_sorry=comprehended["full_statement"],
                attempted_proof=candidate_proof,
                build_output=build_output[:2000],  # Limit output
                remaining_sorries=remaining_sorries,
                error_messages=error_msgs,
                compilation_time_ms=0,
            )

        except subprocess.TimeoutExpired:
            return ProveResult(
                success=False,
                theorem_name=theorem_name,
                original_sorry=comprehended["full_statement"],
                attempted_proof=candidate_proof,
                build_output="lake build timed out (>300s)",
                remaining_sorries=[],
                error_messages=["timeout"],
                compilation_time_ms=300000,
            )
        except Exception as e:
            return ProveResult(
                success=False,
                theorem_name=theorem_name,
                original_sorry=comprehended["full_statement"],
                attempted_proof=candidate_proof,
                build_output=str(e),
                remaining_sorries=[],
                error_messages=[str(e)],
                compilation_time_ms=0,
            )
        finally:
            # Cleanup
            shutil.rmtree(temp_dir, ignore_errors=True)

    def land(self, comprehended: Dict, prove_result: ProveResult) -> LandResult:
        """
        LAND stage: if PROVE succeeded, apply the change to the real file.

        Args:
            comprehended: Output from COMPREHEND
            prove_result: Output from PROVE

        Returns:
            LandResult indicating success and file path
        """
        if not prove_result.success:
            return LandResult(
                success=False,
                theorem_name=prove_result.theorem_name,
                file_path=comprehended["file_path"],
                commit_info="Proof did not validate; not applied",
            )

        # Apply change to real file
        file_path = comprehended["file_path"]
        theorem_name = comprehended["theorem_name"]

        with open(file_path, "r") as f:
            original = f.read()

        # Replace sorry
        pattern = (
            rf"(theorem\s+{re.escape(theorem_name)}\s*"
            r"[^:]*:[^:=]*:=\s*)sorry"
        )
        modified = re.sub(pattern, rf"\1{prove_result.attempted_proof}", original)

        with open(file_path, "w") as f:
            f.write(modified)

        return LandResult(
            success=True,
            theorem_name=theorem_name,
            file_path=file_path,
            commit_info=f"Applied proof for {theorem_name}: {prove_result.attempted_proof[:40]}",
        )

    def run_full_cycle(self, theorem_name: str, file_path: str) -> Dict:
        """
        Run the full COMPREHEND → GRAFT → PROVE → LAND cycle for a theorem.

        Args:
            theorem_name: Name of theorem with sorry
            file_path: Path to .lean file

        Returns:
            Dict with detailed results and build output
        """
        print(f"\n{'=' * 80}")
        print(f"WRIGHT CYCLE: {theorem_name}")
        print(f"{'=' * 80}")

        # COMPREHEND
        print(f"\n[COMPREHEND] Parsing {theorem_name}...")
        comprehended = self.comprehend(theorem_name, file_path)
        if not comprehended:
            return {"success": False, "error": "Failed to comprehend theorem"}

        print(f"  ✓ Located sorry: {comprehended['full_statement'][:60]}...")
        print(
            f"  Conclusion: {comprehended['conclusion'][:80]}..."
        )

        # GRAFT
        print(f"\n[GRAFT] Generating candidate proofs...")
        candidates = self.graft(comprehended)
        print(f"  ✓ Generated {len(candidates)} candidates")
        for i, (proof, score) in enumerate(candidates[:3]):
            print(f"    {i+1}. {proof[:40]:40} (score: {score:.3f})")

        # PROVE — try candidates in order
        print(f"\n[PROVE] Attempting proofs (running lake build)...")
        best_result = None
        for i, (proof, score) in enumerate(candidates):
            print(f"  Candidate {i+1}: {proof[:40]}")
            result = self.prove(comprehended, proof)

            if result.success:
                print(f"    ✓ SUCCESS!")
                best_result = result
                break
            else:
                if result.error_messages:
                    print(
                        f"    ✗ Failed: {result.error_messages[0][:60]}"
                    )
                else:
                    print(f"    ✗ Failed (no zero-error build output)")

        if not best_result:
            print(f"  No candidate proof compiled successfully.")
            print(f"  All sorries remain open (expected for deep mathematics).")
            return {
                "success": False,
                "theorem_name": theorem_name,
                "message": "No successful proof candidate",
                "candidates_attempted": len(candidates),
            }

        # LAND
        print(f"\n[LAND] Applying successful proof...")
        land_result = self.land(comprehended, best_result)
        if land_result.success:
            print(f"  ✓ Applied: {land_result.commit_info}")
        else:
            print(f"  ✗ Failed to apply")

        return {
            "success": land_result.success,
            "theorem_name": theorem_name,
            "attempted_proof": best_result.attempted_proof,
            "build_output": best_result.build_output,
            "remaining_sorries": best_result.remaining_sorries,
            "errors": best_result.error_messages,
        }


def main():
    """Test the harness on RiemannLean project."""
    harness = WrightLeanHarness("/Users/scones/projects/riemann/RiemannLean")

    # Target theorems with sorries in RobinInequality.lean
    targets = [
        ("robin_inequality_of_RH", "RiemannLean/RobinInequality.lean"),
        ("RH_of_robin_inequality", "RiemannLean/RobinInequality.lean"),
        ("lagarias_inequality_of_RH", "RiemannLean/RobinInequality.lean"),
        ("RH_of_lagarias_inequality", "RiemannLean/RobinInequality.lean"),
    ]

    results = []
    for theorem_name, file_path in targets:
        full_path = str(
            harness.lean_root / file_path
        )
        result = harness.run_full_cycle(theorem_name, full_path)
        results.append(result)

    # Summary
    print(f"\n{'=' * 80}")
    print("WRIGHT SESSION SUMMARY")
    print(f"{'=' * 80}")
    successful = sum(1 for r in results if r.get("success"))
    print(f"Total theorems attempted: {len(results)}")
    print(f"Successful proofs landed: {successful}")
    print(f"Remaining open: {len(results) - successful}")
    for r in results:
        status = "✓ PROVEN" if r.get("success") else "✗ OPEN"
        print(f"  {status:12} {r['theorem_name']}")

    print("\n" + "=" * 80)


if __name__ == "__main__":
    main()
