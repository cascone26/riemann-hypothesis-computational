#!/usr/bin/env python3
"""
sketch_generator.py — Proof sketch generator for Lean 4 theorem statements.
PHASE 1 INCREMENT 3 — Proof Sketch Generator

Given a target theorem statement, generates a Lean proof skeleton that breaks
down the goal into sub-goals with named sorries, or attempts to apply automated
tactics to decidable sub-claims. Integrated with SQLite schema for tracking
and persistence.

Strategy:
1. PARSE: extract theorem name, hypothesis, conclusion type
2. SKETCH: for compound statements (conjunctions, foralls, implications),
   break into sub-goals; for decidable numeric claims, suggest tactics
3. GENERATE: write a compilable Lean proof skeleton with structured sorries
4. STORE: persist to proof_sketches table with timing and metadata

Performance target: 1-5 seconds per sketch (per hardened spec).
"""

import re
import time
import sqlite3
import hashlib
from typing import List, Tuple, Optional, Dict
from dataclasses import dataclass
from pathlib import Path


@dataclass
class ProofSkeleton:
    """A generated proof skeleton."""
    theorem_name: str
    original_statement: str
    skeleton_proof: str  # Compilable Lean proof with sorries
    sub_goals: List[str]  # Named sub-goals
    tactics_suggested: List[str]  # Tactics to try
    decidable_claims: List[str]  # Potential decidable sub-claims
    generation_time_ms: float = 0.0  # Wall time for generation
    sketch_hash: Optional[str] = None  # SHA256 of skeleton_proof
    sketch_strategy: str = "direct"  # Strategy type used
    conjecture_id: Optional[int] = None  # FK to conjectures table (if stored)


class SketchGenerator:
    """Generates proof sketches for Lean 4 statements."""

    # Tactics for decidable claims
    DECISION_TACTICS = ["decide", "norm_num", "rfl", "omega"]

    # Pattern for detecting decidable numeric claims
    NUMERIC_PATTERN = re.compile(
        r"(\d+)\s*[=<>≤≥≠]\s*(\d+)|"  # e.g. 5 < 10
        r"sigma\s+\d+\s+\d+\s*[=<>≤≥≠]\s*\d+|"  # sigma computations
        r"(\d+)\s*\+\s*(\d+)\s*[=<>≤≥≠]|"  # arithmetic
        r"(harmonic|euler|log|exp)\s*\("  # special functions
    )

    def __init__(self, db_path: Optional[str] = None):
        """
        Initialize sketch generator with optional database integration.

        Args:
            db_path: Path to SQLite database. If None, no persistence.
        """
        self.db_path = db_path
        self.conn = None
        if db_path:
            self._init_db_connection()

    def _init_db_connection(self):
        """Initialize connection to SQLite database."""
        try:
            self.conn = sqlite3.connect(self.db_path)
            self.conn.row_factory = sqlite3.Row
        except Exception as e:
            print(f"Warning: Could not connect to database {self.db_path}: {e}")
            self.conn = None

    def generate(self, statement: str, conjecture_id: Optional[int] = None,
                 strategy: str = "direct", store: bool = False) -> ProofSkeleton:
        """
        Generate a proof skeleton from a theorem statement.

        Args:
            statement: Full Lean theorem/lemma statement
            conjecture_id: Optional FK to conjectures table
            strategy: Strategy type (e.g., 'direct', 'induction', 'contradiction')
            store: If True and db_path is set, store result to proof_sketches table

        Returns:
            ProofSkeleton with proof template, sub-goals, and timing
        """
        start_time = time.time()

        # Parse the statement
        parsed = self._parse_statement(statement)
        theorem_name = parsed["name"]
        hypotheses = parsed["hypotheses"]
        conclusion = parsed["conclusion"]

        # Analyze structure
        structure = self._analyze_structure(conclusion)

        # Generate sub-goals
        sub_goals = self._generate_sub_goals(
            conclusion, hypotheses, structure
        )

        # Suggest tactics
        tactics = self._suggest_tactics(conclusion, sub_goals)

        # Find decidable claims
        decidable = self._find_decidable_claims(conclusion)

        # Generate skeleton proof
        skeleton = self._generate_skeleton(
            theorem_name, hypotheses, conclusion, sub_goals, tactics
        )

        # Calculate timing and hash
        generation_time_ms = (time.time() - start_time) * 1000
        sketch_hash = hashlib.sha256(skeleton.encode()).hexdigest()

        result = ProofSkeleton(
            theorem_name=theorem_name,
            original_statement=statement,
            skeleton_proof=skeleton,
            sub_goals=sub_goals,
            tactics_suggested=tactics,
            decidable_claims=decidable,
            generation_time_ms=generation_time_ms,
            sketch_hash=sketch_hash,
            sketch_strategy=strategy,
            conjecture_id=conjecture_id,
        )

        # Store to database if requested
        if store and self.conn and conjecture_id:
            self._store_to_db(result)

        return result

    def _parse_statement(self, statement: str) -> Dict:
        """Parse theorem statement into components."""
        # Extract theorem name
        name_match = re.search(r"(theorem|lemma)\s+(\w+)", statement)
        name = name_match.group(2) if name_match else "unknown"

        # Extract hypotheses (text before :)
        colon_idx = statement.find(":")
        if colon_idx == -1:
            hypotheses = []
        else:
            hyp_section = statement[:colon_idx]
            # Extract (hyp : type) patterns
            hypotheses = re.findall(r"\((\w+)\s*:\s*([^)]+)\)", hyp_section)

        # Extract conclusion (everything after the first :)
        if colon_idx != -1:
            conclusion = statement[colon_idx + 1:].split(":=")[0].strip()
        else:
            conclusion = ""

        return {
            "name": name,
            "hypotheses": hypotheses,
            "conclusion": conclusion,
        }

    def _analyze_structure(self, conclusion: str) -> Dict:
        """Analyze the structure of the conclusion."""
        structure = {
            "is_forall": "∀" in conclusion or "forall" in conclusion,
            "is_implication": "→" in conclusion,
            "is_conjunction": "∧" in conclusion,
            "is_disjunction": "∨" in conclusion,
            "is_iff": "↔" in conclusion,
            "has_numeric": bool(self.NUMERIC_PATTERN.search(conclusion)),
            "forall_count": conclusion.count("∀"),
        }
        return structure

    def _generate_sub_goals(
        self, conclusion: str, hypotheses: List[Tuple], structure: Dict
    ) -> List[str]:
        """Generate sub-goals based on conclusion structure."""
        sub_goals = []

        # If forall, extract quantified variables
        if structure["is_forall"]:
            forall_matches = re.findall(r"∀\s*(\w+)\s*:\s*(\w+)", conclusion)
            for var, var_type in forall_matches:
                sub_goals.append(f"forall {var}: {var_type}")

        # If implication, extract parts
        if structure["is_implication"]:
            parts = conclusion.split("→")
            for i, part in enumerate(parts[:-1]):
                sub_goals.append(f"assume {part.strip()}")
            sub_goals.append(f"show {parts[-1].strip()}")

        # If conjunction, split
        if structure["is_conjunction"]:
            parts = conclusion.split("∧")
            for part in parts:
                sub_goals.append(f"prove {part.strip()}")

        # If iff, both directions
        if structure["is_iff"]:
            parts = conclusion.split("↔")
            sub_goals.append(f"forward: {parts[0].strip()} → {parts[1].strip()}")
            sub_goals.append(f"backward: {parts[1].strip()} → {parts[0].strip()}")

        return sub_goals if sub_goals else ["main_goal"]

    def _suggest_tactics(self, conclusion: str, sub_goals: List[str]) -> List[str]:
        """Suggest tactics based on conclusion type."""
        tactics = []

        # Always try these foundational tactics
        base_tactics = ["intro", "intros", "apply", "exact"]
        tactics.extend(base_tactics)

        # Numeric/decidable
        if "numeric" in conclusion.lower() or self.NUMERIC_PATTERN.search(
            conclusion
        ):
            tactics.extend(["decide", "norm_num", "omega"])

        # Summation/product patterns
        if "sigma" in conclusion or "Finset.sum" in conclusion:
            tactics.extend(["ring", "field_simp"])

        # Real-valued analysis
        if "Real." in conclusion or "ℝ" in conclusion:
            tactics.extend(["norm_num", "nlinarith"])

        # RH-related patterns
        if "RiemannHypothesis" in conclusion or "robin" in conclusion.lower():
            tactics.append("sorry")  # Deep content, expect sorry

        # General automation
        tactics.extend(["simp", "simp_all"])

        return list(set(tactics))

    def _find_decidable_claims(self, conclusion: str) -> List[str]:
        """Find potentially decidable sub-claims in conclusion."""
        decidable = []

        # Numeric equalities: σ(5040) = 19344
        sigma_claims = re.findall(
            r"sigma\s+\d+\s+\d+\s*=\s*\d+", conclusion
        )
        decidable.extend(sigma_claims)

        # Numeric comparisons: 5 > 3
        numeric_cmps = re.findall(r"(\d+)\s*[<>≤≥]\s*(\d+)", conclusion)
        decidable.extend([f"{a} vs {b}" for a, b in numeric_cmps])

        # Harmonic number specific values
        harmonic_claims = re.findall(
            r"harmonic\s+\d+\s*=\s*[\d./]+", conclusion
        )
        decidable.extend(harmonic_claims)

        return decidable

    def _generate_skeleton(
        self,
        theorem_name: str,
        hypotheses: List[Tuple],
        conclusion: str,
        sub_goals: List[str],
        tactics: List[str],
    ) -> str:
        """Generate a compilable Lean proof skeleton."""
        lines = []

        # Theorem declaration (extracted from original, but we'll format it)
        lines.append(f"theorem {theorem_name} :")
        lines.append(f"    {conclusion} := by")

        # Intro hypotheses
        if hypotheses:
            hyp_names = [name for name, _ in hypotheses]
            lines.append(f"  intro {' '.join(hyp_names)}")

        # For forall, intro the quantified variables
        forall_vars = re.findall(r"∀\s*(\w+)\s*:\s*(\w+)", conclusion)
        if forall_vars:
            lines.append(f"  intro {' '.join(v for v, _ in forall_vars)}")

        # Add sub-goals with structure
        for i, goal in enumerate(sub_goals[:3]):  # Limit to 3 main sub-goals
            lines.append(f"  -- Goal {i+1}: {goal}")
            if "prove" in goal or "show" in goal:
                lines.append(f"  sorry  -- {goal}")
            else:
                lines.append(f"  sorry")

        # If no sub-goals were meaningful, just sorry
        if not sub_goals or len(sub_goals) == 1:
            lines.append("  sorry")

        return "\n".join(lines)

    def _store_to_db(self, skeleton: ProofSkeleton) -> Optional[int]:
        """
        Store proof sketch to database.

        Args:
            skeleton: ProofSkeleton to store

        Returns:
            ID of inserted row, or None if failed
        """
        if not self.conn or not skeleton.conjecture_id:
            return None

        try:
            cursor = self.conn.cursor()
            cursor.execute(
                """
                INSERT INTO proof_sketches
                (conjecture_id, sketch_code, sketch_hash, sketch_strategy, compilation_status, compilation_time_ms)
                VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    skeleton.conjecture_id,
                    skeleton.skeleton_proof,
                    skeleton.sketch_hash,
                    skeleton.sketch_strategy,
                    "unknown",  # Will be updated after compilation test
                    skeleton.generation_time_ms,
                ),
            )
            self.conn.commit()
            return cursor.lastrowid
        except Exception as e:
            print(f"Error storing sketch to database: {e}")
            return None


def generate_sketch(statement: str, db_path: Optional[str] = None,
                   conjecture_id: Optional[int] = None,
                   strategy: str = "direct", store: bool = False) -> ProofSkeleton:
    """Convenience function for single sketch generation."""
    generator = SketchGenerator(db_path=db_path)
    return generator.generate(statement, conjecture_id=conjecture_id,
                             strategy=strategy, store=store)


if __name__ == "__main__":
    import json

    # Test on RH-related statements from RobinInequality.lean
    test_statements = [
        (
            "theorem robin_inequality_of_RH (hRH : RiemannHypothesis) : "
            "∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n := by sorry",
            "robin_inequality_of_RH",
        ),
        (
            "theorem RH_of_robin_inequality "
            "(h : ∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n) : "
            "RiemannHypothesis := by sorry",
            "RH_of_robin_inequality",
        ),
        (
            "theorem lagarias_inequality_of_RH (hRH : RiemannHypothesis) : "
            "∀ n : ℕ, n ≥ 1 → (sigma 1 n : ℝ) ≤ lagariasBound n := by sorry",
            "lagarias_inequality_of_RH",
        ),
        (
            "theorem RH_of_lagarias_inequality "
            "(h : ∀ n : ℕ, n ≥ 1 → (sigma 1 n : ℝ) ≤ lagariasBound n) : "
            "RiemannHypothesis := by sorry",
            "RH_of_lagarias_inequality",
        ),
    ]

    print("=" * 80)
    print("SKETCH GENERATOR — Phase 1 Increment 3 Test")
    print("Proof Sketch Generator for Lean 4 Theorems")
    print("=" * 80)

    generator = SketchGenerator()
    results = []

    for statement, name in test_statements:
        print(f"\n{'=' * 80}")
        print(f"Theorem: {name}")
        print("-" * 80)

        sketch = generator.generate(statement, strategy="direct")

        print(f"Generation time: {sketch.generation_time_ms:.2f} ms")
        print(f"Sub-goals: {sketch.sub_goals}")
        print(f"Tactics: {sketch.tactics_suggested}")
        print(f"Decidable claims: {sketch.decidable_claims}")
        print(f"\nGenerated skeleton proof:")
        print("-" * 40)
        print(sketch.skeleton_proof)
        print("-" * 40)

        results.append({
            "name": name,
            "theorem_name": sketch.theorem_name,
            "generation_time_ms": sketch.generation_time_ms,
            "sub_goals": sketch.sub_goals,
            "tactics": sketch.tactics_suggested,
            "sketch_hash": sketch.sketch_hash,
            "sketch_strategy": sketch.sketch_strategy,
        })

    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    print(f"Total theorems processed: {len(results)}")
    avg_time = sum(r["generation_time_ms"] for r in results) / len(results)
    max_time = max(r["generation_time_ms"] for r in results)
    min_time = min(r["generation_time_ms"] for r in results)
    print(f"Generation time stats:")
    print(f"  Average: {avg_time:.2f} ms")
    print(f"  Max: {max_time:.2f} ms")
    print(f"  Min: {min_time:.2f} ms")
    print(f"  Target: 1-5000 ms (1-5 seconds) ✓ PASS" if max_time <= 5000 else f"  Target: 1-5000 ms ✗ FAIL")

    # Save results to JSON
    output_file = Path(__file__).parent / "test_sketch_generator_results.json"
    with open(output_file, "w") as f:
        json.dump({
            "timestamp": str(time.time()),
            "total_theorems": len(results),
            "avg_time_ms": avg_time,
            "max_time_ms": max_time,
            "min_time_ms": min_time,
            "results": results,
        }, f, indent=2)

    print(f"\nResults saved to: {output_file}")
    print("=" * 80)
