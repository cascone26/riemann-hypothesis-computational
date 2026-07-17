#!/usr/bin/env python3
"""
test_real_data.py — Test database with realistic Lean formalization examples.

This script:
1. Loads actual theorem statements from RiemannLean/RobinInequality.lean
2. Populates the database with realistic conjectures
3. Simulates the 5-stage pipeline (formalizability, sketches, attempts, results)
4. Validates query functionality
5. Produces a report showing Phase 1 schema works end-to-end
"""

import sys
import json
import re
from pathlib import Path

# Add bot directory to path
sys.path.insert(0, Path(__file__).parent)
from database import RiemannBotDB


def extract_lean_theorems(lean_file: str) -> list[dict]:
    """Extract theorem statements from a Lean file."""
    theorems = []

    with open(lean_file, 'r') as f:
        content = f.read()

    # Extract theorem/def declarations with names and types
    pattern = r'(?:theorem|def)\s+(\w+)\s*(?:\([^)]*\))?\s*:\s*([^:=]+?)(?:\s*:=\s*by|:=)'
    matches = re.finditer(pattern, content, re.MULTILINE | re.DOTALL)

    for match in matches:
        name = match.group(1)
        statement_raw = match.group(2).strip()
        # Clean up multiline statements
        statement = ' '.join(statement_raw.split())

        theorems.append({
            "name": name,
            "statement": statement,
            "full_statement": f"theorem {name} : {statement}"
        })

    return theorems


def simulate_predictor_scoring(statement: str) -> tuple[float, float, list[str]]:
    """
    Simulate formalizability predictor output.

    Returns: (score, confidence, issues)
    """
    # Heuristics for scoring (in reality, this is a neural network)
    issues = []

    # Check for common issues
    if "sorry" in statement:
        issues.append("contains_sorry")
    if len(statement) > 500:
        issues.append("statement_too_long")
    if "∀" in statement and "∃" in statement:
        issues.append("complex_quantifier_nesting")
    if "ArithmeticFunction" in statement:
        issues.append("needs_specialized_import")

    # Base score: longer/more complex = lower score
    base = 0.7
    complexity_penalty = len(statement) / 1000.0  # Rough heuristic
    issue_penalty = len(issues) * 0.1

    score = max(0.1, min(1.0, base - complexity_penalty - issue_penalty))
    confidence = 0.7 + 0.3 * (1.0 - issue_penalty)  # More issues = lower confidence

    return score, confidence, issues


def generate_proof_sketch(statement: str, strategy: str = "direct") -> str:
    """Generate a simple proof sketch (skeleton with sorry)."""
    theorem_name = statement.split('(')[0].replace('theorem ', '').strip().split()[0]
    body = statement[statement.find(':') + 1:].strip()

    sketch = f"""theorem {theorem_name} {body} := by
  -- Strategy: {strategy}
  -- Step 1: Setup and introduce variables
  intro h₁
  -- Step 2: Apply key lemmas
  sorry  -- TODO: complete the proof
"""
    return sketch


def main():
    # Load RiemannLean theorems (try multiple paths)
    possible_paths = [
        "/Users/scones/projects/riemann/RiemannLean/RiemannLean/RobinInequality.lean",
        Path.home() / "projects/riemann/RiemannLean/RiemannLean/RobinInequality.lean",
    ]

    riemann_lean_path = None
    for path in possible_paths:
        if Path(path).exists():
            riemann_lean_path = str(path)
            break

    if not riemann_lean_path:
        # Fallback: use synthetic theorems
        print("Warning: RiemannLean file not found. Using synthetic theorems for testing.\n")
        theorems = [
            {
                "name": "example1",
                "statement": "∀ n : ℕ, n > 0 → n + 1 > 1",
                "full_statement": "theorem example1 : ∀ n : ℕ, n > 0 → n + 1 > 1"
            },
            {
                "name": "example2",
                "statement": "∀ x y : ℝ, x + y = y + x",
                "full_statement": "theorem example2 : ∀ x y : ℝ, x + y = y + x"
            }
        ]
    else:
        print(f"Loading theorems from {riemann_lean_path}\n")
        theorems = None

    print("=== Phase 1 Database Test: Real Lean Data ===\n")

    # Extract theorems
    if theorems is None:
        theorems = extract_lean_theorems(riemann_lean_path)
    print(f"Found {len(theorems)} theorems\n")

    # Initialize database
    db_path = "/tmp/riemann_bot_test.db"
    db = RiemannBotDB(db_path=db_path)
    print(f"Initialized database: {db_path}\n")

    # Stage 1: Insert conjectures
    print("=== Stage 1: Conjecture Input ===")
    conjecture_ids = []

    for i, theorem in enumerate(theorems):
        conj_id = db.insert_conjecture(
            statement=theorem["full_statement"],
            strategy_type=["dfs", "bfs", "best_first"][i % 3],
            difficulty_estimate=5.0 + (i % 3) * 1.5,
            source="riemann_lean_extracted"
        )
        conjecture_ids.append(conj_id)
        print(f"  [{conj_id}] {theorem['name']}: {theorem['full_statement'][:60]}...")

    print(f"\nInserted {len(conjecture_ids)} conjectures\n")

    # Stage 2: Formalizability scoring
    print("=== Stage 2: Formalizability Predictor ===")
    scored_count = 0

    for conj_id, theorem in zip(conjecture_ids, theorems):
        score, confidence, issues = simulate_predictor_scoring(theorem["full_statement"])

        db.update_formalizability_score(
            conj_id,
            score=score,
            confidence=confidence,
            reason_fail=issues[0] if issues else None,
            issues=issues,
            processing_time_ms=1200 + len(issues) * 100
        )

        status = "✓ PASS" if score >= 0.5 else "✗ FAIL"
        print(f"  [{conj_id}] {status} score={score:.3f} conf={confidence:.2f}")
        if score >= 0.5:
            scored_count += 1

    print(f"\nFormalizability results: {scored_count}/{len(conjecture_ids)} passed filter\n")

    # Stage 3: Proof sketch generation
    print("=== Stage 3: Proof Sketch Generator ===")
    sketch_ids = []

    for conj_id, theorem in zip(conjecture_ids, theorems):
        # Only generate sketches for conjectures that passed filter
        cursor = db.conn.cursor()
        cursor.execute("SELECT formalizability_score FROM conjectures WHERE id = ?", (conj_id,))
        score = cursor.fetchone()[0]

        if score >= 0.5:
            sketch_code = generate_proof_sketch(theorem["full_statement"], strategy="direct")
            sketch_id = db.insert_proof_sketch(
                conj_id,
                sketch_code=sketch_code,
                sketch_strategy="direct",
                generated_from_model="claude_3_5_sonnet"
            )
            sketch_ids.append((conj_id, sketch_id))
            print(f"  [{conj_id}] Generated sketch {sketch_id}")

    print(f"\nGenerated {len(sketch_ids)} proof sketches\n")

    # Stage 4: Simulate formal attempts
    print("=== Stage 4: Formal Verification ===")
    successful_attempts = 0

    for conj_id, sketch_id in sketch_ids:
        cursor = db.conn.cursor()
        cursor.execute("SELECT statement FROM conjectures WHERE id = ?", (conj_id,))
        statement = cursor.fetchone()[0]

        # Simulate attempt
        attempt_id = db.insert_formal_attempt(
            conjecture_id=conj_id,
            lean_code=statement,
            strategy_id=1,  # dfs strategy
            gpu_device="cpu",
            attempt_type="sketch_only",
            sketch_id=sketch_id
        )

        # Simulate outcome (50% success for demo)
        import random
        success = random.random() > 0.5

        if success:
            db.finalize_formal_attempt(
                attempt_id,
                result="success",
                wall_time_ms=2500.0,
                lemma_hits=3,
                lemma_misses=1
            )
            db.update_strategy_result(1, success=True)
            successful_attempts += 1
            print(f"  [{conj_id}] Attempt {attempt_id}: SUCCESS ✓")
        else:
            error_msg = "elaboration elaborator failed"
            db.finalize_formal_attempt(
                attempt_id,
                result="elaboration_error",
                error_message=error_msg,
                error_type="elaboration",
                wall_time_ms=5000.0
            )
            db.update_strategy_result(1, success=False)

            # Record feedback for Claude
            db.insert_feedback(
                conjecture_id=conj_id,
                feedback_type="elaboration_error",
                message=f"Lean elaboration failed: {error_msg}",
                severity="high",
                suggested_fix="Try adding explicit type hints",
                attempt_id=attempt_id
            )
            print(f"  [{conj_id}] Attempt {attempt_id}: FAILED ✗ (elaboration)")

    print(f"\nFormal verification results: {successful_attempts}/{len(sketch_ids)} succeeded\n")

    # Query results
    print("=== Query Results ===")

    # Session summary
    summary = db.get_session_summary()
    print(f"Session Summary:")
    for key, value in summary.items():
        if isinstance(value, float):
            print(f"  {key}: {value:.2f}")
        else:
            print(f"  {key}: {value}")

    # Error summary
    errors = db.get_error_summary()
    print(f"\nError Summary:")
    for error_type, count in errors.items():
        print(f"  {error_type}: {count}")

    # Learning feedback
    feedback = db.get_unused_feedback(limit=3)
    print(f"\nRecent Feedback ({len(feedback)} unused):")
    for fb in feedback:
        print(f"  [{fb['conjecture_id']}] {fb['feedback_type']}: {fb['message'][:60]}...")

    # Strategy info
    print(f"\nStrategy Performance:")
    cursor = db.conn.cursor()
    cursor.execute("SELECT name, success_count, failure_count, win_rate FROM proof_strategies ORDER BY win_rate DESC")
    for row in cursor.fetchall():
        name, wins, fails, wr = row
        print(f"  {name}: {wins}W-{fails}L (rate={wr:.1%})" if wr else f"  {name}: {wins}W-{fails}L")

    db.close()

    print(f"\n✓ Test complete. Database: {db_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
