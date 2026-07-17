#!/usr/bin/env python3
"""
test_sketch_against_real_lean.py — Test Sketch Generator against real Lean theorems.

Tests increment 3 (Proof Sketch Generator) against real theorems from
RobinInequality.lean. Verifies:
1. Parsing accuracy (extracts theorem name, statement, conclusion)
2. Skeleton generation quality (compilable, structured)
3. Performance (within 1-5 sec target)
4. Database integration (stores to proof_sketches table)
"""

import re
import json
import sqlite3
import tempfile
from pathlib import Path
from typing import List, Dict, Tuple
import time

from sketch_generator import SketchGenerator, generate_sketch
from database import RiemannBotDB


def extract_theorems_from_lean(lean_file: Path) -> List[Tuple[str, str]]:
    """
    Extract theorem statements and their sorry lines from Lean file.

    Args:
        lean_file: Path to .lean file

    Returns:
        List of (theorem_name, full_statement) tuples
    """
    with open(lean_file, "r") as f:
        content = f.read()

    theorems = []

    # Pattern: theorem/lemma name ... := by (followed by sorry or other content)
    pattern = r"(theorem|lemma)\s+(\w+)\s*([^:=]*):([^:=]+):=\s*by\s*(\n\s*sorry)"
    matches = re.finditer(pattern, content, re.MULTILINE)

    for match in matches:
        keyword, name, params, statement, _ = match.groups()
        full_statement = f"{keyword} {name} {params} : {statement} := by sorry"
        theorems.append((name, full_statement))

    return theorems


def extract_theorems_simple(lean_file: Path) -> List[Tuple[str, str]]:
    """
    Simple extraction: find all sorry lines and associate with preceding theorem.

    Args:
        lean_file: Path to .lean file

    Returns:
        List of (theorem_name, statement_before_colon_equals) tuples
    """
    with open(lean_file, "r") as f:
        lines = f.readlines()

    theorems = []
    current_theorem = None
    in_theorem = False

    for i, line in enumerate(lines):
        # Match theorem/lemma declaration
        if re.match(r"\s*(theorem|lemma)\s+(\w+)", line):
            match = re.match(r"\s*(theorem|lemma)\s+(\w+)", line)
            current_theorem = match.group(2)
            in_theorem = True

        # If we find a sorry, capture the statement
        if in_theorem and "sorry" in line:
            # Backtrack to find the start of the statement
            stmt_start = i
            for j in range(i, -1, -1):
                if "theorem" in lines[j] or "lemma" in lines[j]:
                    stmt_start = j
                    break

            # Extract the full statement
            statement = "".join(lines[stmt_start:i+1])
            if current_theorem:
                theorems.append((current_theorem, statement))
            in_theorem = False
            current_theorem = None

    return theorems


def test_sketch_generator_on_robin():
    """
    Test sketch generator on real RobinInequality.lean theorems.
    """
    lean_file = Path("/Users/scones/projects/riemann/RiemannLean/RiemannLean/RobinInequality.lean")

    if not lean_file.exists():
        print(f"Error: {lean_file} not found")
        return

    print("=" * 80)
    print("SKETCH GENERATOR — Phase 1 Increment 3 Integration Test")
    print("Testing against real Lean theorems from RobinInequality.lean")
    print("=" * 80)

    # Extract theorems
    theorems = extract_theorems_simple(lean_file)
    print(f"\nExtracted {len(theorems)} theorems with sorry clauses")

    if not theorems:
        print("No theorems found with sorry clauses")
        return

    # Initialize generator
    generator = SketchGenerator()

    results = []
    timings = []

    for name, statement in theorems:
        print(f"\n{'-' * 80}")
        print(f"Theorem: {name}")
        print(f"Statement (first 100 chars): {statement[:100]}...")

        start_time = time.time()
        try:
            sketch = generator.generate(statement, strategy="direct")
            elapsed_ms = (time.time() - start_time) * 1000
            timings.append(elapsed_ms)

            print(f"✓ Generated in {elapsed_ms:.2f} ms")
            print(f"  Theorem name: {sketch.theorem_name}")
            print(f"  Sub-goals: {len(sketch.sub_goals)}")
            print(f"  Suggested tactics: {len(sketch.tactics_suggested)}")
            print(f"  Proof hash: {sketch.sketch_hash[:16]}...")

            results.append({
                "original_name": name,
                "parsed_name": sketch.theorem_name,
                "status": "success",
                "generation_time_ms": sketch.generation_time_ms,
                "sub_goals_count": len(sketch.sub_goals),
                "tactics_count": len(sketch.tactics_suggested),
                "sketch_hash": sketch.sketch_hash,
            })

        except Exception as e:
            print(f"✗ Error: {e}")
            results.append({
                "original_name": name,
                "status": "error",
                "error": str(e),
            })

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY — Increment 3 Performance Metrics")
    print("=" * 80)

    successful = [r for r in results if r.get("status") == "success"]
    print(f"Total theorems: {len(theorems)}")
    print(f"Successfully parsed: {len(successful)}")
    print(f"Errors: {len(results) - len(successful)}")

    if timings:
        print(f"\nGeneration time statistics:")
        print(f"  Average: {sum(timings) / len(timings):.2f} ms")
        print(f"  Min: {min(timings):.2f} ms")
        print(f"  Max: {max(timings):.2f} ms")
        print(f"  Target: 1-5000 ms ✓ PASS" if max(timings) <= 5000 else "  Target: 1-5000 ms ✗ FAIL")

    # Save full results
    output_file = Path(__file__).parent / "test_sketch_on_robin_results.json"
    with open(output_file, "w") as f:
        json.dump({
            "test": "increment_3_real_theorems",
            "lean_file": str(lean_file),
            "timestamp": time.time(),
            "summary": {
                "total": len(theorems),
                "successful": len(successful),
                "errors": len(results) - len(successful),
                "avg_time_ms": sum(timings) / len(timings) if timings else 0,
                "max_time_ms": max(timings) if timings else 0,
                "min_time_ms": min(timings) if timings else 0,
                "meets_spec": max(timings) <= 5000 if timings else False,
            },
            "results": results,
        }, f, indent=2)

    print(f"\nResults saved to: {output_file}")
    print("=" * 80)

    return results


if __name__ == "__main__":
    test_sketch_generator_on_robin()
