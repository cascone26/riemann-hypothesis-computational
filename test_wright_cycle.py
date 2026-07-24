#!/usr/bin/env python3
"""
Test script for the Wright harness on RiemannLean RobinInequality theorems.

This script runs the full COMPREHEND → GRAFT → PROVE → LAND cycle on each
sorry in RobinInequality.lean and captures real lake build output for the report.
"""

import sys
import json
from pathlib import Path

# Add bot directory to path
sys.path.insert(0, str(Path(__file__).parent / "bot"))

from wright_lean_harness import WrightLeanHarness


def main():
    """Run the Wright cycle on RiemannLean theorems."""

    print("\n" + "=" * 100)
    print("WRIGHT LEAN HARNESS — Testing RiemannLean.RobinInequality")
    print("=" * 100)

    # Initialize harness
    harness = WrightLeanHarness("/Users/scones/projects/riemann/RiemannLean")

    # Theorems with sorries in RobinInequality.lean
    targets = [
        ("robin_inequality_of_RH", "RiemannLean/RobinInequality.lean"),
        ("RH_of_robin_inequality", "RiemannLean/RobinInequality.lean"),
        ("lagarias_inequality_of_RH", "RiemannLean/RobinInequality.lean"),
        ("RH_of_lagarias_inequality", "RiemannLean/RobinInequality.lean"),
    ]

    # Run cycles
    all_results = []
    for theorem_name, file_path_rel in targets:
        full_path = str(harness.lean_root / file_path_rel)
        result = harness.run_full_cycle(theorem_name, full_path)
        all_results.append(result)

    # Print summary
    print("\n" + "=" * 100)
    print("WRIGHT SESSION SUMMARY")
    print("=" * 100)

    successful = sum(1 for r in all_results if r.get("success"))
    open_count = len(all_results) - successful

    print(f"\nTheorems attempted: {len(all_results)}")
    print(f"Successful proofs landed: {successful}")
    print(f"Remaining open: {open_count}")
    print(f"\nDetailed results:")

    for r in all_results:
        status = "✓ PROVEN" if r.get("success") else "✗ OPEN"
        print(f"\n  {status:12} {r['theorem_name']}")
        if not r.get("success"):
            if r.get("message"):
                print(f"    Message: {r['message']}")
            if r.get("candidates_attempted"):
                print(f"    Candidates attempted: {r['candidates_attempted']}")
        else:
            print(f"    Applied proof: {r.get('attempted_proof', 'unknown')[:50]}")

    # Save results to JSON for the report
    output_file = Path(__file__).parent / "wright_test_results.json"
    with open(output_file, "w") as f:
        json.dump(all_results, f, indent=2)

    print(f"\n\nDetailed results saved to: {output_file}")
    print("\n" + "=" * 100)

    return 0 if successful > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
