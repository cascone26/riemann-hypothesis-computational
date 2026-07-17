#!/usr/bin/env python3
"""
test_formalizability_integration.py — End-to-end test of formalizability predictor.

Tests:
1. Predictor accuracy on calibrated test cases
2. Database integration (insert, score, retrieve)
3. Real-world validation on RobinInequality.lean
4. Performance benchmarks
"""

import time
import sys
from pathlib import Path

from database import RiemannBotDB
from formalizability_predictor import FormalizabilityPredictor


def test_predictor_accuracy():
    """Test predictor on known calibration cases."""
    print("\n" + "=" * 80)
    print("TEST 1: Predictor Accuracy")
    print("=" * 80)

    predictor = FormalizabilityPredictor()

    test_cases = [
        # High-difficulty formalizable (expect 0.65-0.85)
        (
            "theorem robin_inequality_of_RH (hRH : RiemannHypothesis) : "
            "∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n := by sorry",
            0.65,
            0.15
        ),
        # Easy/trivial (expect 0.80+)
        (
            "theorem simple_eq : 1 + 1 = 2 := by norm_num",
            0.80,
            0.10
        ),
        # Questionable (expect 0.35-0.60)
        (
            "theorem undefined_mystery : mystery_func 42 = mystery_result := by sorry",
            0.45,
            0.15
        ),
        # Low/broken (expect <0.40)
        (
            "theorem broken_syntax : unknown_mystery_function ( ( ( := by sorry",
            0.25,
            0.15
        ),
    ]

    passed = 0
    failed = 0

    for statement, target, tolerance in test_cases:
        result = predictor.predict(statement)
        diff = abs(result.score - target)
        status = "PASS" if diff <= tolerance else "FAIL"

        if status == "PASS":
            passed += 1
        else:
            failed += 1

        print(f"\n  {status}: Score {result.score:.3f} vs target {target:.3f} (±{tolerance:.3f})")
        print(f"    {statement[:60]}...")
        print(f"    Confidence: {result.confidence:.3f}, Time: {result.processing_time_ms:.1f}ms")

    print(f"\nAccuracy: {passed}/{len(test_cases)} passed")
    return failed == 0


def test_database_integration():
    """Test database insert, score, and retrieval."""
    print("\n" + "=" * 80)
    print("TEST 2: Database Integration")
    print("=" * 80)

    # Use in-memory database for testing
    db = RiemannBotDB(db_path=":memory:")
    predictor = FormalizabilityPredictor()

    # Test 1: Insert conjecture
    conj_id = db.insert_conjecture(
        statement="theorem test : 2 + 2 = 4 := by norm_num",
        strategy_type="dfs",
        difficulty_estimate=1.0,
        source="test"
    )
    print(f"  ✓ Inserted conjecture ID: {conj_id}")

    # Test 2: Score conjecture
    result = predictor.predict("theorem test : 2 + 2 = 4 := by norm_num")
    db.update_formalizability_score(
        conjecture_id=conj_id,
        score=result.score,
        confidence=result.confidence,
        reason_fail=None,
        issues=result.top_3_issues,
        processing_time_ms=result.processing_time_ms
    )
    print(f"  ✓ Scored conjecture: {result.score:.3f}")

    # Test 3: Query passed filter
    passed = db.get_passed_formalizability(limit=5)
    assert len(passed) == 1, f"Expected 1 passed conjecture, got {len(passed)}"
    print(f"  ✓ Retrieved passed conjecture: {passed[0]['id']}")

    # Test 4: Get session summary
    summary = db.get_session_summary()
    assert summary["total_conjectures"] == 1, "Should have 1 conjecture"
    assert summary["passed_filter"] == 1, "Should have passed filter"
    print(f"  ✓ Session summary: {summary['passed_filter']} passed, {summary['total_conjectures']} total")

    db.close()
    return True


def test_robin_inequality():
    """Test on real RobinInequality.lean file."""
    print("\n" + "=" * 80)
    print("TEST 3: Real-World Validation (RobinInequality.lean)")
    print("=" * 80)

    db = RiemannBotDB(db_path=":memory:")
    predictor = FormalizabilityPredictor()

    # Real theorems from RobinInequality.lean
    robin_theorems = [
        "theorem robin_inequality_of_RH (hRH : RiemannHypothesis) : ∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n := by sorry",
        "theorem RH_of_robin_inequality (h : ∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n) : RiemannHypothesis := by sorry",
        "theorem robin_iff_RH : RiemannHypothesis ↔ ∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n := by sorry",
    ]

    scores = []
    for stmt in robin_theorems:
        result = predictor.predict(stmt)
        scores.append(result.score)
        print(f"  Score: {result.score:.3f} | {stmt[:60]}...")

    avg_score = sum(scores) / len(scores)
    print(f"\nAverage score: {avg_score:.3f}")
    print(f"All passed filter (>0.5): {all(s >= 0.5 for s in scores)}")

    db.close()
    return all(s >= 0.5 for s in scores)


def test_performance():
    """Benchmark predictor performance."""
    print("\n" + "=" * 80)
    print("TEST 4: Performance Benchmark")
    print("=" * 80)

    predictor = FormalizabilityPredictor()

    statements = [
        "theorem simple : 1 = 1 := rfl",
        "theorem medium : ∀ x y : ℕ, x + y = y + x := by ring",
        "theorem complex : ∀ a b c d e : ℕ, a + b = c + d → e + a + b = e + c + d := by simp",
    ] * 10  # Run 30 total predictions

    times = []
    start = time.time()

    for stmt in statements:
        result = predictor.predict(stmt)
        times.append(result.processing_time_ms)

    elapsed = time.time() - start
    avg_time = sum(times) / len(times)
    max_time = max(times)
    min_time = min(times)

    print(f"\n  Total time for {len(statements)} predictions: {elapsed:.2f}s")
    print(f"  Average: {avg_time:.2f}ms")
    print(f"  Min: {min_time:.2f}ms")
    print(f"  Max: {max_time:.2f}ms")
    print(f"  Target: <2000ms per statement, <1-2s total")

    return avg_time < 2.0 and max_time < 5.0


def main():
    """Run all tests."""
    print("\n")
    print("╔" + "=" * 78 + "╗")
    print("║" + "FORMALIZABILITY PREDICTOR — INTEGRATION TEST SUITE".center(78) + "║")
    print("╚" + "=" * 78 + "╝")

    tests = [
        ("Predictor Accuracy", test_predictor_accuracy),
        ("Database Integration", test_database_integration),
        ("Real-World Validation", test_robin_inequality),
        ("Performance", test_performance),
    ]

    results = {}
    for name, test_func in tests:
        try:
            result = test_func()
            results[name] = "PASS" if result else "FAIL"
        except Exception as e:
            print(f"\n  ERROR: {e}")
            results[name] = "ERROR"

    # Summary
    print("\n" + "=" * 80)
    print("SUMMARY")
    print("=" * 80)
    for name, result in results.items():
        status = "✓" if result == "PASS" else "✗"
        print(f"  {status} {name:.<40} {result}")

    print("\n" + "=" * 80)

    all_passed = all(r == "PASS" for r in results.values())
    if all_passed:
        print("ALL TESTS PASSED")
        return 0
    else:
        print("SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    sys.exit(main())
