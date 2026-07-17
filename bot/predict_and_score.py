#!/usr/bin/env python3
"""
predict_and_score.py — Main driver for formalizability prediction and database scoring.

Usage:
    python3 predict_and_score.py --statement "theorem ..." [--db /path/to/db]
    python3 predict_and_score.py --file file.lean [--db /path/to/db]
"""

import argparse
import sys
import json
from pathlib import Path
from typing import List

from database import RiemannBotDB
from formalizability_predictor import FormalizabilityPredictor


def score_and_insert(
    db: RiemannBotDB,
    statements: List[str],
    strategy_type: str = "dfs",
    source: str = "manual"
) -> None:
    """Score statements and insert into database."""
    predictor = FormalizabilityPredictor()

    for statement in statements:
        print(f"\nProcessing statement ({len(statement)} chars):")
        print(f"  {statement[:80]}..." if len(statement) > 80 else f"  {statement}")

        # Insert conjecture
        conj_id = db.insert_conjecture(
            statement=statement,
            strategy_type=strategy_type,
            difficulty_estimate=5.0,
            source=source
        )
        print(f"  Inserted conjecture ID: {conj_id}")

        # Predict formalizability
        result = predictor.predict(statement)
        print(f"  Score: {result.score:.3f} | Confidence: {result.confidence:.3f}")
        print(f"  Time: {result.processing_time_ms:.1f}ms")

        if result.top_3_issues:
            print(f"  Top issues: {', '.join(result.top_3_issues[:2])}")

        # Update database with score
        db.update_formalizability_score(
            conjecture_id=conj_id,
            score=result.score,
            confidence=result.confidence,
            reason_fail=result.reason_fail,
            issues=result.top_3_issues,
            processing_time_ms=result.processing_time_ms
        )
        print(f"  Status: {'PASS' if result.score >= 0.5 else 'FAIL'} filter")

    # Print session summary
    print("\n" + "=" * 80)
    print("SESSION SUMMARY")
    print("=" * 80)
    summary = db.get_session_summary()
    for key, value in summary.items():
        if isinstance(value, float):
            print(f"  {key:.<40} {value:.2f}")
        else:
            print(f"  {key:.<40} {value}")


def extract_statements_from_lean(filepath: str) -> List[str]:
    """Extract theorem statements from a Lean file."""
    statements = []
    with open(filepath, 'r') as f:
        lines = f.readlines()

    import re

    i = 0
    while i < len(lines):
        line = lines[i]

        # Check if this line starts a definition/theorem/lemma
        if re.match(r"^\s*(theorem|lemma|def|noncomputable)\s+", line):
            # Collect lines until we find := or :=
            stmt_lines = [line]
            i += 1

            while i < len(lines):
                stmt_lines.append(lines[i])
                combined = "".join(stmt_lines)

                # Check if we have a complete statement (found :=)
                if ":=" in combined:
                    # Extract just up to := or the full theorem
                    stmt = combined.strip()
                    if stmt and len(stmt) > 10:  # Filter out trivially small statements
                        statements.append(stmt)
                    break

                i += 1
        else:
            i += 1

    return statements


def main():
    parser = argparse.ArgumentParser(
        description="Score Lean statements for formalizability"
    )
    parser.add_argument(
        "--statement",
        help="Single Lean statement to score"
    )
    parser.add_argument(
        "--file",
        help="Lean file to extract and score statements from"
    )
    parser.add_argument(
        "--db",
        default="/tmp/riemann_bot.db",
        help="Path to SQLite database (default: /tmp/riemann_bot.db)"
    )
    parser.add_argument(
        "--strategy",
        default="dfs",
        help="Proof strategy type (default: dfs)"
    )
    parser.add_argument(
        "--source",
        default="manual",
        help="Source annotation (default: manual)"
    )

    args = parser.parse_args()

    if not args.statement and not args.file:
        parser.print_help()
        sys.exit(1)

    # Open database
    db = RiemannBotDB(db_path=args.db)

    # Collect statements
    statements = []
    if args.statement:
        statements.append(args.statement)
    if args.file:
        extracted = extract_statements_from_lean(args.file)
        statements.extend(extracted)
        print(f"Extracted {len(extracted)} statements from {args.file}")

    # Score and insert
    score_and_insert(
        db=db,
        statements=statements,
        strategy_type=args.strategy,
        source=args.source
    )

    db.close()
    print(f"\nDatabase saved to: {args.db}")


if __name__ == "__main__":
    main()
