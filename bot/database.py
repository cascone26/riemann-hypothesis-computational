#!/usr/bin/env python3
"""
database.py — SQLite database setup and management for Riemann hardened bot.

Provides:
- Database initialization and schema creation
- Connection pooling
- Utility functions for bot stages (conjecture insertion, filtering, attempt logging)
- Query helpers for learning feedback loops
"""

import sqlite3
import os
import json
import hashlib
from datetime import datetime, timedelta
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple


class RiemannBotDB:
    """SQLite database manager for hardened bot pipeline."""

    def __init__(self, db_path: str = "/tmp/riemann_bot.db", schema_path: Optional[str] = None):
        """
        Initialize database connection.

        Args:
            db_path: Path to SQLite database file
            schema_path: Path to schema.sql (if None, use default relative path)
        """
        self.db_path = db_path
        self.schema_path = schema_path or Path(__file__).parent / "schema.sql"
        self.conn = None
        self._init_db()

    def _init_db(self):
        """Create database and apply schema if needed."""
        self.conn = sqlite3.connect(self.db_path)
        self.conn.row_factory = sqlite3.Row

        # Check if tables already exist
        cursor = self.conn.cursor()
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='conjectures'")
        exists = cursor.fetchone() is not None

        if not exists:
            self._apply_schema()
            self._init_default_strategies()
            print(f"Database initialized: {self.db_path}")
        else:
            print(f"Database exists: {self.db_path}")

    def _apply_schema(self):
        """Apply SQL schema from file."""
        with open(self.schema_path, 'r') as f:
            schema_sql = f.read()

        cursor = self.conn.cursor()
        cursor.executescript(schema_sql)
        self.conn.commit()

    def _init_default_strategies(self):
        """Initialize default proof strategies."""
        strategies = [
            ("dfs", "Depth-First Search: exhaustive backward chaining"),
            ("bfs", "Breadth-First Search: systematic forward exploration"),
            ("best_first", "Best-First: greedy heuristic using tactic success rate"),
            ("mcts", "Monte Carlo Tree Search: random playouts with UCB1"),
            ("iterative_deepening", "Iterative Deepening: bounded depth increase"),
        ]

        cursor = self.conn.cursor()
        for name, desc in strategies:
            cursor.execute(
                "INSERT OR IGNORE INTO proof_strategies (name, description, gpu_allocation_fraction) "
                "VALUES (?, ?, 0.2)",
                (name, desc)
            )
        self.conn.commit()

    def close(self):
        """Close database connection."""
        if self.conn:
            self.conn.close()

    # ========================
    # Stage 1: Conjecture Input
    # ========================

    def insert_conjecture(
        self,
        statement: str,
        strategy_type: str,
        difficulty_estimate: float = 5.0,
        source: str = "claude_generated"
    ) -> int:
        """
        Insert a new conjecture from Claude.

        Args:
            statement: Lean theorem statement
            strategy_type: Proof strategy type ('dfs', 'bfs', etc.)
            difficulty_estimate: 1-10 difficulty score
            source: Where this came from

        Returns:
            conjecture_id
        """
        stmt_hash = hashlib.sha256(statement.encode()).hexdigest()

        cursor = self.conn.cursor()
        cursor.execute(
            """
            INSERT INTO conjectures
            (statement, statement_hash, strategy_type, difficulty_estimate, source, status)
            VALUES (?, ?, ?, ?, ?, 'pending')
            """,
            (statement, stmt_hash, strategy_type, difficulty_estimate, source)
        )
        self.conn.commit()
        return cursor.lastrowid

    def get_pending_conjectures(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get pending conjectures ready for filtering."""
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT * FROM conjectures
            WHERE status = 'pending' AND formalizability_score IS NULL
            LIMIT ?
            """,
            (limit,)
        )
        return [dict(row) for row in cursor.fetchall()]

    # ================================
    # Stage 2: Formalizability Scoring
    # ================================

    def update_formalizability_score(
        self,
        conjecture_id: int,
        score: float,
        confidence: float,
        reason_fail: Optional[str] = None,
        issues: Optional[List[str]] = None,
        processing_time_ms: float = 0.0
    ) -> None:
        """Record formalizability predictor output."""
        top_3 = json.dumps(issues[:3]) if issues else "[]"

        cursor = self.conn.cursor()
        # Update conjecture status
        new_status = "filtered" if score < 0.5 else "sketched"
        cursor.execute(
            "UPDATE conjectures SET status = ?, formalizability_score = ?, formalizability_confidence = ? "
            "WHERE id = ?",
            (new_status, score, confidence, conjecture_id)
        )

        # Insert score record
        cursor.execute(
            """
            INSERT INTO formalizability_scores
            (conjecture_id, score, confidence, reason_fail, top_3_issues, processing_time_ms)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (conjecture_id, score, confidence, reason_fail, top_3, processing_time_ms)
        )
        self.conn.commit()

    def get_passed_formalizability(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get conjectures that passed formalizability filter (score >= 0.5)."""
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT c.*, f.score as formalizability_score
            FROM conjectures c
            LEFT JOIN formalizability_scores f ON c.id = f.conjecture_id
            WHERE f.score >= 0.5 AND c.status IN ('sketched', 'attempted')
            ORDER BY f.score DESC
            LIMIT ?
            """,
            (limit,)
        )
        return [dict(row) for row in cursor.fetchall()]

    # =============================
    # Stage 3: Proof Sketch Generator
    # =============================

    def insert_proof_sketch(
        self,
        conjecture_id: int,
        sketch_code: str,
        sketch_strategy: str,
        generated_from_model: str = "claude_3_5_sonnet"
    ) -> int:
        """Insert a proof sketch (Lean skeleton with `sorry` placeholders)."""
        sketch_hash = hashlib.sha256(sketch_code.encode()).hexdigest()

        cursor = self.conn.cursor()
        cursor.execute(
            """
            INSERT INTO proof_sketches
            (conjecture_id, sketch_code, sketch_hash, sketch_strategy, generated_from_model)
            VALUES (?, ?, ?, ?, ?)
            """,
            (conjecture_id, sketch_code, sketch_hash, sketch_strategy, generated_from_model)
        )
        self.conn.commit()

        # Update conjecture to point to this sketch
        cursor.execute("UPDATE conjectures SET proof_sketch_id = ? WHERE id = ?",
                      (cursor.lastrowid, conjecture_id))
        self.conn.commit()

        return cursor.lastrowid

    def update_sketch_compilation_status(
        self,
        sketch_id: int,
        status: str,
        compilation_time_ms: float
    ) -> None:
        """Update sketch compilation result."""
        cursor = self.conn.cursor()
        cursor.execute(
            """
            UPDATE proof_sketches
            SET compilation_status = ?, compilation_time_ms = ?
            WHERE id = ?
            """,
            (status, compilation_time_ms, sketch_id)
        )
        self.conn.commit()

    # =============================
    # Stage 4: Formal Verification
    # =============================

    def insert_formal_attempt(
        self,
        conjecture_id: int,
        lean_code: str,
        strategy_id: int,
        gpu_device: str = "cpu",
        attempt_type: str = "full_proof",
        sketch_id: Optional[int] = None
    ) -> int:
        """Record the start of a formal verification attempt."""
        cursor = self.conn.cursor()
        cursor.execute(
            """
            INSERT INTO formal_attempts
            (conjecture_id, sketch_id, attempt_type, lean_code, strategy_id, gpu_device, started_at)
            VALUES (?, ?, ?, ?, ?, ?, CURRENT_TIMESTAMP)
            """,
            (conjecture_id, sketch_id, attempt_type, lean_code, strategy_id, gpu_device)
        )
        self.conn.commit()

        # Update conjecture status
        cursor.execute("UPDATE conjectures SET status = 'attempted' WHERE id = ?", (conjecture_id,))
        self.conn.commit()

        return cursor.lastrowid

    def finalize_formal_attempt(
        self,
        attempt_id: int,
        result: str,
        error_message: Optional[str] = None,
        error_type: Optional[str] = None,
        error_line: Optional[int] = None,
        wall_time_ms: float = 0.0,
        lemma_hits: int = 0,
        lemma_misses: int = 0
    ) -> None:
        """Record the result of a formal verification attempt."""
        cursor = self.conn.cursor()
        cursor.execute(
            """
            UPDATE formal_attempts
            SET ended_at = CURRENT_TIMESTAMP, result = ?, error_message = ?,
                error_type = ?, error_line = ?, wall_time_ms = ?,
                lemma_cache_hits = ?, lemma_cache_misses = ?
            WHERE id = ?
            """,
            (result, error_message, error_type, error_line, wall_time_ms, lemma_hits, lemma_misses, attempt_id)
        )
        self.conn.commit()

        # Update conjecture final result
        cursor.execute(
            "SELECT conjecture_id FROM formal_attempts WHERE id = ?",
            (attempt_id,)
        )
        conj_id = cursor.fetchone()[0]

        if result == "success":
            cursor.execute("UPDATE conjectures SET status = 'proven', final_result = 'success' WHERE id = ?",
                          (conj_id,))
        else:
            cursor.execute("UPDATE conjectures SET final_result = ? WHERE id = ?", (result, conj_id))

        self.conn.commit()

    # =======================
    # Strategy Management
    # =======================

    def update_strategy_result(self, strategy_id: int, success: bool) -> None:
        """Record attempt outcome for Thompson sampling."""
        cursor = self.conn.cursor()

        if success:
            cursor.execute(
                """
                UPDATE proof_strategies
                SET success_count = success_count + 1, total_attempts = total_attempts + 1
                WHERE id = ?
                """,
                (strategy_id,)
            )
        else:
            cursor.execute(
                """
                UPDATE proof_strategies
                SET failure_count = failure_count + 1, total_attempts = total_attempts + 1
                WHERE id = ?
                """,
                (strategy_id,)
            )

        self.conn.commit()
        self._update_strategy_win_rate(strategy_id)

    def _update_strategy_win_rate(self, strategy_id: int) -> None:
        """Recompute win rate for a strategy."""
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT success_count, total_attempts FROM proof_strategies WHERE id = ?",
            (strategy_id,)
        )
        row = cursor.fetchone()

        if row and row[1] > 0:
            win_rate = row[0] / row[1]
            cursor.execute("UPDATE proof_strategies SET win_rate = ? WHERE id = ?", (win_rate, strategy_id))
            self.conn.commit()

    def get_strategy_allocations(self) -> Dict[str, float]:
        """Get current GPU allocation fractions for all strategies."""
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT name, gpu_allocation_fraction FROM proof_strategies ORDER BY gpu_allocation_fraction DESC"
        )
        return {row[0]: row[1] for row in cursor.fetchall()}

    # =======================
    # Learning Feedback
    # =======================

    def insert_feedback(
        self,
        conjecture_id: int,
        feedback_type: str,
        message: str,
        severity: str = "medium",
        suggested_fix: Optional[str] = None,
        attempt_id: Optional[int] = None
    ) -> None:
        """Record learning feedback for Claude to process."""
        cursor = self.conn.cursor()
        cursor.execute(
            """
            INSERT INTO learning_feedback
            (conjecture_id, attempt_id, feedback_type, message, suggested_fix, severity)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (conjecture_id, attempt_id, feedback_type, message, suggested_fix, severity)
        )
        self.conn.commit()

    def get_unused_feedback(self, limit: int = 10) -> List[Dict[str, Any]]:
        """Get recent feedback not yet used in Claude prompts."""
        cursor = self.conn.cursor()
        cursor.execute(
            """
            SELECT lf.*, c.statement, fa.lean_code, fa.error_message
            FROM learning_feedback lf
            JOIN conjectures c ON lf.conjecture_id = c.id
            LEFT JOIN formal_attempts fa ON lf.attempt_id = fa.id
            WHERE lf.used_in_prompt = FALSE
            ORDER BY lf.severity DESC, lf.created_at DESC
            LIMIT ?
            """,
            (limit,)
        )
        return [dict(row) for row in cursor.fetchall()]

    def mark_feedback_used(self, feedback_id: int) -> None:
        """Mark feedback as incorporated into a prompt."""
        cursor = self.conn.cursor()
        cursor.execute("UPDATE learning_feedback SET used_in_prompt = TRUE WHERE id = ?", (feedback_id,))
        self.conn.commit()

    # =======================
    # Metrics & Reporting
    # =======================

    def get_session_summary(self) -> Dict[str, Any]:
        """Get summary metrics for current session."""
        cursor = self.conn.cursor()

        cursor.execute("SELECT COUNT(*) FROM conjectures")
        total_conj = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM conjectures WHERE final_result = 'success'")
        proven = cursor.fetchone()[0]

        cursor.execute("SELECT COUNT(*) FROM conjectures WHERE formalizability_score >= 0.5")
        passed_filter = cursor.fetchone()[0]

        cursor.execute("SELECT AVG(score) FROM formalizability_scores")
        avg_score = cursor.fetchone()[0]

        cursor.execute("SELECT SUM(wall_time_ms) FROM formal_attempts WHERE result = 'success'")
        total_gpu_ms = cursor.fetchone()[0] or 0

        cursor.execute(
            "SELECT name FROM proof_strategies ORDER BY success_count DESC LIMIT 1"
        )
        best_strategy = cursor.fetchone()

        return {
            "total_conjectures": total_conj,
            "proven": proven,
            "passed_filter": passed_filter,
            "formalizability_pass_rate": (passed_filter / total_conj * 100) if total_conj > 0 else 0,
            "avg_formalizability_score": avg_score,
            "total_gpu_time_ms": total_gpu_ms,
            "best_strategy": best_strategy[0] if best_strategy else None
        }

    def get_error_summary(self) -> Dict[str, int]:
        """Count error types for feedback loops."""
        cursor = self.conn.cursor()
        cursor.execute(
            "SELECT error_type, COUNT(*) FROM formal_attempts WHERE error_type IS NOT NULL GROUP BY error_type"
        )
        return {row[0]: row[1] for row in cursor.fetchall()}


def init_db_for_testing() -> RiemannBotDB:
    """Create in-memory test database."""
    db = RiemannBotDB(db_path=":memory:")
    return db


if __name__ == "__main__":
    # Test database setup
    db = RiemannBotDB()

    # Insert a test conjecture (from RobinInequality.lean)
    conj_id = db.insert_conjecture(
        statement="theorem robin_inequality_of_RH (hRH : RiemannHypothesis) : ∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n",
        strategy_type="dfs",
        difficulty_estimate=7.5,
        source="riemann_mathlib"
    )
    print(f"Inserted conjecture: {conj_id}")

    # Simulate formalizability scoring
    db.update_formalizability_score(
        conj_id,
        score=0.72,
        confidence=0.85,
        reason_fail=None,
        issues=["needs_lemma_sigma_properties", "euler_mascheroni_import"],
        processing_time_ms=1350
    )
    print(f"Scored conjecture {conj_id}: 0.72")

    # Insert a proof sketch
    sketch_id = db.insert_proof_sketch(
        conj_id,
        sketch_code="""theorem robin_inequality_of_RH (hRH : RiemannHypothesis) :
    ∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n := by
  intro n hn
  sorry  -- Use positivity of log(log(n)) for n > 5040
  """,
        sketch_strategy="direct"
    )
    print(f"Inserted sketch: {sketch_id}")

    # Check retrieval
    pending = db.get_pending_conjectures()
    print(f"Pending conjectures: {len(pending)}")

    passed = db.get_passed_formalizability()
    print(f"Passed filter: {len(passed)}")

    # Get session summary
    summary = db.get_session_summary()
    print(f"Session summary: {summary}")

    db.close()
    print(f"\nTest complete. Database: {db.db_path}")
