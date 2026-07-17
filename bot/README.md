# Riemann Hardened Bot — Phase 1: SQL Schema & Database Infrastructure

**Status:** Phase 1 increment complete and tested (2026-07-17)

## What's Implemented

Phase 1 of the hardened bot spec includes the **SQL schema and database infrastructure** for the 5-stage automated Lean 4 proof-search pipeline. This is foundational work that enables all downstream components.

### Schema Overview

10 core tables supporting the complete pipeline:

1. **conjectures** — Claude-generated Lean theorem statements, with formalizability scores and status tracking
2. **formalizability_scores** — Lightweight neural predictor output (1-2 sec filter before expensive Lean compilation)
3. **proof_sketches** — Skeleton Lean proofs (with `sorry` placeholders) for fast feedback
4. **critique_sessions** — Optional adversarial self-review stage
5. **formal_attempts** — Lean compilation & verification attempts, with error tracking
6. **proof_strategies** — GPU allocation management (Thompson sampling for adaptive learning)
7. **cross_validations** — Multi-run agreement tracking for high-confidence results
8. **learning_feedback** — Failure reasons fed back to Claude for prompt refinement
9. **lemma_cache** — Reusable proven lemmas for efficiency and convergence
10. **session_log** — High-level metrics and session summaries

### 5-Stage Pipeline Supported

```
Claude (5 candidates)
  ↓ [Stage 2] Formalizability Predictor (1-2 sec, filter P < 0.50)
  ↓ [Stage 3] Proof Sketch Generator (convert strategy → Lean skeleton)
  ↓ [Stage 4] Formal Verification (Lean compilation + error tracking)
  ↓ [Stage 5] Learning & Thompson Sampling (adapt GPU allocation)
  ↓ [Post] Cross-Validation + Lemma Caching
  ↓ Memory & Active Learning Checkpoints
  → Database (all results archived for analysis)
```

### Key Features

- **Deduplication:** statement_hash uniqueness prevents redundant conjectures
- **Error tracking:** Detailed error_type classification for learning loops
- **GPU allocation:** Thompson sampling over 5 proof strategies
- **Lemma caching:** Reusable components reduce Lean compilation time
- **Learning feedback:** Failed conjectures feed back to Claude with categorized failure reasons
- **Cross-validation:** Multiple independent runs for high-stakes conjectures
- **Session metrics:** Pass rates, GPU utilization, strategy convergence

## Files

- **schema.sql** — Full database schema (11 tables, 50+ columns with constraints and indexes)
- **database.py** — Python ORM-like interface with 40+ methods for bot stages
  - `insert_conjecture()` — Stage 1: add conjectures
  - `update_formalizability_score()` — Stage 2: record predictor output
  - `insert_proof_sketch()` — Stage 3: add skeleton proofs
  - `insert_formal_attempt()` / `finalize_formal_attempt()` — Stage 4: attempt tracking
  - `update_strategy_result()` — Thompson sampling state updates
  - `get_session_summary()` — Metrics reporting
  - And many more query helpers
- **test_real_data.py** — Comprehensive end-to-end test with RobinInequality.lean extraction
- **README.md** — This file

## Testing

Run the end-to-end test:

```bash
python3 test_real_data.py
```

Expected output:
- Extracts 2 theorems from RobinInequality.lean (or 2 synthetic examples if file not found)
- Simulates all 5 stages of the pipeline
- Produces session summary with metrics
- Validates query functionality

**Test results (2026-07-17):**
- ✓ Mac: Database initialized, theorems extracted, all stages completed, queries working
- ✓ Olympus: Database on Linux, synthetic theorems, all stages completed, queries working

## Database Location

- **Mac (dev):** `/tmp/riemann_bot.db`
- **Olympus (prod):** `/tmp/riemann_bot.db` (or configurable via `RiemannBotDB(db_path=...)`)

## API Usage

### Quick Start

```python
from database import RiemannBotDB

# Initialize
db = RiemannBotDB(db_path="/tmp/riemann_bot.db")

# Stage 1: Insert conjecture
conj_id = db.insert_conjecture(
    statement="theorem foo : ∀ x, x + 0 = x",
    strategy_type="dfs",
    difficulty_estimate=3.5
)

# Stage 2: Add formalizability score
db.update_formalizability_score(
    conjecture_id=conj_id,
    score=0.78,
    confidence=0.92,
    issues=["needs_lemma_add_zero"]
)

# Stage 3: Add proof sketch
sketch_id = db.insert_proof_sketch(
    conjecture_id=conj_id,
    sketch_code="theorem foo : ∀ x, x + 0 = x := by\n  intro x\n  simp  -- sorry",
    sketch_strategy="direct"
)

# Stage 4: Log formal attempt
attempt_id = db.insert_formal_attempt(
    conjecture_id=conj_id,
    lean_code="theorem foo : ∀ x, x + 0 = x := by intro x; simp",
    strategy_id=1,  # dfs
    sketch_id=sketch_id
)

# Stage 4 (cont): Record result
db.finalize_formal_attempt(
    attempt_id=attempt_id,
    result="success",
    wall_time_ms=2400.0
)

# Query
summary = db.get_session_summary()
print(f"Proved: {summary['proven']}/{summary['total_conjectures']}")

db.close()
```

## Next Phase 1 Increments

Now that the schema is in place, the remaining Phase 1 work (19h total):

1. **Formalizability Predictor** (8h) — Neural network trained on Mathlib corpus
   - Input: Lean theorem statement (text)
   - Output: Score 0-1 (prob. of successful formalization)
   - Uses embedding index to avoid repeated computations
   
2. **Proof Sketch Template Generator** (4h) — Strategy → Lean skeleton
   - 5 templates: direct, induction, contradiction, case_split, generalization
   - Maps Claude strategy to Lean tactic skeleton
   - Compiles in 1-5 sec for fast feedback
   
3. **Embedding Index Pre-computation** (6h) — Mathlib corpus analysis
   - Extract + vectorize all Mathlib theorems
   - Compute similarity matrix
   - Cache for formalizability predictor

4. **Schema Expansion** (1h) — Already mostly done; minor tweaks as phases 2-5 are built

## Indices & Performance

Optimized indices for common queries:
- `idx_conjectures_status` — Stage filtering
- `idx_formalizability_score` — Pass/fail boundary
- `idx_attempts_result` — Aggregate success rate
- `idx_feedback_type` — Learning loop analysis
- `idx_lemma_hits` — Cache effectiveness

## Integration Points

This schema integrates with:
- **Claude SDK** (conjectures → Claude, feedback → Claude)
- **Lean 4 compiler** (formal_attempts table receives stderr)
- **CuPy GPU batch** (gpu_device tracking in attempts)
- **Olympus 2x RTX 3090** (gpu_device = "gpu0" or "gpu1")
- **Thompson sampling allocator** (proof_strategies → GPU allocation)

## Migration & Deployment

To deploy on Olympus:

```bash
scp bot/schema.sql olympus:/tmp/
scp bot/database.py olympus:/tmp/

ssh olympus
python3 << 'EOF'
from database import RiemannBotDB
db = RiemannBotDB(db_path="/tmp/riemann_bot.db")
print("Database ready on Olympus")
EOF
```

## References

- **Spec:** `~/.claude-cobo/projects/-Users-scones/memory/riemann_bot_hardened_spec.md`
- **Lean Files:** `RiemannLean/RiemannLean/RobinInequality.lean` (real theorems to formalize)
- **PROCESS.md:** Full decision history and design rationale
- **STATUS.md:** Numerical results (spectral depth_ratio, etc.)

---

**Author:** Sam (Jacob's MacBook agent)  
**Date:** 2026-07-17  
**Status:** Phase 1 increment (schema) complete; Phase 1 ready for next increments
