-- Riemann Hardened Bot — Phase 1 SQL Schema
-- Database for automated Lean 4 proof-search pipeline
-- Created: 2026-07-17

-- Core conjecture data: Claude-generated Lean statements
CREATE TABLE IF NOT EXISTS conjectures (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    statement TEXT NOT NULL,              -- The Lean conjecture statement
    statement_hash TEXT UNIQUE,           -- SHA256 of statement (dedup)
    strategy_type TEXT NOT NULL,          -- 'dfs', 'bfs', 'best_first', 'mcts', 'iterative_deepening'
    difficulty_estimate REAL,             -- 1-10 score from Claude
    source TEXT,                          -- 'claude_dfs', 'claude_random', etc.
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    status TEXT DEFAULT 'pending',        -- 'pending', 'filtered', 'sketched', 'attempted', 'proven', 'failed'
    formalizability_score REAL,           -- [0, 1] from predictor (filled by stage 2)
    formalizability_confidence REAL,      -- predictor confidence
    proof_sketch_id INTEGER,              -- FK to proof_sketches
    final_result TEXT,                    -- 'success', 'timeout', 'type_error', 'elaboration_error', 'tactic_failed'
    lean_error_msg TEXT,                  -- Lean compiler stderr (if failed)
    lemma_names TEXT,                     -- Comma-separated lemmas used (for caching)
    iteration_count INTEGER DEFAULT 0     -- How many times Claude revisited this
);
CREATE INDEX IF NOT EXISTS idx_conjectures_status ON conjectures(status);
CREATE INDEX IF NOT EXISTS idx_conjectures_strategy ON conjectures(strategy_type);
CREATE INDEX IF NOT EXISTS idx_conjectures_created ON conjectures(created_at);


-- Formalizability predictor outputs: stage 2 filter
CREATE TABLE IF NOT EXISTS formalizability_scores (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    conjecture_id INTEGER NOT NULL UNIQUE,
    score REAL NOT NULL,                 -- [0, 1], >0.5 passes filter
    confidence REAL NOT NULL,            -- predictor confidence in score
    reason_fail TEXT,                    -- Why it might fail: 'syntax_invalid', 'type_mismatch', 'tactic_inapplicable', etc.
    top_3_issues TEXT,                   -- JSON array of top failure reasons
    processing_time_ms REAL,             -- Wall time for prediction
    model_version TEXT,                  -- predictor version that generated score
    computed_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (conjecture_id) REFERENCES conjectures(id)
);
CREATE INDEX IF NOT EXISTS idx_formalizability_score ON formalizability_scores(score DESC);


-- Proof sketches: lightweight Lean skeletons (compiles in 1-5s)
CREATE TABLE IF NOT EXISTS proof_sketches (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    conjecture_id INTEGER NOT NULL,
    sketch_code TEXT NOT NULL,           -- Lean skeleton with `sorry` placeholders
    sketch_hash TEXT UNIQUE,             -- SHA256 dedup
    sketch_strategy TEXT,                -- Strategy used to generate: 'direct', 'induction', 'contradiction', 'case_split', 'generalization'
    generated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    compilation_status TEXT,             -- 'unknown', 'compiles', 'syntax_error', 'type_error'
    compilation_time_ms REAL,            -- Wall time for compilation
    generated_from_model TEXT,           -- 'claude_3_5_sonnet', etc.
    iteration_count INTEGER DEFAULT 0,   -- Refinement count
    FOREIGN KEY (conjecture_id) REFERENCES conjectures(id)
);
CREATE INDEX IF NOT EXISTS idx_sketches_conjecture ON proof_sketches(conjecture_id);
CREATE INDEX IF NOT EXISTS idx_sketches_status ON proof_sketches(compilation_status);


-- Adversarial critique: optional self-review stage
CREATE TABLE IF NOT EXISTS critique_sessions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    conjecture_id INTEGER NOT NULL,
    sketch_id INTEGER NOT NULL,
    critique_text TEXT,                 -- Claude's self-review of the proof
    issues_found INTEGER DEFAULT 0,     -- Count of logical/structural issues
    confidence_delta REAL,               -- How much critique changed confidence
    actions_suggested TEXT,              -- JSON: ["strengthen_invariant", "add_case", "use_different_tactic", ...]
    resolved BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (conjecture_id) REFERENCES conjectures(id),
    FOREIGN KEY (sketch_id) REFERENCES proof_sketches(id)
);


-- Formal verification attempts: Lean compilation & validation
CREATE TABLE IF NOT EXISTS formal_attempts (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    conjecture_id INTEGER NOT NULL,
    sketch_id INTEGER,                  -- FK to proof_sketches (can be NULL for direct attempts)
    attempt_type TEXT NOT NULL,         -- 'full_proof', 'sketch_only', 'tactic_debug'
    lean_code TEXT,                     -- Full Lean theorem + proof
    strategy_id INTEGER,                -- Which GPU strategy was used
    gpu_device TEXT,                    -- 'gpu0', 'gpu1', or 'cpu'
    started_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    ended_at TIMESTAMP,
    wall_time_ms REAL,                  -- Total compute time
    lean_version TEXT,                  -- e.g., "4.29.0"
    result TEXT DEFAULT 'in_progress',  -- 'success', 'timeout', 'type_error', 'elaboration_error', 'tactic_failed', 'other', 'in_progress'
    error_message TEXT,                 -- Lean stderr
    error_type TEXT,                    -- 'syntax', 'type_mismatch', 'tactic', 'elaboration', 'other'
    error_line INTEGER,                 -- Line number in Lean code
    lemma_cache_hits INTEGER DEFAULT 0,
    lemma_cache_misses INTEGER DEFAULT 0,
    rollback_count INTEGER DEFAULT 0,   -- How many times we had to retry
    FOREIGN KEY (conjecture_id) REFERENCES conjectures(id),
    FOREIGN KEY (sketch_id) REFERENCES proof_sketches(id),
    FOREIGN KEY (strategy_id) REFERENCES proof_strategies(id)
);
CREATE INDEX IF NOT EXISTS idx_attempts_conjecture ON formal_attempts(conjecture_id);
CREATE INDEX IF NOT EXISTS idx_attempts_result ON formal_attempts(result);
CREATE INDEX IF NOT EXISTS idx_attempts_started ON formal_attempts(started_at DESC);


-- Proof strategies: GPU allocation & learning (Thompson sampling)
CREATE TABLE IF NOT EXISTS proof_strategies (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT UNIQUE NOT NULL,          -- 'dfs', 'bfs', 'best_first', 'mcts', 'iterative_deepening'
    description TEXT,
    gpu_allocation_fraction REAL DEFAULT 0.2,  -- [0,1]: how much GPU to allocate
    success_count INTEGER DEFAULT 0,    -- How many successful attempts
    failure_count INTEGER DEFAULT 0,    -- How many failed attempts
    total_attempts INTEGER DEFAULT 0,
    win_rate REAL,                      -- success_count / total_attempts
    thompson_alpha INTEGER DEFAULT 1,   -- Thompson sampling: successes + 1
    thompson_beta INTEGER DEFAULT 1,    -- Thompson sampling: failures + 1
    last_updated TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    estimated_win_prob REAL,            -- Beta-binomial sampled win probability
    next_gpu_reallocation TIMESTAMP     -- When to recompute Thompson estimate
);


-- Cross-validation results: multiple runs on same conjecture
CREATE TABLE IF NOT EXISTS cross_validations (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    conjecture_id INTEGER NOT NULL,
    attempt1_id INTEGER NOT NULL,
    attempt2_id INTEGER,               -- Second independent run (NULL if not yet done)
    attempt3_id INTEGER,               -- Third run (for high-value conjectures)
    agreement BOOLEAN,                 -- Do all attempts agree on result?
    consensus_result TEXT,             -- 'success', 'failure', 'uncertain'
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (conjecture_id) REFERENCES conjectures(id),
    FOREIGN KEY (attempt1_id) REFERENCES formal_attempts(id),
    FOREIGN KEY (attempt2_id) REFERENCES formal_attempts(id),
    FOREIGN KEY (attempt3_id) REFERENCES formal_attempts(id)
);


-- Learning feedback loops: feed failures back to Claude
CREATE TABLE IF NOT EXISTS learning_feedback (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    conjecture_id INTEGER NOT NULL,
    attempt_id INTEGER,                -- FK to formal_attempts (NULL if predictor-stage)
    feedback_type TEXT NOT NULL,       -- 'formalizability_failure', 'type_error', 'tactic_failure', 'sketch_suggests_change'
    message TEXT,                      -- Description of failure
    suggested_fix TEXT,                -- How to fix it (for Claude feedback)
    severity TEXT,                     -- 'low', 'medium', 'high'
    used_in_prompt BOOLEAN DEFAULT FALSE,
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (conjecture_id) REFERENCES conjectures(id),
    FOREIGN KEY (attempt_id) REFERENCES formal_attempts(id)
);
CREATE INDEX IF NOT EXISTS idx_feedback_type ON learning_feedback(feedback_type);
CREATE INDEX IF NOT EXISTS idx_feedback_used ON learning_feedback(used_in_prompt);


-- Lemma cache: reusable lemmas discovered during proofs
CREATE TABLE IF NOT EXISTS lemma_cache (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    lemma_name TEXT UNIQUE NOT NULL,   -- e.g., "robin_bound_positive"
    lemma_statement TEXT NOT NULL,     -- The Lean statement
    lemma_proof TEXT,                  -- The full Lean proof (if proven)
    proof_attempt_id INTEGER,          -- Where it was discovered
    hit_count INTEGER DEFAULT 0,       -- How many times reused
    last_used TIMESTAMP,
    mathlib_reference TEXT,            -- If it's a Mathlib lemma
    domain TEXT,                       -- 'number_theory', 'algebra', 'analysis', etc.
    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (proof_attempt_id) REFERENCES formal_attempts(id)
);
CREATE INDEX IF NOT EXISTS idx_lemma_hits ON lemma_cache(hit_count DESC);


-- Session log: high-level metrics for monitoring
CREATE TABLE IF NOT EXISTS session_log (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    session_start TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
    session_end TIMESTAMP,
    conjectures_generated INTEGER,
    conjectures_passed_filter INTEGER,
    conjectures_proven INTEGER,
    total_gpu_time_ms REAL,
    avg_formalizability_score REAL,
    successful_strategy_id INTEGER,    -- Which strategy won this session?
    notes TEXT,
    FOREIGN KEY (successful_strategy_id) REFERENCES proof_strategies(id)
);
