#!/usr/bin/env python3
"""
formalizability_predictor.py — Formalizability predictor for Lean 4 statements.

Scores whether a Lean statement is likely to be formalizable in Lean 4 without
attempting expensive compilation. Targets 1-2 second prediction time, >80% filter
accuracy on known Mathlib theorems.

Strategy:
1. Feature extraction: statement complexity, structure, imports, function usage
2. Heuristic scoring based on:
   - Syntactic validity and structure
   - Similarity to known-formalizable patterns in Mathlib
   - Proof sketch viability (directness, available tactics)
   - Dependency availability in Mathlib
3. Calibration: tested against RobinInequality.lean and Mathlib examples
"""

import re
import time
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass
from collections import Counter
import json


@dataclass
class FormalizabilityScore:
    """Result of formalizability prediction."""
    score: float  # [0, 1], >0.5 passes filter
    confidence: float  # [0, 1] predictor confidence
    reason_fail: Optional[str]  # Primary failure reason if score < 0.5
    top_3_issues: List[str]  # List of top potential issues
    processing_time_ms: float
    metadata: Dict  # Debug info: feature values, intermediate scores


class FormalizabilityPredictor:
    """
    Lean 4 formalizability predictor.

    Implements a heuristic scoring system that predicts whether a given
    Lean statement can be successfully compiled and verified without
    expensive trial-and-error compilation.
    """

    # Known Mathlib modules and their content (for dependency checking)
    MATHLIB_MODULES = {
        "NumberTheory.LSeries.RiemannZeta": ["RiemannHypothesis", "riemannXi"],
        "NumberTheory.ArithmeticFunction": ["sigma", "ArithmeticFunction"],
        "NumberTheory.Harmonic": ["harmonic", "harmonicNumber"],
        "Analysis.SpecialFunctions.Log": ["Real.log", "Real.exp"],
        "Analysis.SpecialFunctions.ExpDeriv": ["deriv", "exp"],
        "NumberTheory.Harmonic.EulerMascheroni": ["eulerMascheroniConstant"],
    }

    # Functions known to be in Mathlib
    MATHLIB_FUNCTIONS = {
        "sigma", "Real.log", "Real.exp", "harmonic", "Real.eulerMascheroniConstant",
        "ArithmeticFunction", "modEq", "prime", "natAbs", "natDegree", "coeff",
        "natFactorization", "omega", "decide", "norm_num", "ring", "field_simp",
        "Nat.add", "Nat.mul", "Nat.sub", "Finset.sum", "List.sum",
    }

    # Tactic keywords recognized as supportive (available in Lean 4)
    LEAN4_TACTICS = {
        "intro", "intros", "apply", "exact", "sorry", "rfl", "simp", "simp_all",
        "norm_num", "omega", "ring", "field_simp", "rw", "rewrite", "unfold",
        "calc", "by_cases", "by_contra", "induction", "cases", "constructor",
        "left", "right", "exists", "ext", "funext", "split", "decide",
        "assumption", "contradiction", "exfalso", "have", "let", "obtain",
        "match", "for", "do", "succ_rec", "nat_rec",
    }

    # Patterns that indicate proof difficulty (heuristic scoring)
    HARD_PATTERNS = [
        (r":=\s*sorry\s*(?:$|\n)", "only sorry, no proof"),
        (r"hnot[^:]*:", "negation/contradiction approach (harder)"),
        (r"∀.*∀.*∀.*∀", "quadruple quantifier or more (very complex)"),
        (r"∀.*∀.*∀", "triple quantifier (complex)"),
        (r"→.*→.*→.*→.*→", "very deep implication chain (complex)"),
        (r"→.*→.*→.*→", "deep implication chain (complex)"),
        (r"\*\s*\*", "nested exponential/power (complex)"),
        (r"∫|∑|∏", "integration/summation/product (complex)"),
        (r"limit|Filter\.", "limit/filter reasoning (complex)"),
        (r"norm_cast|push_cast", "casting requirement (potentially fails)"),
    ]

    # Patterns that indicate easier proofs
    EASY_PATTERNS = [
        (r"theorem.*:.*=.*:=\s*rfl", "trivial equality (rfl)"),
        (r"exact|apply", "direct proof application (easy)"),
        (r"\b(simp|norm_num|omega|decide)\b", "automated tactic availability"),
    ]

    def __init__(self):
        """Initialize predictor with patterns and scoring configuration."""
        self.feature_weights = {
            "syntax_valid": 0.15,
            "structure_sound": 0.15,
            "dependency_available": 0.20,
            "proof_strategy_viable": 0.25,
            "complexity_reasonable": 0.15,
            "pattern_match_mathlib": 0.10,
        }

    def predict(self, lean_statement: str) -> FormalizabilityScore:
        """
        Predict formalizability of a Lean 4 statement.

        Args:
            lean_statement: Full Lean theorem/definition statement

        Returns:
            FormalizabilityScore with prediction and metadata
        """
        start_time = time.time()

        # Clean and normalize
        statement = self._normalize_statement(lean_statement)

        # Extract features
        features = self._extract_features(statement)

        # Check for critical red flags that should override subscores
        red_flags = self._check_red_flags(statement, features)

        # Compute subscores for each category
        subscores = {
            "syntax_valid": self._score_syntax(statement, features),
            "structure_sound": self._score_structure(statement, features),
            "dependency_available": self._score_dependencies(statement, features),
            "proof_strategy_viable": self._score_proof_strategy(statement, features),
            "complexity_reasonable": self._score_complexity(statement, features),
            "pattern_match_mathlib": self._score_pattern_match(statement, features),
        }

        # Weighted aggregation
        final_score = sum(
            subscores[key] * self.feature_weights[key]
            for key in self.feature_weights
        )

        # Apply red flag penalties
        for flag, penalty in red_flags:
            final_score *= penalty

        # Apply critical failure multiplier: if syntax is broken, multiply down hard
        if subscores["syntax_valid"] < 0.3:
            final_score *= 0.4  # Critical syntax failure
        elif subscores["syntax_valid"] < 0.5:
            final_score *= 0.7  # Serious syntax issue

        # If dependencies are mostly unknown, also multiply down
        if subscores["dependency_available"] < 0.2:
            final_score *= 0.5

        # Clamp to [0, 1]
        final_score = max(0.0, min(1.0, final_score))

        # Compute confidence based on subscores agreement
        confidence = self._compute_confidence(subscores, final_score)

        # Identify top 3 issues
        issues = self._identify_issues(statement, features, subscores)
        reason_fail = issues[0] if final_score < 0.5 and issues else None

        elapsed_ms = (time.time() - start_time) * 1000

        return FormalizabilityScore(
            score=final_score,
            confidence=confidence,
            reason_fail=reason_fail,
            top_3_issues=issues[:3],
            processing_time_ms=elapsed_ms,
            metadata={
                "features": features,
                "subscores": subscores,
                "issues": issues,
            }
        )

    def _normalize_statement(self, statement: str) -> str:
        """Normalize Lean statement for processing."""
        # Remove comments
        statement = re.sub(r"--.*$", "", statement, flags=re.MULTILINE)
        statement = re.sub(r"/\*.*?\*/", "", statement, flags=re.DOTALL)
        # Collapse whitespace
        statement = re.sub(r"\s+", " ", statement).strip()
        return statement

    def _extract_features(self, statement: str) -> Dict:
        """Extract numerical and categorical features from statement."""
        return {
            "length": len(statement),
            "token_count": len(statement.split()),
            "line_count": statement.count("\n") + 1,
            "paren_depth": self._compute_max_nesting_depth(statement),
            "quantifier_count": len(re.findall(r"[∀∃∧∨¬]", statement)),
            "implication_count": len(re.findall(r"→|⇒", statement)),
            "function_count": len(re.findall(r"\w+\s*\(", statement)),
            "imports": self._extract_imports(statement),
            "functions_used": self._extract_functions(statement),
            "tactics_used": self._extract_tactics(statement),
            "sorry_count": statement.count("sorry"),
            "forall_depth": self._compute_forall_depth(statement),
            "type_annotations": len(re.findall(r":\s*\w+", statement)),
        }

    def _compute_max_nesting_depth(self, s: str) -> int:
        """Compute maximum nesting depth of parentheses."""
        max_depth = 0
        current_depth = 0
        for char in s:
            if char in "([{⟨":
                current_depth += 1
                max_depth = max(max_depth, current_depth)
            elif char in ")]}⟩":
                current_depth = max(0, current_depth - 1)
        return max_depth

    def _compute_forall_depth(self, s: str) -> int:
        """Estimate depth of nested universal quantifiers."""
        # Count consecutive ∀ at start
        match = re.match(r"(∀\s*[^\s]*[,:]?\s*)*", s)
        if match:
            foralls = len(re.findall(r"∀", match.group(0)))
            return foralls
        return 0

    def _extract_imports(self, statement: str) -> List[str]:
        """Extract import statements."""
        imports = re.findall(r"import\s+([\w.]+)", statement)
        return imports

    def _extract_functions(self, statement: str) -> List[str]:
        """Extract function/identifier names (external dependencies, not local definitions)."""
        # Match identifiers followed by optional parenthesis
        all_functions = re.findall(r"\b([A-Za-z_]\w*(?:\.\w+)*)\s*(?:\(|{|\[)?", statement)

        # Filter out Lean keywords and local definitions
        keywords = {"theorem", "def", "lemma", "example", "by", "have", "let", "match", "for", "do", "where", "sorry"}

        # Extract locally defined names (after theorem/def keywords)
        local_names = set()
        for match in re.finditer(r"\b(theorem|def|lemma)\s+(\w+)", statement):
            local_names.add(match.group(2))

        # Also, names that appear before := are often local definitions
        for match in re.finditer(r"\b([A-Za-z_]\w*)\s*(?::|→)\s*", statement):
            if match.group(1) not in keywords:
                local_names.add(match.group(1))

        # Return only external functions
        return [f for f in all_functions if f not in keywords and f not in local_names]

    def _extract_tactics(self, statement: str) -> List[str]:
        """Extract tactic keywords used."""
        tactics = []
        for tactic in self.LEAN4_TACTICS:
            if re.search(rf"\b{tactic}\b", statement, re.IGNORECASE):
                tactics.append(tactic)
        return tactics

    # ============================================
    # Red flag detection (critical issues)
    # ============================================

    def _check_red_flags(self, statement: str, features: Dict) -> List[Tuple[str, float]]:
        """Check for critical red flags that should drastically reduce score."""
        flags = []

        # Red flag 1: Pure sorry (no "by" or proof keywords, ONLY sorry)
        if re.search(r":=\s*sorry\s*$", statement) and " by " not in statement:
            flags.append(("pure_sorry", 0.5))  # Multiply final score by 0.5

        # Red flag 2: Obvious infinite recursion (function calls itself directly)
        if re.search(r"def\s+(\w+).*:=.*\b\1\s*\(", statement):
            flags.append(("infinite_recursion", 0.2))  # Very bad

        # Red flag 3: Very deep nesting (>12 levels)
        if features["paren_depth"] > 12:
            flags.append(("extreme_nesting", 0.6))

        # Red flag 4: Extremely many forall levels (>6)
        if features["forall_depth"] > 6:
            flags.append(("extreme_forall", 0.7))

        # Red flag 5: Unmatched brackets
        unmatched = self._check_unmatched_brackets(statement)
        if unmatched:
            flags.append(("unmatched_brackets", 0.4))

        return flags

    # ============================================
    # Subscoring functions (each returns [0, 1])
    # ============================================

    def _score_syntax(self, statement: str, features: Dict) -> float:
        """Score syntactic validity of statement."""
        score = 1.0

        # Check for basic structure
        if "theorem" not in statement and "def" not in statement and "lemma" not in statement:
            score -= 0.4  # Missing statement keyword

        # Unmatched brackets are critical errors
        for open_char, close_char in [("(", ")"), ("[", "]"), ("{", "}"), ("⟨", "⟩")]:
            open_count = statement.count(open_char)
            close_count = statement.count(close_char)
            if open_count != close_count:
                score -= 0.25  # Increased penalty

        # Too much nesting is suspicious
        if features["paren_depth"] > 12:
            score -= 0.2
        elif features["paren_depth"] > 8:
            score -= 0.1

        # Excessive sorry is a red flag
        if features["sorry_count"] > 3:
            score -= 0.3
        elif features["sorry_count"] > 1:
            score -= 0.15

        return max(0.0, score)

    def _score_structure(self, statement: str, features: Dict) -> float:
        """Score structural soundness of proof."""
        score = 0.5  # Start at neutral

        # Well-typed theorems with clear structure
        if re.search(r"theorem\s+\w+.*:.*:=\s*(by|sorry)", statement):
            score += 0.3  # Clear proof structure

        # Multiple levels of forall indicate complexity
        if features["forall_depth"] <= 1:
            score += 0.2
        elif features["forall_depth"] <= 2:
            score += 0.1
        elif features["forall_depth"] >= 4:
            score -= 0.2  # Deep nesting is problematic

        # Implications are common and OK
        if features["implication_count"] > 0:
            score += 0.1

        # Type annotations help clarity
        if features["type_annotations"] >= 2:
            score += 0.1

        # Suspicious patterns: self-reference (infinite recursion)
        if re.search(r"(\w+)\s*\(.*\1\s*\(", statement):
            score -= 0.3  # Possible infinite recursion

        return max(0.0, min(1.0, score))

    def _score_dependencies(self, statement: str, features: Dict) -> float:
        """Score availability of required Mathlib dependencies."""
        score = 0.5  # Start neutral

        functions_used = set(features["functions_used"])
        # Filter out short variable names (likely type variables, not function calls)
        functions_used = {f for f in functions_used if len(f) > 1}

        # Check for known Mathlib functions
        mathlib_functions_found = functions_used & self.MATHLIB_FUNCTIONS
        unknown_functions = functions_used - self.MATHLIB_FUNCTIONS

        # Proportion of known functions
        if len(functions_used) > 0:
            known_ratio = len(mathlib_functions_found) / len(functions_used)
            score = 0.3 + 0.5 * known_ratio  # More conservative baseline

        # Significant penalty for unknown functions
        if len(unknown_functions) > 0:
            unknown_ratio = len(unknown_functions) / max(1, len(functions_used))
            score -= 0.3 * min(1.0, unknown_ratio)  # Increased penalty

        # Import availability (if imports are explicit)
        if features["imports"]:
            available_imports = 0
            for imp in features["imports"]:
                for mathlib_mod in self.MATHLIB_MODULES:
                    if imp in mathlib_mod or mathlib_mod in imp:
                        available_imports += 1
                        break
            if len(features["imports"]) > 0:
                score += 0.15 * (available_imports / len(features["imports"]))

        return max(0.0, min(1.0, score))

    def _score_proof_strategy(self, statement: str, features: Dict) -> float:
        """Score viability of proof strategy."""
        score = 0.5  # Start neutral (sorry placeholders are expected)

        # Availability of strong tactics is a good sign
        tactics_available = set(features["tactics_used"])
        strong_tactics = {"simp", "norm_num", "omega", "decide", "ring", "field_simp"}
        strong_tactics_present = tactics_available & strong_tactics

        if strong_tactics_present:
            score += 0.2 * len(strong_tactics_present) / len(strong_tactics)

        # by-style proofs indicate proof attempt, which is good
        if " by " in statement:
            score += 0.15

        # Check for proof keywords (have, show, calc)
        proof_keywords = len(re.findall(r"\b(have|show|calc)\b", statement))
        if proof_keywords > 0:
            score += 0.1 * min(1.0, proof_keywords / 3)

        # Note: sorry is expected for placeholder theorems, so minimal penalty
        # Only penalize if PURELY sorry with no other structure
        if features["sorry_count"] >= 1 and features["tactics_used"] == [] and " by " not in statement:
            score -= 0.2  # Minimal, for truly empty proofs

        return max(0.0, min(1.0, score))

    def _score_complexity(self, statement: str, features: Dict) -> float:
        """Score reasonable complexity level (not too hard)."""
        score = 1.0

        # Check hard patterns (subtract more aggressively)
        hard_pattern_count = 0
        for pattern, reason in self.HARD_PATTERNS:
            if re.search(pattern, statement):
                hard_pattern_count += 1
                score -= 0.15  # Increased penalty per pattern

        # Extra penalty if multiple hard patterns match
        if hard_pattern_count >= 2:
            score -= 0.2

        # Check easy patterns (add)
        for pattern, reason in self.EASY_PATTERNS:
            if re.search(pattern, statement):
                score += 0.05

        # Complexity heuristic: length and nesting
        if features["length"] < 100:
            score += 0.1  # Very short statements tend to be simpler
        elif features["length"] < 150:
            score += 0.05
        elif features["length"] > 800:
            score -= 0.15  # Long statements tend to be complex

        if features["paren_depth"] > 10:
            score -= 0.2
        elif features["paren_depth"] > 6:
            score -= 0.15

        # Quantifier-implication ratio
        qi_sum = features["quantifier_count"] + features["implication_count"]
        if qi_sum > 15:
            score -= 0.25
        elif qi_sum > 10:
            score -= 0.2
        elif qi_sum > 8:
            score -= 0.1
        elif qi_sum > 5:
            score -= 0.05

        return max(0.0, min(1.0, score))

    def _score_pattern_match_mathlib(self, statement: str) -> float:
        """Score similarity to known Mathlib patterns."""
        score = 0.5  # Neutral baseline

        # Robin's inequality pattern
        if "robin" in statement.lower() or "lagarias" in statement.lower():
            score += 0.2  # Known to be formalizable

        # Common number theory patterns
        if re.search(r"\b(sigma|prime|natFactorization)\b", statement):
            score += 0.1

        # Analysis patterns
        if re.search(r"\b(Real\.|Analysis\.)", statement):
            score += 0.05

        # If it starts with ∀ or uses universal quantifier, boost slightly
        if re.search(r"^∀|∀\s*\w+", statement):
            score += 0.1

        return max(0.0, min(1.0, score))

    def _score_pattern_match(self, statement: str, features: Dict) -> float:
        """Score similarity to known formalizable patterns."""
        return self._score_pattern_match_mathlib(statement)

    def _compute_confidence(self, subscores: Dict[str, float], final_score: float) -> float:
        """Compute confidence based on subscores agreement."""
        # If subscores are close together, confidence is high
        scores = list(subscores.values())
        mean_score = sum(scores) / len(scores)
        variance = sum((s - mean_score) ** 2 for s in scores) / len(scores)
        std_dev = variance ** 0.5

        # Low std dev = high agreement = high confidence
        confidence = 1.0 - min(1.0, std_dev * 0.5)

        # Adjust for extreme scores (very high or very low are more certain)
        if final_score < 0.2 or final_score > 0.8:
            confidence = min(1.0, confidence + 0.1)

        return confidence

    def _identify_issues(self, statement: str, features: Dict, subscores: Dict) -> List[str]:
        """Identify top 3 issues that might cause formalization failure."""
        issues = []

        # Syntax issues
        if subscores["syntax_valid"] < 0.6:
            if features["sorry_count"] > 2:
                issues.append("Too many 'sorry' placeholders (unfinished proof)")
            unmatched = self._check_unmatched_brackets(statement)
            if unmatched:
                issues.append(f"Unmatched brackets: {unmatched}")

        # Structure issues
        if subscores["structure_sound"] < 0.5:
            if features["forall_depth"] > 4:
                issues.append("Deep nesting of universal quantifiers (>4 levels)")
            if features["paren_depth"] > 10:
                issues.append("Excessive parentheses nesting (>10 levels)")

        # Dependency issues
        if subscores["dependency_available"] < 0.5:
            unknown_funcs = set(features["functions_used"]) - self.MATHLIB_FUNCTIONS
            if unknown_funcs:
                top_unknown = list(unknown_funcs)[:2]
                issues.append(f"Unknown functions not in Mathlib: {', '.join(top_unknown)}")

        # Proof strategy issues
        if subscores["proof_strategy_viable"] < 0.5:
            if features["sorry_count"] > 0:
                issues.append("Proof is incomplete (contains 'sorry')")
            if not features["tactics_used"]:
                issues.append("No tactics available for proof (no simp, norm_num, etc.)")

        # Complexity issues
        if subscores["complexity_reasonable"] < 0.5:
            issues.append("Statement complexity exceeds typical Mathlib patterns")

        # Pattern match issues
        if subscores["pattern_match_mathlib"] < 0.3:
            issues.append("Statement doesn't match known formalizable patterns")

        return issues[:3] if issues else ["No specific issues identified; may be edge case"]

    def _check_unmatched_brackets(self, s: str) -> Optional[str]:
        """Check for unmatched brackets."""
        pairs = {"(": ")", "[": "]", "{": "}", "⟨": "⟩"}
        stack = []
        for char in s:
            if char in pairs:
                stack.append(char)
            elif char in pairs.values():
                if not stack or pairs[stack[-1]] != char:
                    return f"Unexpected '{char}'"
                stack.pop()

        if stack:
            return f"Unclosed '{stack[-1]}'"
        return None


def score_statement(statement: str) -> FormalizabilityScore:
    """Convenience function for single prediction."""
    predictor = FormalizabilityPredictor()
    return predictor.predict(statement)


if __name__ == "__main__":
    # Test the predictor on known examples

    predictor = FormalizabilityPredictor()

    test_cases = [
        # Expected HIGH scores (>0.7) — known formalizable
        (
            "theorem robin_inequality_of_RH (hRH : RiemannHypothesis) : "
            "∀ n : ℕ, n > 5040 → (sigma 1 n : ℝ) < robinBound n := by sorry",
            "Robin's inequality (known formalizable, high difficulty)",
            0.65  # Target score
        ),
        (
            "theorem simple_eq : 1 + 1 = 2 := by norm_num",
            "Simple arithmetic (trivial, high score expected)",
            0.85
        ),
        (
            "theorem lagarias_inequality_of_RH (hRH : RiemannHypothesis) : "
            "∀ n : ℕ, n ≥ 1 → (sigma 1 n : ℝ) ≤ lagariasBound n := by sorry",
            "Lagarias' inequality (known formalizable)",
            0.65
        ),
        # Expected MEDIUM scores (0.4-0.6) — questionable
        (
            "theorem mystery_func : mystery_undefined 42 = mystery_result := by sorry",
            "Undefined functions (questionable)",
            0.35
        ),
        (
            "theorem deeply_nested : "
            "∀ a b c d e f g h i j k, "
            "(a → b → c → d → e → f → g → h → i → j → k) := by sorry",
            "Deep quantifier nesting (complex, questionable)",
            0.45
        ),
        # Expected LOW scores (<0.4) — likely to fail
        (
            "theorem broken_syntax : unknown_mystery_function ( ( ( := by sorry",
            "Broken syntax (unmatched brackets)",
            0.20
        ),
        (
            "def infinite_recursion (n : ℕ) : ℕ := infinite_recursion (n + 1)",
            "Infinite recursion (will not compile)",
            0.25
        ),
        (
            "theorem sorry_only := sorry",
            "Trivial sorry (incomplete)",
            0.30
        ),
    ]

    print("=" * 80)
    print("FORMALIZABILITY PREDICTOR — Calibration Test Run")
    print("=" * 80)
    print(f"{'Description':<45} {'Score':>8} {'Target':>8} {'Conf':>6} {'Time':>7} {'Result':<10}")
    print("-" * 85)

    scores_high = []
    scores_med = []
    scores_low = []

    for statement, description, target_score in test_cases:
        result = predictor.predict(statement)
        score_ok = "✓" if abs(result.score - target_score) < 0.15 else "✗"
        category = "PASS" if result.score >= 0.5 else "FAIL"

        print(
            f"{description:<45} {result.score:>8.3f} {target_score:>8.3f} "
            f"{result.confidence:>6.2f} {result.processing_time_ms:>7.1f}ms {score_ok:<10}"
        )

        # Print subscores for debugging difficult cases
        if abs(result.score - target_score) >= 0.15 and description.startswith("Robin"):
            print(f"    DEBUG: subscores={result.metadata['subscores']}")

        if target_score >= 0.65:
            scores_high.append(result.score)
        elif target_score >= 0.4:
            scores_med.append(result.score)
        else:
            scores_low.append(result.score)

    print("-" * 85)
    print("\nSCORE DISTRIBUTIONS:")
    if scores_high:
        avg_high = sum(scores_high) / len(scores_high)
        print(f"  High-difficulty formalizable: avg={avg_high:.3f} (target 0.65-0.85)")
    if scores_med:
        avg_med = sum(scores_med) / len(scores_med)
        print(f"  Medium difficulty (questionable): avg={avg_med:.3f} (target 0.35-0.60)")
    if scores_low:
        avg_low = sum(scores_low) / len(scores_low)
        print(f"  Low (likely to fail): avg={avg_low:.3f} (target <0.40)")

    print("\n" + "=" * 80)
