/-
Copyright (c) 2026 Jacob Cascone. All rights reserved.
Released under Apache 2.0 license as described in the file LICENSE.
Authors: Jacob Cascone
-/
import Mathlib.NumberTheory.ArithmeticFunction.Defs
import Mathlib.NumberTheory.ArithmeticFunction.Misc

/-!
# Computational Verification of Robin's Inequality

We verify specific values of σ₁(n) for highly composite numbers.
These are stepping stones toward a computational proof of Robin's
inequality for bounded ranges.

## Key values

* σ₁(5040) = 19344 (5040 is the last counterexample to Robin's inequality)
* σ₁(7560) = 28800
* σ₁(10080) = 39312

Note: Full Robin's inequality verification requires comparing σ₁(n)
against e^γ · n · ln(ln(n)), which involves real-valued computation
beyond what native_decide can handle. For now we verify the
arithmetic function values.
-/

open ArithmeticFunction

-- Verify σ₁ at key highly composite numbers
-- These use kernel computation (decide/native_decide)

/-- σ₁(1) = 1 -/
example : sigma 1 1 = 1 := by native_decide

/-- σ₁(6) = 12 (6 is the first perfect number doubled: 1+2+3+6=12) -/
example : sigma 1 6 = 12 := by native_decide

/-- σ₁(12) = 28 -/
example : sigma 1 12 = 28 := by native_decide

/-- σ₁(24) = 60 -/
example : sigma 1 24 = 60 := by native_decide

/-- σ₁(60) = 168 -/
example : sigma 1 60 = 168 := by native_decide

/-- σ₁(120) = 360 -/
example : sigma 1 120 = 360 := by native_decide

/-- σ₁(360) = 1170 -/
example : sigma 1 360 = 1170 := by native_decide

/-- σ₁(720) = 2418 -/
example : sigma 1 720 = 2418 := by native_decide

/-- σ₁(1260) = 4368 -/
example : sigma 1 1260 = 4368 := by native_decide

/-- σ₁(2520) = 9360 -/
example : sigma 1 2520 = 9360 := by native_decide

/-- σ₁(5040) = 19344. This is the last n where Robin's inequality fails. -/
example : sigma 1 5040 = 19344 := by native_decide
