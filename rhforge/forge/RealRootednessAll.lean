import HankelPowerSumSOS
import WeightedHankelInterval
import MultisetHankel
import PolynomialHankel
import ConjPairNeg
import ComplexHankel
import ComplexHankelReal
import HankelConverseCore
import RealLinearSurj
import HankelPairConverse
import ConjClosedProdReal
import ProdEvalNeZero
import RealPolyCoeffVector
import HankelFullConverse
import ComplexHankelNonnegReal
import ConjClosedDecomp
import HankelIff
import IntervalConverseCore
import RealRootIsolation
import IntervalConverse
import IntervalWeightedNonneg
import IntervalIff
import PalindromicAllRoots
import QpolyDef
import PalindromicQpolyRoot
import DicksonMapConj
import QpolyMapConj
import QpolyHankelIff

/-!
# Real-Rootedness + On-Circle module — aggregator (§1790–1817, 28 bricks)

Hermite real-rootedness criterion (both iffs §1806), interval localization iff (§1811), Chebyshev on-circle
criterion (§1781–1789), and their wiring (§1812–1817).  Compiling this verifies the whole DAG builds coherently.
All axiom-clean (⊆{propext,Classical.choice,Quot.sound}).  Mathlib-absent, PR-worthy.
-/
