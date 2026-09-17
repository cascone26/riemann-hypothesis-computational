#!/usr/bin/env python3
"""
Simplified K5 generator: copy an existing proven file and adapt for new degree.

The key insight is that K5 coefficient files differ ONLY in:
1. Theorem name
2. The goal formula (theorem conclusion)
3. The coefficient degree being extracted (hc19 vs hc18 etc)
4. The antidiagonal pairs used in simp

We generate by:
1. Loading CrossK5FifthCoeff (degree 15, proven)
2. Copying structure to new file
3. Adapting minimal set of things based on target degree
4. Manually constructing goal formula from antidiagonal pairs
"""

import sys
import re
from pathlib import Path
from typing import List, Tuple

ORDINAL_NAMES = {
    19: "LeadingCoeff",
    18: "SecondCoeff",
    17: "ThirdCoeff",
    16: "FourthCoeff",
    15: "FifthCoeff",
    14: "SixthCoeff",
    13: "SeventhCoeff",
    12: "EighthCoeff",
    11: "NinthCoeff",
    10: "TenthCoeff",
    9:  "EleventhCoeff",
    8:  "TwelfthCoeff",
    7:  "ThirteenthCoeff",
    6:  "FourteenthCoeff",
    5:  "FifteenthCoeff",
    4:  "SixteenthCoeff",
    3:  "SeventeenthCoeff",
    2:  "EighteenthCoeff",
    1:  "NineteenthCoeff",
}

def sin_lean(i: int) -> str:
    """Lean syntax for sqrt(2)*sin(L*z_i/2)"""
    return f"(((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z{i} / 2))"

def csin_lean(i: int) -> str:
    """Lean syntax for conj(sqrt(2)*sin(L*z_i/2))"""
    return f"(starRingEnd ℂ (((Real.sqrt 2 : ℝ) : ℂ) * Complex.sin (L * z{i} / 2)))"

def coeff_0_lean(i: int) -> str:
    """coeff_0 of cubic_i: (p+q)*sin*z*zc^2 + (r+s)*conj*z^2*zc"""
    p, q, r, s = f"p{i}", f"q{i}", f"r{i}", f"s{i}"
    return f"(({p} + {q}) * {sin_lean(i)} * (z{i} * (starRingEnd ℂ z{i}) ^ 2) + ({r} + {s}) * {csin_lean(i)} * (z{i} ^ 2 * (starRingEnd ℂ z{i})))"

def coeff_2_lean(i: int) -> str:
    """coeff_2 of cubic_i: -(p+q)*sin*z - (r+s)*conj*zc"""
    p, q, r, s = f"p{i}", f"q{i}", f"r{i}", f"s{i}"
    return f"(-(({p} + {q}) * {sin_lean(i)} * z{i} + ({r} + {s}) * {csin_lean(i)} * (starRingEnd ℂ z{i})))"

def coeff_3_lean(i: int) -> str:
    """coeff_3 of cubic_i: (p-q)*sin + (r-s)*conj*... wait, this one is different"""
    # Looking at the actual files, coeff_3 is just (p-q)*sin + (r-s)*conj*sin
    # No z dependencies!
    p, q, r, s = f"p{i}", f"q{i}", f"r{i}", f"s{i}"
    return f"(({p} - {q}) * {sin_lean(i)} + ({r} - {s}) * {csin_lean(i)})"

def others(i: int) -> List[int]:
    """Return indices other than i"""
    return [j for j in range(1, 6) if j != i]

def antidiagonal_pairs(d: int) -> List[Tuple[int, int]]:
    """Valid (a,b) pairs for degree d where a+b=d, a in {0,1,2,3}, b even in [0,16]"""
    pairs = []
    for a in range(4):
        b = d - a
        if 0 <= b <= 16 and b % 2 == 0:
            pairs.append((a, b))
    return pairs

def magnitude_poly_description(b: int, indices: List[int]) -> str:
    """
    Describe what the degree-b coefficient of quadruple product looks like.
    indices: the 4 indices (other than the current block's index).

    For simplicity, we return a symbolic description. In the actual Lean code,
    this will appear as one of the terms from the expanded quartic products.
    """
    # This is tricky because we need to enumerate all the B and A combinations
    # that produce degree b. For now, return a symbolic name.
    zlist = "_".join(f"z{i}" for i in indices)
    return f"(magnitude_{b}_of_{zlist})"

def generate_goal_formula(d: int, include_block_breakdown: bool = False) -> str:
    """
    Generate the goal formula for degree d.
    Returns Lean code for: (sum of 5 block contributions) = 0
    """
    pairs = antidiagonal_pairs(d)
    if not pairs:
        return f"(* NO VALID PAIRS FOR DEGREE {d} *) = 0"

    block_strs = []
    for i in range(1, 6):
        js, ks, ls, ms = others(i)
        terms = []

        for a, b in pairs:
            if a == 0:
                # coeff_0 term: (p+q)*sin*z*zc^2 + (r+s)*conj*z^2*zc
                coeff_part = coeff_0_lean(i)
            elif a == 1:
                # coeff_1: (q-p)*zc^2 + (s-r)*z^2
                # This is a coefficient of X in the cubic
                p, q, r, s = f"p{i}", f"q{i}", f"r{i}", f"s{i}"
                coeff_part = f"(({q} - {p}) * (starRingEnd ℂ z{i}) ^ 2 + ({s} - {r}) * z{i} ^ 2)"
            elif a == 2:
                # coeff_2 term: -(p+q)*sin*z - (r+s)*conj*zc
                coeff_part = coeff_2_lean(i)
            elif a == 3:
                # coeff_3 term: (p-q)*sin + (r-s)*conj
                coeff_part = coeff_3_lean(i)

            mag_desc = magnitude_poly_description(b, [js, ks, ls, ms])
            terms.append(f"{coeff_part} * {mag_desc}")

        block_str = "(" + " + ".join(terms) + ")"
        block_strs.append(block_str)

    goal = "\n      + ".join(block_strs) + " = 0"
    return goal

def load_template(source_file: Path) -> str:
    """Load existing K5 file as template"""
    with open(source_file) as f:
        return f.read()

def extract_goal_from_file(content: str) -> Tuple[int, int, str]:
    """
    Extract the goal formula from existing Lean file.
    Returns (start_pos, end_pos, goal_text) for the conclusion part.
    """
    # The goal formula is between ") :" and " := by"
    # Example: ") : ((FORMULA)) = 0 := by"
    pattern = r"\) : (.*?) := by"
    m = re.search(pattern, content, re.DOTALL)
    if m:
        return m.start(1), m.end(1), m.group(1)
    return None

def adapt_for_degree(template: str, old_degree: int, new_degree: int) -> str:
    """
    Adapt template file for new degree.

    This is INCOMPLETE because we don't yet have the correct goal formula derivation.
    For now, we just rename and mark with TODO.
    """
    old_ordinal = ORDINAL_NAMES.get(old_degree)
    new_ordinal = ORDINAL_NAMES.get(new_degree)

    if not old_ordinal or not new_ordinal:
        raise ValueError(f"Invalid degrees: {old_degree}, {new_degree}")

    result = template

    # 1. Rename theorem
    result = result.replace(f"crossK5{old_ordinal}", f"crossK5{new_ordinal}")

    # 2. Update comment
    result = re.sub(
        rf"# K=5.*?degree-{old_degree}.*?\)",
        f"# K=5 coefficient extraction: the degree-{new_degree} coefficient vanishes ({new_ordinal})",
        result,
        flags=re.DOTALL
    )

    # 3. Update hc degree
    result = result.replace(f"hc{old_degree}", f"hc{new_degree}")

    return result

def main():
    if len(sys.argv) < 2:
        print("Usage: gen_k5_from_template.py <target_degree> [source] [output]", file=sys.stderr)
        sys.exit(1)

    target_degree = int(sys.argv[1])
    source_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path(__file__).parent / "CrossK5FifthCoeff.lean"
    output_path = Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).parent / f"CrossK5{ORDINAL_NAMES[target_degree]}.lean"

    if not source_path.exists():
        print(f"Error: source file not found: {source_path}", file=sys.stderr)
        return 1

    pairs = antidiagonal_pairs(target_degree)
    print(f"Generating degree {target_degree}", file=sys.stderr)
    print(f"  Antidiagonal pairs: {pairs}", file=sys.stderr)
    print(f"  Source: {source_path}", file=sys.stderr)
    print(f"  Output: {output_path}", file=sys.stderr)

    # Extract source degree from file
    with open(source_path) as f:
        source_content = f.read()
    
    source_match = re.search(r"degree-(\d+)", source_content)
    if not source_match:
        print("Error: could not determine source degree", file=sys.stderr)
        return 1
    
    source_degree = int(source_match.group(1))
    print(f"  Source degree: {source_degree}", file=sys.stderr)

    try:
        adapted = adapt_for_degree(source_content, source_degree, target_degree)
        with open(output_path, 'w') as f:
            f.write(adapted)
        print(f"✓ Generated {output_path}")
        return 0
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
