#!/usr/bin/env python3
"""
Generator for K=5 coefficient extraction files (CrossK5<Ordinal>Coeff.lean).

Extension of build_k4reduced_generic.py logic to 5 blocks.
For a given target coefficient degree d, derives antidiagonal pairs and generates
the full Lean file by transforming an existing K5 coefficient file.

KEY INSIGHT: K5 files have identical structure - only the goal formula and
coefficient degree vary. We generate by:
1. Loading a template (e.g., CrossK5ThirdCoeff)
2. Deriving the correct antidiagonal pairs for the target degree
3. Transforming the theorem statement's goal formula
4. Replacing the coefficient extraction degree number
"""

import sys
import re
from pathlib import Path

# === ORDINAL NAMES ===
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

def derive_antidiagonal_pairs(d):
    """
    Derive valid antidiagonal pairs (a,b) for degree d in the coefficient extraction.

    For K=5:
    - cubic_i has degree 3 (coefficients at indices 0, 1, 2, 3)
    - quadruple product Q_i has degree 16 (from 4×4=degree 16 quartic product)
    - valid b indices: {0, 2, 4, 6, 8, 10, 12, 14, 16} (even degrees only)

    Returns list of (a, b) pairs where a ∈ {0,1,2,3}, b ∈ {even in [0,16]}, a+b=d
    """
    pairs = []
    for a in range(4):
        b = d - a
        if 0 <= b <= 16 and b % 2 == 0:
            pairs.append((a, b))
    return pairs

def extract_section_from_file(filepath, start_pattern, end_pattern):
    """
    Extract a section of file content between two regex patterns.
    Returns the matched region including start/end markers.
    """
    with open(filepath, 'r') as f:
        content = f.read()

    start_match = re.search(start_pattern, content)
    end_match = re.search(end_pattern, content)

    if not start_match or not end_match:
        return None

    return content[start_match.start():end_match.end()]

def transform_theorem_statement(template, old_degree, new_degree, ordinal_name):
    """
    Transform the theorem statement from old_degree to new_degree.
    This involves:
    1. Renaming the theorem
    2. Updating the goal formula (THIS IS PLACEHOLDER - REAL VERSION NEEDS SYMBOLIC DERIVATION)
    3. Keeping everything else the same
    """
    # Extract theorem statement from template
    theorem_pattern = rf"theorem cross.*?:=\s*by"

    # For now, we'll use a simpler approach: copy from source and manually adjust
    # In a production version, we'd symbolically derive the correct goal formula
    # But that requires understanding the full antidiagonal structure for each degree

    return template  # Placeholder return

def generate_file_by_transformation(source_degree, target_degree, source_file, output_file):
    """
    Generate a new K5 coefficient file by copying an existing one and transforming it.

    Strategy: The files differ ONLY in:
    1. Theorem name
    2. Goal formula (antidiagonal structure)
    3. Coefficient degree extracted (hc19 vs hc18 vs etc)
    4. Antidiagonal pairs in simp_only step

    For now, we implement a simplified version that:
    - Validates the structure
    - Copies the file
    - Marks with TODO for manual goal-formula derivation

    PROPER IMPLEMENTATION: Use sympy to derive the goal formula symbolically,
    similar to build_k4reduced_generic.py, extended to 5 blocks.
    """
    if not source_file.exists():
        raise FileNotFoundError(f"Source file not found: {source_file}")

    with open(source_file, 'r') as f:
        template = f.read()

    target_ordinal = ORDINAL_NAMES.get(target_degree)
    if not target_ordinal:
        raise ValueError(f"Unknown degree: {target_degree}")

    source_ordinal = ORDINAL_NAMES.get(source_degree)
    if not source_ordinal:
        raise ValueError(f"Unknown source degree: {source_degree}")

    # Basic transformations (safe ones that don't change logic)
    transformed = template

    # 1. Rename theorem
    old_theorem_name = f"crossK5{source_ordinal}"
    new_theorem_name = f"crossK5{target_ordinal}"
    transformed = transformed.replace(old_theorem_name, new_theorem_name)

    # 2. Update doc comment degree
    transformed = re.sub(
        rf"# K=5.*?degree-{source_degree}.*",
        f"# K=5 coefficient extraction: the degree-{target_degree} coefficient vanishes ({target_ordinal})",
        transformed,
        flags=re.DOTALL
    )

    # 3. Update coefficient extraction degree (hc19 -> hc18 etc)
    transformed = transformed.replace(f"hc{source_degree}", f"hc{target_degree}")

    # TODO: Update goal formula (requires symbolic derivation)
    # This is the hard part and needs proper antidiagonal derivation
    # For now, we'll mark it as needing manual review
    transformed += f"\n\n(* TODO: Goal formula for degree {target_degree} needs symbolic verification *)\n"

    with open(output_file, 'w') as f:
        f.write(transformed)

    return transformed

def main():
    if len(sys.argv) < 2:
        print("Usage: gen_k5_coeff.py <target_degree> [source_degree] [source_file] [output_file]",
              file=sys.stderr)
        print("  target_degree: desired coefficient degree (19, 18, 17, ...)", file=sys.stderr)
        print("  source_degree: source file degree (default: closest available)", file=sys.stderr)
        print("  source_file: template file path (default: CrossK5ThirdCoeff.lean)", file=sys.stderr)
        print("  output_file: output path (default: ./CrossK5<Ordinal>Coeff.lean)", file=sys.stderr)
        sys.exit(1)

    target_degree = int(sys.argv[1])

    # Infer source from target if not provided
    # Use degree 17 (ThirdCoeff) as default source - it's proven and well-structured
    default_source_degree = 17
    source_degree = int(sys.argv[2]) if len(sys.argv) > 2 else default_source_degree

    work_dir = Path(__file__).parent
    default_source_file = work_dir / f"CrossK5{ORDINAL_NAMES.get(source_degree)}.lean"
    source_file = Path(sys.argv[3]) if len(sys.argv) > 3 else default_source_file

    target_ordinal = ORDINAL_NAMES.get(target_degree)
    if not target_ordinal:
        print(f"Error: no ordinal name for degree {target_degree}", file=sys.stderr)
        sys.exit(1)

    default_output_file = work_dir / f"CrossK5{target_ordinal}.lean"
    output_file = Path(sys.argv[4]) if len(sys.argv) > 4 else default_output_file

    print(f"Generating K=5 coefficient file for degree {target_degree}", file=sys.stderr)
    print(f"  Source: degree {source_degree} ({source_file})", file=sys.stderr)
    print(f"  Target: degree {target_degree} ({output_file})", file=sys.stderr)

    try:
        pairs = derive_antidiagonal_pairs(target_degree)
        print(f"  Antidiagonal pairs: {pairs}", file=sys.stderr)

        generate_file_by_transformation(source_degree, target_degree, source_file, output_file)
        print(f"✓ Generated {output_file}")
        return 0
    except Exception as e:
        print(f"✗ Error: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())
