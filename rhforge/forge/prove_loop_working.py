#!/usr/bin/env python3
"""
Automated Lean proof loop (v2): read a file with `by sorry`, generate tactic candidates
from both a rule-based ladder AND the local AI gateway (vllm on http://localhost:8000),
splice and verify until one compiles with exit 0.

Falls back gracefully if the AI gateway is unreachable.
"""

import re
import sys
import subprocess
import tempfile
import requests
import json
from pathlib import Path


# Tactics to try in order (rule-based ladder)
TACTICS = [
    "rfl",
    "trivial",
    "exact rfl",
    "decide",
    "norm_num",
    "simp",
    "simp_all",
    "ring",
    "ring_nf",
    "omega",
    "linarith",
    "nlinarith",
]

# AI gateway endpoints to try (in order)
AI_GATEWAY_ENDPOINTS = [
    ("http://localhost:8000/v1/chat/completions", "qwen3-32b-awq"),  # vllm primary
    ("http://localhost:7792/v1/chat/completions", "gpt-3.5-turbo"),  # fallback
]
AI_GATEWAY_TIMEOUT = 90  # 32B reasoning model needs headroom  # seconds

# Detected working endpoint (cached)
_detected_endpoint = None
_detected_model = None


def find_all_sorries(content: str) -> list[tuple[int, int]]:
    """
    Find all occurrences of `by sorry` in the content.
    Returns list of (start, end) positions.
    """
    sorries = []
    for match in re.finditer(r"by\s+sorry\b", content):
        sorries.append((match.start(), match.end()))
    return sorries


def replace_sorry_at_index(content: str, sorry_index: int, tactic: str) -> str:
    """
    Replace the sorry_index-th occurrence of `by sorry` with `by <tactic>`.
    """
    sorries = find_all_sorries(content)
    if sorry_index >= len(sorries):
        return content

    start, end = sorries[sorry_index]
    before = content[:start]
    after = content[end:]
    return before + f"by {tactic}" + after


def verify_file(file_path: str) -> tuple[bool, str]:
    """
    Run `lake env lean` on the file. Return (success, stdout+stderr_output).
    """
    try:
        result = subprocess.run(
            ["bash", "-c", f"export PATH=\"$HOME/.elan/bin:$PATH\"; cd ~/rhforge/mathlib4; lake env lean '{file_path}' 2>&1"],
            capture_output=True,
            timeout=30,
            text=True
        )
        return (result.returncode == 0, result.stdout + result.stderr)
    except subprocess.TimeoutExpired:
        return (False, "Timeout")


FORGE_DIR = Path.home() / "rhforge" / "forge"


def gather_import_lemmas(content: str, budget: int = 1500) -> str:
    """Read the imported FORGE modules (not Mathlib) and list their theorem/lemma signatures,
    so the AI knows which project lemmas it may apply. Budget-limited for the 2048-token window."""
    sigs = []
    for m in re.finditer(r'^import\s+([A-Za-z_]\w*)\s*$', content, flags=re.M):
        mod = m.group(1)
        f = FORGE_DIR / f"{mod}.lean"
        if not f.exists():
            continue  # Mathlib / stdlib import -> skip
        try:
            txt = f.read_text()
        except Exception:
            continue
        for tm in re.finditer(r'^(?:theorem|lemma)\s+(\w+)', txt, flags=re.M):
            name = tm.group(1)
            tail = txt[tm.start():]
            cut = re.search(r':=|\n\s*by\b|:= by', tail)
            sig = tail[:cut.start()] if cut else tail[:200]
            sig = re.sub(r'\s+', ' ', sig).strip()[:200]
            sigs.append(sig)
    if not sigs:
        return ""
    body = "\n".join(f"- {x}" for x in sigs)
    return ("Available project lemmas from the imports (APPLY THESE by name if they match):\n"
            + body)[:budget]


def extract_goal_context(content: str, sorry_index: int) -> str:
    """
    Extract the theorem/example context around the sorry.
    Returns a string with the theorem signature and some surrounding context.
    """
    sorries = find_all_sorries(content)
    if sorry_index >= len(sorries):
        return ""

    start, _ = sorries[sorry_index]
    before = content[:start]

    # Find the last `theorem`, `example`, `lemma` keyword
    lines = before.split('\n')
    for i in range(len(lines) - 1, -1, -1):
        if re.search(r'\b(theorem|example|lemma)\b', lines[i]):
            # Grab this line and the next few (up to `by sorry`)
            snippet = '\n'.join(lines[i:])[:500]
            lemctx = gather_import_lemmas(content)
            return (lemctx + "\n\nGOAL:\n" + snippet) if lemctx else snippet

    return ""


def detect_working_endpoint():
    """
    Test all endpoints and find the first one that works.
    Returns (url, model) tuple or (None, None) if none work.
    """
    global _detected_endpoint, _detected_model

    if _detected_endpoint is not None:
        return (_detected_endpoint, _detected_model)

    for url, model in AI_GATEWAY_ENDPOINTS:
        try:
            payload = {
                "model": model,
                "messages": [{"role": "user", "content": "test"}],
                "temperature": 0.7,
                "max_tokens": 50
            }
            response = requests.post(
                url,
                json=payload,
                timeout=3
            )
            if response.status_code == 200:
                _detected_endpoint = url
                _detected_model = model
                return (url, model)
        except:
            continue

    return (None, None)


def parse_ai_response(text: str) -> list[str]:
    """
    Parse tactic suggestions from AI response, handling thinking blocks and formatting.
    """
    # Remove thinking blocks (handle both complete and incomplete tags)
    text = re.sub(r'<think>.*?</think>', '', text, flags=re.DOTALL)
    text = re.sub(r'<think>.*$', '', text, flags=re.DOTALL)  # Handle unclosed think tags

    tactics = []
    for line in text.split('\n'):
        line = line.strip()

        # Skip empty lines and comments
        if not line or line.startswith('#') or line.startswith('--') or line.startswith('```'):
            continue

        # Remove common list prefixes
        if line.startswith('- '):
            line = line[2:].strip()
        elif line.startswith('* '):
            line = line[2:].strip()
        elif re.match(r'^\d+\.\s+', line):
            line = re.sub(r'^\d+\.\s+', '', line).strip()

        # Skip obvious non-tactics
        if line.startswith('**') or line.startswith('Here') or line.startswith('The') or line.startswith('Given') or 'e.g.' in line or len(line) > 100:
            continue

        # Accept any line starting with a lowercase tactic keyword; Lean lemma args may be
        # UpperCamelCase, contain dots/brackets/parens, so DON'T over-restrict the arg chars.
        if line and re.match(r'^[a-z][a-zA-Z_0-9]*', line) and not line.endswith(':') \
                and '=' not in line.split(' ')[0]:
            tactics.append(line)

    return tactics[:8]


def call_ai_gateway(goal_context: str) -> list[str]:
    """
    Call the AI gateway (OpenAI-compatible) to generate tactic suggestions.
    Returns a list of suggested tactics, or empty list if unreachable.
    """
    url, model = detect_working_endpoint()

    if url is None or model is None:
        return []

    try:
        # Specialized prompt for Lean 4 tactics - very explicit
        prompt = f"""/no_think
You are a Lean 4 + Mathlib proof assistant. Below is a theorem whose proof is by sorry.
Output ONLY candidate proof tactics, one COMPLETE tactic per line (for example: exact foo bar,
or simpa using baz, or rw [h] then ring). Use lemmas visible from the imports. No prose.

Lean code:
{goal_context}

Candidate tactics (one complete tactic per line):"""

        payload = {
            "model": model,
            "messages": [
                {"role": "user", "content": prompt}
            ],
            "temperature": 0.2,
            "max_tokens": 512
        }

        response = requests.post(
            url,
            json=payload,
            timeout=AI_GATEWAY_TIMEOUT
        )

        if response.status_code == 200:
            data = response.json()
            text = data.get("choices", [{}])[0].get("message", {}).get("content", "")
            return parse_ai_response(text)
        else:
            return []

    except (requests.RequestException, requests.Timeout, json.JSONDecodeError, KeyError, Exception):
        return []


def prove_file(file_path: str, use_ai_gateway: bool = True) -> tuple[str, bool, list[str], bool]:
    """
    Attempt to prove all theorems with `by sorry` in the file.
    Returns (final_content, overall_success, proven_tactics, gateway_reachable).
    """
    with open(file_path, 'r') as f:
        original_content = f.read()

    # Find all sorries
    sorries = find_all_sorries(original_content)
    if not sorries:
        print(f"No `by sorry` found in {file_path}.")
        return (original_content, True, [], None)

    print(f"Found {len(sorries)} theorem(s) with `by sorry`\n")

    content = original_content
    proven_tactics = []
    gateway_reachable = None  # Will be set on first attempt

    for i, (start, end) in enumerate(sorries):
        print(f"[Theorem {i+1}/{len(sorries)}]")

        # Extract theorem name (simple heuristic)
        before_sorry = content[:start]
        lines_before = before_sorry.split('\n')
        theorem_line = lines_before[-1] if lines_before else ""
        theorem_name = re.search(r"(theorem|example|lemma)\s+(\w+)?", theorem_line)
        name_str = theorem_name.group(2) if theorem_name and theorem_name.group(2) else "unnamed"
        print(f"  Theorem: {name_str}")

        # Build candidate list: ladder tactics first, then AI suggestions
        candidates = list(TACTICS)

        if use_ai_gateway:
            goal_context = extract_goal_context(content, i)
            if goal_context:
                print(f"  Querying AI gateway...", end=" ", flush=True)
                ai_tactics = call_ai_gateway(goal_context)
                if ai_tactics is not None:
                    if gateway_reachable is None:
                        gateway_reachable = len(ai_tactics) > 0
                    if ai_tactics:
                        print(f"got {len(ai_tactics)} suggestions: {', '.join(ai_tactics[:3])}")
                        # Expand: a bare lemma-term suggestion needs `exact`/`apply`; try both forms.
                        _kw = ("exact","apply","refine","simp","rw","simpa","rfl","ring","linarith",
                               "nlinarith","omega","constructor","intro","convert","norm_num","decide",
                               "unfold","gcongr","positivity","field_simp","aesop","trivial","assumption")
                        expanded = []
                        for c in ai_tactics:
                            expanded.append(c)
                            if not c.split(" ")[0] in _kw:
                                expanded.append("exact " + c)
                                expanded.append("apply " + c)
                        # Add AI suggestions to the front (try them early)
                        candidates = expanded + candidates
                    else:
                        print("(no response)")
                else:
                    if gateway_reachable is None:
                        gateway_reachable = False
                    print("(unreachable)")

        if gateway_reachable is None:
            gateway_reachable = False

        found = False
        for tactic in candidates:
            tactic = tactic.strip()
            if not tactic:
                continue

            print(f"    Trying: {tactic}...", end=" ", flush=True)

            # Try this tactic
            test_content = replace_sorry_at_index(content, i, tactic)

            with open(file_path, 'w') as f:
                f.write(test_content)

            success, error_msg = verify_file(file_path)
            if success:
                print("✓ SUCCESS")
                content = test_content
                proven_tactics.append(f"{name_str}: {tactic}")
                found = True
                break
            else:
                print("✗", end=" ", flush=True)

        if not found:
            print(f"\n    → Could not prove. Leaving sorry.\n")
            # Restore the sorry for this one
            content = replace_sorry_at_index(content, i, "sorry")
            with open(file_path, 'w') as f:
                f.write(content)
        else:
            print()

    # Final verification
    with open(file_path, 'r') as f:
        final_content = f.read()

    success, _ = verify_file(file_path)

    return (final_content, success, proven_tactics, gateway_reachable)


def count_sorrys(content: str) -> int:
    """Count remaining sorrys."""
    return len(find_all_sorries(content))


def main():
    if len(sys.argv) < 2:
        print("Usage: prove_loop2.py <lean_file> [--no-ai-gateway]")
        sys.exit(1)

    file_path = sys.argv[1]
    use_ai = "--no-ai-gateway" not in sys.argv
    abs_path = str(Path(file_path).resolve())

    if not Path(abs_path).exists():
        print(f"File not found: {abs_path}")
        sys.exit(1)

    print(f"Starting proof loop (v2) on {abs_path}")
    if use_ai:
        url, model = detect_working_endpoint()
        if url and model:
            print(f"AI gateway: {url} (model: {model})")
        else:
            print(f"AI gateway: DISABLED (no endpoint found)")
            use_ai = False
    else:
        print(f"AI gateway: DISABLED")
    print("=" * 70 + "\n")

    final_content, success, proven_tactics, gateway_reachable = prove_file(abs_path, use_ai_gateway=use_ai)

    print("\n" + "=" * 70)
    print(f"\nFinal result: {'COMPLETE ✓' if success else 'INCOMPLETE (has sorry)'}")
    print(f"File: {abs_path}")

    if gateway_reachable is not None:
        print(f"AI Gateway reachable: {'YES' if gateway_reachable else 'NO'}")

    if proven_tactics:
        print(f"\nProven tactics ({len(proven_tactics)}):")
        for tactic in proven_tactics:
            print(f"  • {tactic}")

    remaining = count_sorrys(final_content)
    print(f"\nRemaining sorries: {remaining}")

    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
