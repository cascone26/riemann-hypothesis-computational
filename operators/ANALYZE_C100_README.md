# analyze_c100.py — CvS High-Precision Results Analysis

## Purpose

Post-processing script for Riemann zeta CvS (Connes-Vareille-Séries) computation results from Olympus (after c=100, T=800, dps=150 job completes).

**Workflow**:
1. SSH fetch `/home/cobo/riemann/results/cvs_highprec_c100.json` from Olympus
2. Extract `lam_min` (ground state eigenvalue) and extracted Riemann zeros
3. Compare against known reference zeros (γ₁..γ₁₀)
4. Compute digit precision for each zero
5. Detect outliers (>2σ below mean precision)
6. Special analysis of k=7 zero: compute |sin(γ₇ × log(100) / 2)|
7. Verify chain log shows T=1200 and extraction jobs launched
8. Append findings to `SESSION_FINDINGS_2026_07_11.md` as section §49

## Usage

```bash
cd /Users/scones/projects/riemann/
python3 operators/analyze_c100.py
```

### Prerequisites
- Python 3.8+
- `mpmath` (for high-precision arithmetic)
- SSH access to `olympus` (configured in `~/.ssh/config` or Tailscale alias)
- Results file exists at `/home/cobo/riemann/results/cvs_highprec_c100.json`

### Install mpmath
```bash
pip install mpmath
```

## Output

**Console**: Step-by-step progress report with metrics:
```
[2026-07-12T...] Starting c=100 analysis...
[1/5] Fetching results from olympus:/home/cobo/riemann/results/cvs_highprec_c100.json
[2/5] Inspecting JSON structure...
JSON structure:
  lam_min: float = -1.234e-10
  zeros: list (len=10)
  ...
[3/5] Computing precision metrics...
  ✓ Analyzed 10 zeros
[4/5] Detecting outliers...
  ✓ Found 2 outliers
[5/5] Fetching chain log for verification...
[6/6] Appending to /Users/scones/projects/riemann/results/SESSION_FINDINGS_2026_07_11.md
  ✓ Findings appended

====================================================================
ANALYSIS COMPLETE
====================================================================
lam_min (sign): -1.234e-10 (-)
Zeros extracted: 10
Mean digit precision: 12.34
Outliers: 2
Findings written to: /Users/scones/projects/riemann/results/SESSION_FINDINGS_2026_07_11.md
====================================================================
```

**File Output**: Appends section §49 to `SESSION_FINDINGS_2026_07_11.md`:

```markdown
## §49. CvS High-Precision Analysis (c=100, T=800, dps=150)

**Date**: 2026-07-12T...

### lam_min
- **Value**: -1.234e-10
- **Sign**: -

### Extracted Riemann Zeros

**Count**: 10 zeros

| k | γ_extracted | γ_ref | error | digit_precision | sign |
|---|---|---|---|---|---|
| 1 | 14.134725142 | 14.134725142 | 1.2e-12 | 11.92 | + |
...

### Precision Summary
- **Mean digit precision**: 12.34 digits
- **Min digit precision**: 9.45 digits (zero k=7)
- **Max digit precision**: 14.12 digits (zero k=3)

### Outliers (>2σ below mean precision)

**Count**: 2 outliers detected.

| k | z_score | digit_precision | error |
|---|---|---|---|
| 7 | 2.45 | 9.45 | 3.6e-10 |
...

### k=7 Special Analysis
γ₇ = 40.918719...
L = log(100) = 4.6051...
|sin(γ₇L/2)| = 1.234e-15
lam_min = -1.234e-10

### Chain Log Verification
Chain log (first 500 chars):
```
[chain log excerpt]
```
```

## JSON Input Format

The script expects `/home/cobo/riemann/results/cvs_highprec_c100.json` to contain (at minimum):

```json
{
  "lam_min": -1.234e-10,
  "zeros": [14.134725..., 21.022039..., ...],
  ...other fields...
}
```

**Supported key variants**:
- `lam_min`, `lambda_min`, or `ground_state_eigenvalue`
- `zeros`, `gamma_k`, `extracted_zeros`, or `riemann_zeros`

The script auto-detects structure; if keys are unexpected, it prints the JSON top-level keys for debugging.

## Known Riemann Zeros (Reference)

Hardcoded 10 zeros used for comparison:

| k | γ_k |
|---|---|
| 1 | 14.134725141734693790457251983562 |
| 2 | 21.022039638771554992628479593896 |
| 3 | 25.010857580145688763213790992562 |
| 4 | 30.424876125859513210311897530584 |
| 5 | 32.935061587739189690662368964074 |
| 6 | 37.586178158825671257217763480705 |
| 7 | 40.918719012147495187398126914633 |
| 8 | 43.327073280914999519496122165406 |
| 9 | 48.005150881167159727942472749427 |
| 10 | 49.773832477672302181916784678563 |

## Error Handling

- **SSH failure**: Reports which host/path failed; suggests SSH config check
- **Invalid JSON**: Reports decode error and raw file path
- **Missing keys**: Prints JSON structure to help debug
- **No chain log**: Notes gracefully in output (log may not exist yet)
- **No zero match**: Silently skips k values not in reference table

## High-Precision Arithmetic

Uses `mpmath` with `dps=200` (200 decimal places) for:
- Error computation: `error = |γ_extracted - γ_ref|`
- Precision: `digit_precision = -log₁₀(error)`
- k=7 analysis: `sin(γ₇ × log(100) / 2)`

## Troubleshooting

### "ssh command not found"
Install OpenSSH (macOS: `brew install openssh`)

### "SSH timeout"
Check Tailscale connection to Olympus:
```bash
ping olympus
# or explicitly: ping 100.98.199.94
```

### "Invalid JSON"
Check result file exists and is complete:
```bash
ssh olympus ls -lh /home/cobo/riemann/results/cvs_highprec_c100.json
ssh olympus head -c 100 /home/cobo/riemann/results/cvs_highprec_c100.json
```

### "Could not extract lam_min or zeros"
File structure differs from expected. Run with `-v` (if implemented) or check raw JSON:
```bash
ssh olympus cat /home/cobo/riemann/results/cvs_highprec_c100.json | jq . | head -50
```

## Author Notes

- Script is idempotent: safe to re-run after fixing issues
- Findings append to SESSION_FINDINGS_2026_07_11.md (does not overwrite)
- All timestamps are local machine time
- High-precision arithmetic ensures digit_precision values are accurate to >100 digits
