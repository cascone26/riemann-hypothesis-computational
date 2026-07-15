#!/usr/bin/env bash
# OLYMPUS_RELAUNCH.sh — Run this on Olympus the moment it comes back online.
# Kills any stale jobs, re-launches all pending experiments.
# Run from Sam: ssh olympus "bash ~/riemann/operators/OLYMPUS_RELAUNCH.sh"

set -e
cd ~/riemann

echo "=== OLYMPUS RELAUNCH $(date) ==="
echo "Killing stale jobs..."
pkill -f "gamma2_resonance.py" 2>/dev/null || true
pkill -f "resonance_theory.py" 2>/dev/null || true
pkill -f "sensitivity_extended.py" 2>/dev/null || true
pkill -f "resonance_lattice_scan.py" 2>/dev/null || true
pkill -f "lsq250_dense.py" 2>/dev/null || true
pkill -f "extended_powerlaw.py" 2>/dev/null || true
pkill -f "gamma3_scan.py" 2>/dev/null || true
sleep 2

echo "Launching all experiments..."

# 1. γ₂ main resonance scan (lsq=380-500, N=min(5*lsq,2000))
nohup python3 operators/gamma2_resonance.py \
    > results/gamma2_resonance.out 2>&1 &
echo "  [1] gamma2_resonance PID=$!"

# 2. γ₃ resonance scan (lsq=400-800, N=min(5*lsq,3200)) — NEW
nohup python3 operators/gamma3_scan.py \
    > results/gamma3_scan.out 2>&1 &
echo "  [2] gamma3_scan PID=$!"

# 3. Extended power-law fit (lsq=100-350, N=min(5*lsq,2000)) — NEW
nohup python3 operators/extended_powerlaw.py \
    > results/extended_powerlaw.out 2>&1 &
echo "  [3] extended_powerlaw PID=$!"

# 4. Resonance lattice scan (lsq=240-350 with labels)
nohup python3 operators/resonance_lattice_scan.py \
    > results/resonance_lattice_scan.out 2>&1 &
echo "  [4] resonance_lattice_scan PID=$!"

# 5. Dense lsq=230-270 scan (N=min(5*lsq,1500))
nohup python3 operators/lsq250_dense.py \
    > results/lsq250_dense.out 2>&1 &
echo "  [5] lsq250_dense PID=$!"

# 6. Sensitivity extended: k=16-30 at lsq=50
nohup python3 operators/sensitivity_extended.py \
    > results/sensitivity_extended.out 2>&1 &
echo "  [6] sensitivity_extended PID=$!"

# 7. Per-zero resonance theory
nohup python3 operators/resonance_theory.py \
    > results/resonance_theory.out 2>&1 &
echo "  [7] resonance_theory PID=$!"

# 8. Check c=100 chain status
echo ""
echo "=== Chain status ==="
ls -lh results/chain_log.txt 2>/dev/null && tail -5 results/chain_log.txt || echo "  No chain_log.txt"
ls -lh results/lam_min_c100.json 2>/dev/null && python3 -c "import json; d=json.load(open('results/lam_min_c100.json')); print(f'  c=100: lam_min={d[\"lam_min\"]:.4e}')" || echo "  c=100 not done yet"

echo ""
echo "=== Jobs running ==="
jobs
echo ""
echo "=== DONE launching. Monitor with: tail -f results/*.out ==="
