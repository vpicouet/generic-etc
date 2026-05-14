#!/usr/bin/env python3
"""
Quick test to verify dark current variation works correctly.
"""

import sys
import os
import numpy as np
from pathlib import Path

# Setup paths
script_dir = "/Users/Vincent/Github/generic-etc/test/cross-validation"
notebooks_dir = str(Path(script_dir).parent.parent / "notebooks")
os.chdir(notebooks_dir)

sys.path.insert(0, ".")
sys.path.insert(0, "/Users/Vincent/Github/exposureTimeCalculator")

from plot_kcwi_snr_ratios_fast import calculate_snr_ratio

# Test dark variation
print("="*60)
print("Testing Dark Current Variation")
print("="*60)

dark_values = [0.0, 5.0, 10.0]
for dark in dark_values:
    pct_diff, snr_gen, snr_keck, sig, sky = calculate_snr_ratio(
        "Small", "BL", 4550, 0.75, 3600.0,
        5.6e-19, 8e-18, 2.7, dark, verbose=False)

    print(f"\nDark = {dark:.1f} e⁻/pix/hr:")
    print(f"  SNR Generic: {snr_gen:.4f}")
    print(f"  SNR KCWI:    {snr_keck:.4f}")
    print(f"  Ratio:       {snr_gen/snr_keck:.4f}")
    print(f"  Diff:        {pct_diff:.2f}%")

print("\n" + "="*60)
print("Expected: SNR Generic should DECREASE as dark increases")
print("KCWI SNR should stay constant (KCWI has no dark)")
print("="*60)
