#!/usr/bin/env python3
"""
Debug script: Compare Generic ETC vs KCWI ETC SNR calculations
Simplified: 1 config, 3 points per parameter, save CSV for analysis
"""

import sys
import os
import numpy as np
import pandas as pd
from pathlib import Path

# Setup paths
script_dir = "/Users/Vincent/Github/generic-etc/test/cross-validation"
notebooks_dir = str(Path(script_dir).parent.parent / "notebooks")
os.chdir(notebooks_dir)

sys.path.insert(0, ".")
sys.path.insert(0, "/Users/Vincent/Github/exposureTimeCalculator")

from Observation import load_instruments, Observation
import etc_kcwi
from etc_kcwi import calculate_snr as kcwi_calculate_snr, get_params

etc_kcwi.datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/kcwi/"

# Test configuration: Small slicer + BL grating
slicer = "Small"
grating = "BL"
gratwave = 4550
arcsec_per_slice = 0.35
A_per_pixel = 0.5

# Nominal parameters
seeing_nominal = 0.75
t_exp_nominal = 3600.0
flux_nominal = 5.6e-19
sky_nominal = 8e-18
rn_nominal = 2.7
dark_nominal = 0.0  # SET TO ZERO for all tests except dark variation

print("="*80)
print("DEBUG: Generic ETC vs KCWI ETC SNR Comparison")
print("="*80)
print(f"Configuration: {slicer} slicer, {grating} grating, {gratwave}A")
print(f"Dark current: {dark_nominal} e-/pix/hr (set to 0 for debugging)")
print("="*80)

results = []

def compare_snr(param_name, param_value, slicer, grating, gratwave, seeing, t_exp,
                flux_erg, sky_erg, read_noise, dark_current):
    """Compare SNR from both ETCs and return detailed breakdown."""

    # KCWI ETC
    try:
        w_kcwi, snr_kcwi = kcwi_calculate_snr(
            slicer=slicer,
            grating=grating,
            gratwave=gratwave,
            seeing=seeing,
            exptime=t_exp,
            flux=flux_erg,
            read_noise=read_noise,
            ccdbin='1x1'
        )
        idx = np.argmin(np.abs(w_kcwi - gratwave))
        snr_keck = snr_kcwi[idx]

        # Get efficiency
        eff_keck = get_params(w_kcwi, grating)
        eff_keck_at_wave = eff_keck[idx]

    except Exception as e:
        print(f"KCWI error: {e}")
        return None

    # Generic ETC
    try:
        instruments, _ = load_instruments()

        obs = Observation(
            instruments=instruments,
            instrument="KCWI blue",
            exposure_time=t_exp,
            acquisition_time=t_exp/3600.0,
            SNR_res="per Source per 2λpix",
            IFS=True,
            test=True,
            wavelength=gratwave / 10,
            Slitwidth=arcsec_per_slice,
            dispersion=A_per_pixel / 10,
            Size_source=seeing,
            Signal=flux_erg,
            Sky=sky_erg,
            RN=read_noise,
            Dark_current=dark_current,
            Throughput=eff_keck_at_wave,
            QE=1.0,
            Atmosphere=1.0,
        )

        snr_generic = obs.SNR[obs.i]

        # Get detailed breakdown
        signal_el = obs.Signal_el[obs.i]
        sky_el = obs.sky[obs.i]
        rn_final = obs.RN_final[obs.i]
        pixels_total = obs.pixels_total_source
        pixels_used = obs.number_pixels_used

        result = {
            'param_name': param_name,
            'param_value': param_value,
            'SNR_generic': snr_generic,
            'SNR_kcwi': snr_keck,
            'SNR_ratio': snr_generic / snr_keck,
            'SNR_diff_pct': 100.0 * (snr_generic / snr_keck - 1.0),
            # Generic ETC breakdown
            'gen_signal_el': signal_el,
            'gen_sky_el': sky_el,
            'gen_rn_final': rn_final,
            'gen_pixels_total': pixels_total,
            'gen_pixels_used': pixels_used,
            'gen_signal_total': signal_el * pixels_total,
            'gen_sky_total': sky_el * pixels_total,
            'gen_rn2_total': rn_final**2 * pixels_used,
            # Inputs
            'seeing': seeing,
            't_exp': t_exp,
            'flux_erg': flux_erg,
            'sky_erg': sky_erg,
            'read_noise': read_noise,
            'dark_current': dark_current,
            'efficiency': eff_keck_at_wave,
        }

        return result

    except Exception as e:
        print(f"Generic ETC error: {e}")
        import traceback
        traceback.print_exc()
        return None


# Test 1: Exposure time (3 points)
print("\n1. Testing exposure time variation...")
t_exp_values = [1000, 3600, 10000]
for t_exp in t_exp_values:
    res = compare_snr('t_exp', t_exp, slicer, grating, gratwave, seeing_nominal,
                     t_exp, flux_nominal, sky_nominal, rn_nominal, dark_nominal)
    if res:
        results.append(res)
        print(f"  t_exp={t_exp:5.0f}s: SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Test 2: Source flux (3 points)
print("\n2. Testing source flux variation...")
flux_factors = [0.5, 1.0, 2.0]
for factor in flux_factors:
    flux = flux_nominal * factor
    res = compare_snr('flux', flux*1e19, slicer, grating, gratwave, seeing_nominal,
                     t_exp_nominal, flux, sky_nominal, rn_nominal, dark_nominal)
    if res:
        results.append(res)
        print(f"  flux={flux*1e19:.2f}e-19: SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Test 3: Sky background (3 points)
print("\n3. Testing sky background variation...")
sky_factors = [0.5, 1.0, 2.0]
for factor in sky_factors:
    sky = sky_nominal * factor
    res = compare_snr('sky', sky*1e18, slicer, grating, gratwave, seeing_nominal,
                     t_exp_nominal, flux_nominal, sky, rn_nominal, dark_nominal)
    if res:
        results.append(res)
        print(f"  sky={sky*1e18:.2f}e-18: SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Test 4: Read noise (3 points)
print("\n4. Testing read noise variation...")
rn_values = [1.0, 2.7, 5.0]
for rn in rn_values:
    res = compare_snr('read_noise', rn, slicer, grating, gratwave, seeing_nominal,
                     t_exp_nominal, flux_nominal, sky_nominal, rn, dark_nominal)
    if res:
        results.append(res)
        print(f"  RN={rn:.1f}e-: SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Test 5: Dark current (3 points) - NOW we add dark to Generic ETC!
print("\n5. Testing dark current variation (WITH dark in Generic ETC)...")
dark_values = [0.0, 5.0, 10.0]
for dark in dark_values:
    res = compare_snr('dark_current', dark, slicer, grating, gratwave, seeing_nominal,
                     t_exp_nominal, flux_nominal, sky_nominal, rn_nominal, dark)
    if res:
        results.append(res)
        print(f"  Dark={dark:.1f}e-/pix/hr: SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Save to CSV
df = pd.DataFrame(results)
csv_path = '/Users/Vincent/Github/generic-etc/test/cross-validation/debug_snr_comparison.csv'
df.to_csv(csv_path, index=False)
print(f"\n✓ Saved detailed results to: {csv_path}")

# Summary
print("\n" + "="*80)
print("SUMMARY")
print("="*80)
print(f"Total tests: {len(results)}")
print(f"SNR difference range: {df['SNR_diff_pct'].min():.1f}% to {df['SNR_diff_pct'].max():.1f}%")
print(f"Mean SNR difference: {df['SNR_diff_pct'].mean():.1f}%")
print("\nCheck the CSV file for detailed breakdown of all parameters!")
