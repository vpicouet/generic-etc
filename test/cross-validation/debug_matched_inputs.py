#!/usr/bin/env python3
"""
Debug: Compare with MATCHED inputs - use KCWI's sky value, dark=0
Test only parameters that KCWI can actually vary
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
from etc_kcwi import calculate_snr as kcwi_calculate_snr, get_params, sky_cts

etc_kcwi.datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/kcwi/"

# Test configuration
slicer = "Small"
grating = "BL"
gratwave = 4550
arcsec_per_slice = 0.35
A_per_pixel = 0.5

# Get KCWI's sky value at this wavelength
w_sky = np.arange(3000.0, 6000.0, 1.0)
c_s_array = sky_cts(w_sky, grating, 1.0)  # Per second
idx_sky = np.argmin(np.abs(w_sky - gratwave))
sky_per_s_per_A_per_arcsec2 = c_s_array[idx_sky]

# Convert to erg/cm²/s/Å/arcsec²
# KCWI gives electrons, we need to work backwards to get erg
# From KCWI: c_s = sky_cts(w, grating, exptime) which already includes efficiency
# But we need the sky in erg units for Generic ETC
# Let's extract it from KCWI's calculation

# Actually, let's just get the sky in electrons and convert it to Generic ETC's expected format
# Generic ETC expects Sky in erg/cm²/s/Å/arcsec²
# KCWI gives us electrons/s/Å/arcsec² (with efficiency already applied)

# We need to reverse engineer: sky_electrons = efficiency × A_geo × sky_photons
# sky_photons = sky_erg / E_photon
# So: sky_erg = sky_electrons / (efficiency × A_geo) × E_photon

H = 6.62607015e-27
C = 2.99792458e10
A_GEO_KECK = np.pi / 4.0 * (1000.0) ** 2

E_photon = H * C / (gratwave * 1e-8)

eff_keck = get_params(w_sky, grating)
eff_keck_at_wave = eff_keck[idx_sky]

# sky_photons/s/Å/arcsec² = sky_electrons/s/Å/arcsec² / (eff × A_geo)
sky_photons_per_s = sky_per_s_per_A_per_arcsec2 / (eff_keck_at_wave * A_GEO_KECK)
sky_erg = sky_photons_per_s * E_photon

print("="*80)
print("MATCHED INPUT TEST: Using KCWI's sky value, dark=0")
print("="*80)
print(f"Configuration: {slicer} slicer, {grating} grating, {gratwave}A")
print(f"\nKCWI sky value:")
print(f"  Sky (electrons/s/Å/arcsec²) = {sky_per_s_per_A_per_arcsec2:.3f}")
print(f"  Sky (erg/cm²/s/Å/arcsec²) = {sky_erg:.3e}")
print(f"  Dark current = 0.0 (KCWI doesn't have dark)")
print("="*80)

# Nominal parameters
seeing_nominal = 0.75
t_exp_nominal = 3600.0
flux_nominal = 5.6e-19
rn_nominal = 2.7
dark_current = 0.0  # Fixed to 0

results = []

def compare_snr_matched(param_name, param_value, slicer, grating, gratwave, seeing, t_exp,
                        flux_erg, read_noise):
    """Compare SNR with matched sky and dark=0."""

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

    # Generic ETC - USE KCWI's SKY VALUE!
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
            Sky=sky_erg,  # USE KCWI's SKY!
            RN=read_noise,
            Dark_current=0.0,  # NO DARK
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
            'gen_signal_el': signal_el,
            'gen_sky_el': sky_el,
            'gen_rn_final': rn_final,
            'gen_pixels_total': pixels_total,
            'gen_pixels_used': pixels_used,
            'gen_signal_total': signal_el * pixels_total,
            'gen_sky_total': sky_el * pixels_total,
            'gen_rn2_total': rn_final**2 * pixels_used,
            'seeing': seeing,
            't_exp': t_exp,
            'flux_erg': flux_erg,
            'sky_erg': sky_erg,
            'read_noise': read_noise,
            'efficiency': eff_keck_at_wave,
        }

        return result

    except Exception as e:
        print(f"Generic ETC error: {e}")
        import traceback
        traceback.print_exc()
        return None


# Test only parameters that KCWI can vary

# Test 1: Exposure time
print("\n1. Testing exposure time variation...")
t_exp_values = [1000, 3600, 10000]
for t_exp in t_exp_values:
    res = compare_snr_matched('t_exp', t_exp, slicer, grating, gratwave, seeing_nominal,
                              t_exp, flux_nominal, rn_nominal)
    if res:
        results.append(res)
        print(f"  t_exp={t_exp:5.0f}s: SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Test 2: Source flux
print("\n2. Testing source flux variation...")
flux_factors = [0.5, 1.0, 2.0]
for factor in flux_factors:
    flux = flux_nominal * factor
    res = compare_snr_matched('flux', flux*1e19, slicer, grating, gratwave, seeing_nominal,
                              t_exp_nominal, flux, rn_nominal)
    if res:
        results.append(res)
        print(f"  flux={flux*1e19:.2f}e-19: SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Test 3: Seeing
print("\n3. Testing seeing variation...")
seeing_values = [0.5, 0.75, 1.0]
for seeing in seeing_values:
    res = compare_snr_matched('seeing', seeing, slicer, grating, gratwave, seeing,
                              t_exp_nominal, flux_nominal, rn_nominal)
    if res:
        results.append(res)
        print(f"  seeing={seeing:.2f}\": SNR_gen={res['SNR_generic']:.3f}, SNR_kcwi={res['SNR_kcwi']:.3f}, diff={res['SNR_diff_pct']:+.1f}%")

# Save to CSV
df = pd.DataFrame(results)
csv_path = '/Users/Vincent/Github/generic-etc/test/cross-validation/debug_matched_inputs.csv'
df.to_csv(csv_path, index=False)
print(f"\n✓ Saved detailed results to: {csv_path}")

# Summary
print("\n" + "="*80)
print("SUMMARY (with matched sky & dark=0)")
print("="*80)
print(f"Total tests: {len(results)}")
print(f"SNR difference range: {df['SNR_diff_pct'].min():.1f}% to {df['SNR_diff_pct'].max():.1f}%")
print(f"Mean SNR difference: {df['SNR_diff_pct'].mean():.1f}%")
print(f"Std SNR difference: {df['SNR_diff_pct'].std():.1f}%")
print("\nIf differences persist, they come from:")
print("  - Different pixel counting methods")
print("  - Different spatial integration")
print("  - Numerical precision differences")
