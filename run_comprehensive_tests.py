#!/usr/bin/env python3
"""
Comprehensive test to identify cross-validation issues
Tests Generic ETC vs Astropy and identifies the exact problem
"""

import sys
sys.path.insert(0, 'notebooks')

import numpy as np
import pandas as pd

def run_comprehensive_tests():
    """Run all tests and identify problems"""

    try:
        from Observation import Observation, load_instruments
        from astropy.stats import signal_to_noise_oir_ccd
    except ImportError as e:
        print(f"Missing dependencies: {e}")
        print("This script needs to run in an environment with numpy, astropy, pandas")
        return

    # Load instruments
    instruments, database = load_instruments()
    instrument_name = "GALEX FUV"

    # Extract parameters
    params = {}
    for i, charact in enumerate(instruments['Charact.']):
        if charact and not isinstance(charact, np.ma.core.MaskedConstant):
            value = instruments[instrument_name][i]
            if not isinstance(value, np.ma.core.MaskedConstant):
                params[charact] = value

    print("="*80)
    print("COMPREHENSIVE CROSS-VALIDATION TEST")
    print("="*80)

    # Helper function for Astropy calculation
    def calc_astropy_snr(signal_flux, sky_flux, dark_current, t_exp, params):
        """Calculate Astropy SNR for given parameters"""
        wavelength_nm = params['wavelength']
        E_photon = 1.986e-8 / wavelength_nm

        signal_photons = signal_flux / E_photon
        sky_photons = sky_flux / E_photon

        collecting_area_cm2 = params['Collecting_area'] * 1e4
        pixel_area_arcsec2 = params['pixel_scale']**2
        wavelength_range = params['dispersion']
        total_throughput = params['Throughput'] * params['QE'] * params['Atmosphere']

        source_eps = (signal_photons * collecting_area_cm2 *
                     pixel_area_arcsec2 * wavelength_range * total_throughput)
        sky_eps = (sky_photons * collecting_area_cm2 *
                  pixel_area_arcsec2 * wavelength_range * total_throughput)
        dark_eps = dark_current / 3600.0

        snr = signal_to_noise_oir_ccd(
            t=t_exp, source_eps=source_eps, sky_eps=sky_eps,
            dark_eps=dark_eps, rd=params['RN'], npix=1, gain=1.0
        )

        return {
            'snr': snr,
            'source_eps': source_eps,
            'sky_eps': sky_eps,
            'dark_eps': dark_eps,
            'signal_total': source_eps * t_exp,
            'sky_total': sky_eps * t_exp,
            'dark_total': dark_eps * t_exp,
        }

    # TEST 1: Baseline comparison
    print("\n" + "="*80)
    print("TEST 1: BASELINE - Default parameters")
    print("="*80)

    t_exp = 1000.0
    obs = Observation(
        instruments=instruments, instrument=instrument_name,
        exposure_time=t_exp, SNR_res="per pix", IFS=False, test=True
    )

    astropy = calc_astropy_snr(
        params['Signal'], params['Sky'], params['Dark_current'], t_exp, params
    )

    print(f"\nGeneric ETC (brut):")
    print(f"  Signal_el: {obs.Signal_el[obs.i]:.6e} e⁻")
    print(f"  sky: {obs.sky[obs.i]:.6e} e⁻")
    print(f"  Dark_current_f: {obs.Dark_current_f[obs.i]:.6e} e⁻")
    print(f"  SNR: {obs.SNR[obs.i]:.6e}")
    print(f"  pixels_total_source: {obs.pixels_total_source:.2f}")

    print(f"\nAstropy (per pixel):")
    print(f"  signal_total: {astropy['signal_total']:.6e} e⁻")
    print(f"  sky_total: {astropy['sky_total']:.6e} e⁻")
    print(f"  dark_total: {astropy['dark_total']:.6e} e⁻")
    print(f"  SNR: {astropy['snr']:.6e}")

    # Apply correction
    signal_corrected = obs.Signal_el[obs.i] * obs.pixels_total_source
    sky_corrected = obs.sky[obs.i] * obs.pixels_total_source

    print(f"\nGeneric ETC (corrected × {obs.pixels_total_source:.2f}):")
    print(f"  Signal_el: {signal_corrected:.6e} e⁻")
    print(f"  sky: {sky_corrected:.6e} e⁻")

    print(f"\nRatios (corrected/astropy):")
    ratio_signal = signal_corrected / astropy['signal_total']
    ratio_sky = sky_corrected / astropy['sky_total']
    ratio_dark = obs.Dark_current_f[obs.i] / astropy['dark_total']
    print(f"  Signal: {ratio_signal:.4f}")
    print(f"  Sky: {ratio_sky:.4f}")
    print(f"  Dark: {ratio_dark:.4f}")

    baseline_issue = {
        'signal_ok': abs(ratio_signal - 1) < 0.15,
        'sky_ok': abs(ratio_sky - 1) < 0.15,
        'dark_ok': abs(ratio_dark - 1) < 0.15,
    }

    # TEST 2: Varying exposure time
    print("\n" + "="*80)
    print("TEST 2: EXPOSURE TIME VARIATION")
    print("="*80)

    exp_times = np.logspace(1, 4, 10)
    results_exp = []

    for t in exp_times:
        obs = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=t, SNR_res="per pix", IFS=False, test=True
        )
        astropy = calc_astropy_snr(
            params['Signal'], params['Sky'], params['Dark_current'], t, params
        )

        results_exp.append({
            'exp_time': t,
            'snr_generic': obs.SNR[obs.i],
            'snr_astropy': astropy['snr'],
            'ratio': obs.SNR[obs.i] / astropy['snr'],
            'pixels_total_source': obs.pixels_total_source,
        })

    df_exp = pd.DataFrame(results_exp)
    print(f"\nSNR ratio variation: {df_exp['ratio'].std():.6f}")
    print(f"pixels_total_source variation: {df_exp['pixels_total_source'].std():.6f}")
    print(f"SNR ratio range: [{df_exp['ratio'].min():.4f}, {df_exp['ratio'].max():.4f}]")

    # TEST 3: Varying signal flux
    print("\n" + "="*80)
    print("TEST 3: SIGNAL FLUX VARIATION")
    print("="*80)

    signal_factors = np.logspace(-1, 1, 10)
    results_signal = []

    for factor in signal_factors:
        signal_flux = params['Signal'] * factor

        # Modify instrument
        idx = list(instruments['Charact.']).index('Signal')
        original = instruments[instrument_name][idx]
        instruments[instrument_name][idx] = signal_flux

        obs = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=1000.0, SNR_res="per pix", IFS=False, test=True
        )
        astropy = calc_astropy_snr(signal_flux, params['Sky'], params['Dark_current'], 1000.0, params)

        instruments[instrument_name][idx] = original

        results_signal.append({
            'signal_factor': factor,
            'snr_generic': obs.SNR[obs.i],
            'snr_astropy': astropy['snr'],
            'ratio': obs.SNR[obs.i] / astropy['snr'],
            'pixels_total_source': obs.pixels_total_source,
        })

    df_signal = pd.DataFrame(results_signal)
    print(f"\nSNR ratio variation: {df_signal['ratio'].std():.6f}")
    print(f"pixels_total_source variation: {df_signal['pixels_total_source'].std():.6f}")
    print(f"SNR ratio range: [{df_signal['ratio'].min():.4f}, {df_signal['ratio'].max():.4f}]")

    # TEST 4: Varying sky background
    print("\n" + "="*80)
    print("TEST 4: SKY BACKGROUND VARIATION")
    print("="*80)

    sky_factors = np.logspace(-1, 1, 10)
    results_sky = []

    for factor in sky_factors:
        sky_flux = params['Sky'] * factor

        idx = list(instruments['Charact.']).index('Sky')
        original = instruments[instrument_name][idx]
        instruments[instrument_name][idx] = sky_flux

        obs = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=1000.0, SNR_res="per pix", IFS=False, test=True
        )
        astropy = calc_astropy_snr(params['Signal'], sky_flux, params['Dark_current'], 1000.0, params)

        instruments[instrument_name][idx] = original

        sky_corrected = obs.sky[obs.i] * obs.pixels_total_source

        results_sky.append({
            'sky_factor': factor,
            'snr_generic': obs.SNR[obs.i],
            'snr_astropy': astropy['snr'],
            'ratio': obs.SNR[obs.i] / astropy['snr'],
            'pixels_total_source': obs.pixels_total_source,
            'sky_generic': obs.sky[obs.i],
            'sky_corrected': sky_corrected,
            'sky_astropy': astropy['sky_total'],
            'sky_ratio': sky_corrected / astropy['sky_total'],
        })

    df_sky = pd.DataFrame(results_sky)
    print(f"\nSNR ratio variation: {df_sky['ratio'].std():.6f}")
    print(f"pixels_total_source variation: {df_sky['pixels_total_source'].std():.6f}")
    print(f"SNR ratio range: [{df_sky['ratio'].min():.4f}, {df_sky['ratio'].max():.4f}]")
    print(f"Sky ratio (corrected) range: [{df_sky['sky_ratio'].min():.4f}, {df_sky['sky_ratio'].max():.4f}]")

    # TEST 5: Varying dark current
    print("\n" + "="*80)
    print("TEST 5: DARK CURRENT VARIATION")
    print("="*80)

    dark_factors = np.logspace(-1, 1, 10)
    results_dark = []

    for factor in dark_factors:
        dark_current = params['Dark_current'] * factor

        idx = list(instruments['Charact.']).index('Dark_current')
        original = instruments[instrument_name][idx]
        instruments[instrument_name][idx] = dark_current

        obs = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=1000.0, SNR_res="per pix", IFS=False, test=True
        )
        astropy = calc_astropy_snr(params['Signal'], params['Sky'], dark_current, 1000.0, params)

        instruments[instrument_name][idx] = original

        results_dark.append({
            'dark_factor': factor,
            'snr_generic': obs.SNR[obs.i],
            'snr_astropy': astropy['snr'],
            'ratio': obs.SNR[obs.i] / astropy['snr'],
            'pixels_total_source': obs.pixels_total_source,
            'dark_generic': obs.Dark_current_f[obs.i],
            'dark_astropy': astropy['dark_total'],
            'dark_ratio': obs.Dark_current_f[obs.i] / astropy['dark_total'],
        })

    df_dark = pd.DataFrame(results_dark)
    print(f"\nSNR ratio variation: {df_dark['ratio'].std():.6f}")
    print(f"pixels_total_source variation: {df_dark['pixels_total_source'].std():.6f}")
    print(f"SNR ratio range: [{df_dark['ratio'].min():.4f}, {df_dark['ratio'].max():.4f}]")
    print(f"Dark ratio range: [{df_dark['dark_ratio'].min():.4f}, {df_dark['dark_ratio'].max():.4f}]")

    # FINAL DIAGNOSIS
    print("\n" + "="*80)
    print("FINAL DIAGNOSIS")
    print("="*80)

    print(f"\n1. BASELINE TEST:")
    if baseline_issue['signal_ok']:
        print(f"   ✓ Signal OK after correction (ratio={ratio_signal:.4f})")
    else:
        print(f"   ✗ Signal still wrong after correction (ratio={ratio_signal:.4f})")

    if baseline_issue['sky_ok']:
        print(f"   ✓ Sky OK after correction (ratio={ratio_sky:.4f})")
    else:
        print(f"   ✗ Sky still wrong after correction (ratio={ratio_sky:.4f})")

    if baseline_issue['dark_ok']:
        print(f"   ✓ Dark OK (ratio={ratio_dark:.4f})")
    else:
        print(f"   ✗ Dark wrong (ratio={ratio_dark:.4f})")

    print(f"\n2. VARIATION TESTS:")

    if df_exp['ratio'].std() < 0.01:
        print(f"   ✓ Exposure time: SNR ratio constant (std={df_exp['ratio'].std():.6f})")
    else:
        print(f"   ✗ Exposure time: SNR ratio varies (std={df_exp['ratio'].std():.6f})")

    if df_signal['ratio'].std() < 0.01:
        print(f"   ✓ Signal flux: SNR ratio constant (std={df_signal['ratio'].std():.6f})")
    else:
        print(f"   ✗ Signal flux: SNR ratio varies (std={df_signal['ratio'].std():.6f})")

    if df_sky['ratio'].std() < 0.01:
        print(f"   ✓ Sky background: SNR ratio constant (std={df_sky['ratio'].std():.6f})")
    else:
        print(f"   ✗ Sky background: SNR ratio varies (std={df_sky['ratio'].std():.6f})")
        print(f"      Sky ratio (corrected) varies: [{df_sky['sky_ratio'].min():.4f}, {df_sky['sky_ratio'].max():.4f}]")

    if df_dark['ratio'].std() < 0.01:
        print(f"   ✓ Dark current: SNR ratio constant (std={df_dark['ratio'].std():.6f})")
    else:
        print(f"   ✗ Dark current: SNR ratio varies (std={df_dark['ratio'].std():.6f})")

    print(f"\n3. pixels_total_source STABILITY:")
    for name, df in [('exposure', df_exp), ('signal', df_signal), ('sky', df_sky), ('dark', df_dark)]:
        pct_var = df['pixels_total_source'].std() / df['pixels_total_source'].mean() * 100
        if pct_var < 0.1:
            print(f"   ✓ {name}: constant (var={pct_var:.4f}%)")
        else:
            print(f"   ✗ {name}: varies (var={pct_var:.4f}%)")

    # Export results
    df_exp.to_csv('test_exposure_variation.csv', index=False)
    df_signal.to_csv('test_signal_variation.csv', index=False)
    df_sky.to_csv('test_sky_variation.csv', index=False)
    df_dark.to_csv('test_dark_variation.csv', index=False)
    print(f"\n✓ Results exported to CSV files")

    return {
        'baseline': baseline_issue,
        'df_exp': df_exp,
        'df_signal': df_signal,
        'df_sky': df_sky,
        'df_dark': df_dark,
    }

if __name__ == "__main__":
    results = run_comprehensive_tests()
