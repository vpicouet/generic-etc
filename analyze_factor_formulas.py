#!/usr/bin/env python3
"""
Analyze which factor formula matches Astropy
The code has 3 different formulas that overwrite each other!
"""

import sys
sys.path.insert(0, 'notebooks')

def analyze_formulas():
    """
    The code in Observation.py lines 2481-2527 has a confusing structure:

    Line 2481: factor_CU2el_tot = ... / pixels_total_source
    Line 2484: factor_CU2el = factor_CU2el_tot
    Line 2488: factor_CU2el = ... * pixel_scale  (OVERWRITES!)
    Line 2504: factor_CU2el_average = ... / pixels_total_source
    Line 2522-2527: if test: use _tot, else use _average (OVERWRITES AGAIN!)

    So which formula is actually used?
    """

    try:
        import numpy as np
        from Observation import Observation, load_instruments

        instruments, database = load_instruments()
        instrument_name = "GALEX FUV"

        params = {}
        for i, charact in enumerate(instruments['Charact.']):
            if charact and not isinstance(charact, np.ma.core.MaskedConstant):
                value = instruments[instrument_name][i]
                if not isinstance(value, np.ma.core.MaskedConstant):
                    params[charact] = value

        print("="*80)
        print("ANALYZING FACTOR FORMULAS")
        print("="*80)

        # Create observation with test=True
        obs_test_true = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=1000.0, SNR_res="per pix", IFS=False, test=True
        )

        # Create observation with test=False
        obs_test_false = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=1000.0, SNR_res="per pix", IFS=False, test=False
        )

        print(f"\nInstrument parameters:")
        print(f"  Slitwidth: {params['Slitwidth']} arcsec")
        print(f"  Slitlength: {params['Slitlength']} arcsec")
        print(f"  Size_source: {params['Size_source']} arcsec")
        print(f"  pixel_scale: {params['pixel_scale']} arcsec/pix")
        print(f"  dispersion: {params['dispersion']} Å/pix")
        print(f"  Line_width: {params['Line_width']} Å")
        print(f"  Bandwidth: {params['Bandwidth']} Å")

        # Compute what Astropy expects (per pixel)
        collecting_area_cm2 = params['Collecting_area'] * 1e4
        pixel_area_arcsec2 = params['pixel_scale']**2
        wavelength_range = params['dispersion']
        total_throughput = params['Throughput'] * params['QE'] * params['Atmosphere']
        effective_area = collecting_area_cm2 * total_throughput
        arcsec2str = (np.pi/180/3600)**2

        factor_astropy_per_pixel = effective_area * arcsec2str * pixel_area_arcsec2 * wavelength_range

        print(f"\n{'='*80}")
        print(f"ASTROPY EXPECTED (per pixel):")
        print(f"{'='*80}")
        print(f"  factor_per_pixel = effective_area × arcsec2str × pixel_scale² × dispersion")
        print(f"  factor_per_pixel = {effective_area:.2e} × {arcsec2str:.4e} × {pixel_area_arcsec2:.2f} × {wavelength_range:.2f}")
        print(f"  factor_per_pixel = {factor_astropy_per_pixel:.6e}")

        print(f"\n{'='*80}")
        print(f"GENERIC ETC FORMULAS:")
        print(f"{'='*80}")

        print(f"\ntest=True:")
        print(f"  factor_CU2el: {obs_test_true.factor_CU2el:.6e}")
        print(f"  factor_CU2el_sky: {obs_test_true.factor_CU2el_sky:.6e}")
        print(f"  pixels_total_source: {obs_test_true.pixels_total_source:.2f}")
        print(f"  factor_CU2el × pixels_total_source: {obs_test_true.factor_CU2el * obs_test_true.pixels_total_source:.6e}")

        print(f"\ntest=False:")
        print(f"  factor_CU2el: {obs_test_false.factor_CU2el:.6e}")
        print(f"  factor_CU2el_sky: {obs_test_false.factor_CU2el_sky:.6e}")
        print(f"  pixels_total_source: {obs_test_false.pixels_total_source:.2f}")
        print(f"  factor_CU2el × pixels_total_source: {obs_test_false.factor_CU2el * obs_test_false.pixels_total_source:.6e}")

        # Manually compute the 3 different formulas
        print(f"\n{'='*80}")
        print(f"MANUAL CALCULATION OF THE 3 FORMULAS:")
        print(f"{'='*80}")

        # Formula 1: Line 2481 (_tot)
        source_size_arcsec_after_slit = np.minimum(params['Size_source'], params['Slitlength']) * \
                                       np.minimum(params['Size_source'], params['Slitwidth'])
        pixels_total_source = obs_test_true.pixels_total_source

        formula1_signal = (effective_area * arcsec2str *
                          np.minimum(params['Line_width'], params['Bandwidth']) *
                          source_size_arcsec_after_slit / pixels_total_source)

        slit_size_arcsec_after_slit = np.minimum(params['Size_source'], params['Slitlength']) * params['Slitwidth']

        formula1_sky = (effective_area * arcsec2str *
                       np.maximum(np.minimum(params['Line_width'], params['Bandwidth']), params['dispersion']) *
                       slit_size_arcsec_after_slit / pixels_total_source)

        print(f"\nFormula 1 (_tot, line 2481-2483):")
        print(f"  Signal: {formula1_signal:.6e}")
        print(f"  Sky: {formula1_sky:.6e}")
        print(f"  Signal × pixels_total_source: {formula1_signal * pixels_total_source:.6e}")
        print(f"  Sky × pixels_total_source: {formula1_sky * pixels_total_source:.6e}")

        # Formula 2: Line 2488-2492 (middle version that gets overwritten)
        sky_spectral_coverage = np.minimum(params['dispersion'],
                                          np.minimum(params['Line_width'], params['Bandwidth'])/2.35)

        formula2_signal = (effective_area * arcsec2str *
                          np.minimum(params['Slitwidth'], params['Size_source']) *
                          np.minimum(params['dispersion'], np.minimum(params['Line_width'], params['Bandwidth'])/2.35) *
                          params['pixel_scale'])

        formula2_sky = (effective_area * arcsec2str *
                       params['Slitwidth'] *
                       sky_spectral_coverage *
                       params['pixel_scale'])

        print(f"\nFormula 2 (line 2488-2492, gets overwritten):")
        print(f"  Signal: {formula2_signal:.6e}")
        print(f"  Sky: {formula2_sky:.6e}")

        # Formula 3: Line 2504-2514 (_average)
        source_spatial_arcsec = np.minimum(params['Slitwidth'], params['Size_source'])
        source_spectral_angstrom = np.minimum(params['Line_width'], params['Bandwidth'])

        formula3_signal = (effective_area * arcsec2str *
                          source_spatial_arcsec *
                          source_spectral_angstrom /
                          pixels_total_source)

        sky_spatial_arcsec = params['Slitwidth']
        sky_spectral_angstrom = np.minimum(params['Line_width'], params['Bandwidth'])

        formula3_sky = (effective_area * arcsec2str *
                       sky_spatial_arcsec *
                       sky_spectral_angstrom /
                       pixels_total_source)

        print(f"\nFormula 3 (_average, line 2504-2514):")
        print(f"  Signal: {formula3_signal:.6e}")
        print(f"  Sky: {formula3_sky:.6e}")
        print(f"  Signal × pixels_total_source: {formula3_signal * pixels_total_source:.6e}")
        print(f"  Sky × pixels_total_source: {formula3_sky * pixels_total_source:.6e}")

        # Compare with Astropy
        print(f"\n{'='*80}")
        print(f"COMPARISON WITH ASTROPY (ratio = Generic/Astropy):")
        print(f"{'='*80}")

        print(f"\nFormula 1 (_tot) × pixels_total_source:")
        print(f"  Signal ratio: {formula1_signal * pixels_total_source / factor_astropy_per_pixel:.4f}")
        print(f"  Sky ratio: {formula1_sky * pixels_total_source / factor_astropy_per_pixel:.4f}")

        print(f"\nFormula 2 (line 2488-2492):")
        print(f"  Signal ratio: {formula2_signal / factor_astropy_per_pixel:.4f}")
        print(f"  Sky ratio: {formula2_sky / factor_astropy_per_pixel:.4f}")

        print(f"\nFormula 3 (_average) × pixels_total_source:")
        print(f"  Signal ratio: {formula3_signal * pixels_total_source / factor_astropy_per_pixel:.4f}")
        print(f"  Sky ratio: {formula3_sky * pixels_total_source / factor_astropy_per_pixel:.4f}")

        # What is actually used?
        print(f"\n{'='*80}")
        print(f"WHAT IS ACTUALLY USED?")
        print(f"{'='*80}")

        ratio_test_true = obs_test_true.factor_CU2el / formula1_signal
        ratio_test_false = obs_test_false.factor_CU2el / formula3_signal

        if abs(ratio_test_true - 1) < 0.01:
            print(f"\n✓ test=True uses Formula 1 (_tot)")
        elif abs(obs_test_true.factor_CU2el / formula2_signal - 1) < 0.01:
            print(f"\n✗ test=True uses Formula 2 (line 2488-2492) - unexpected!")
        else:
            print(f"\n✗ test=True uses unknown formula!")
            print(f"  Ratio to Formula 1: {ratio_test_true:.4f}")
            print(f"  Ratio to Formula 2: {obs_test_true.factor_CU2el / formula2_signal:.4f}")

        if abs(ratio_test_false - 1) < 0.01:
            print(f"✓ test=False uses Formula 3 (_average)")
        elif abs(obs_test_false.factor_CU2el / formula2_signal - 1) < 0.01:
            print(f"✗ test=False uses Formula 2 (line 2488-2492) - unexpected!")
        else:
            print(f"✗ test=False uses unknown formula!")
            print(f"  Ratio to Formula 3: {ratio_test_false:.4f}")
            print(f"  Ratio to Formula 2: {obs_test_false.factor_CU2el / formula2_signal:.4f}")

        # CONCLUSION
        print(f"\n{'='*80}")
        print(f"CONCLUSION")
        print(f"{'='*80}")

        print(f"\nThe code has overlapping formulas:")
        print(f"1. Lines 2481-2485 set factor_CU2el = factor_CU2el_tot")
        print(f"2. Lines 2488-2492 OVERWRITE factor_CU2el (different formula)")
        print(f"3. Lines 2504-2514 compute factor_CU2el_average")
        print(f"4. Lines 2522-2527 OVERWRITE AGAIN based on test flag")

        print(f"\nFor Astropy compatibility (per pixel):")
        print(f"  Expected: effective_area × arcsec2str × pixel_scale² × dispersion")
        print(f"  Expected value: {factor_astropy_per_pixel:.6e}")

        # Check which formula is closest
        diffs = {
            'Formula 1 × npix': abs(formula1_signal * pixels_total_source - factor_astropy_per_pixel),
            'Formula 2': abs(formula2_signal - factor_astropy_per_pixel),
            'Formula 3 × npix': abs(formula3_signal * pixels_total_source - factor_astropy_per_pixel),
        }

        best = min(diffs, key=diffs.get)
        print(f"\nClosest formula: {best}")
        print(f"  Difference: {diffs[best]:.6e}")
        print(f"  Ratio: {(factor_astropy_per_pixel - diffs[best]) / factor_astropy_per_pixel if factor_astropy_per_pixel != 0 else np.nan:.4f}")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    analyze_formulas()
