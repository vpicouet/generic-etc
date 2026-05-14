#!/usr/bin/env python3
"""
Compare Generic ETC with Keck KCWI ETC across all configurations

Tests all valid combinations of:
- Slicers: Large, Medium, Small
- Gratings: BL, BM, BH1, BH2, BH3
- Binning: 1x1, 2x2 (except 2x2 with Small slicer)

MODIFICATIONS:
- Fixed numpy array formatting error: factor_CU2el now indexed with [obs.i]
- Pixel calculation matching: To match Keck ETC's pixel calculation, we use:
  * IFS=True (KCWI is an IFS!)
  * Size_source = seeing (source fills the seeing disk)
  * Slitwidth = seeing (to avoid quadratic explosion in pixels_total_source calculation)
  * Line_width = default (1.25nm, source property)
  Note: Using Slitwidth=seeing instead of actual slicer width (0.35-1.35") is a workaround
  for the Generic ETC formula: pixels_total_source = source_size × ceil(sqrt(Size_source² + PSF²) × 2.35 / Slitwidth)
  which explodes when Slitwidth is small.
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

from Observation import load_instruments, Observation
import etc_kcwi
from etc_kcwi import get_params

# Fix datadir for Keck ETC
etc_kcwi.datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/kcwi/"

# KCWI Configuration parameters
SLICER_PARAMS = {
    "Large": {"arcsec_per_slice": 1.35, "pixels_per_spectral_element": 8},
    "Medium": {"arcsec_per_slice": 0.69, "pixels_per_spectral_element": 4},
    "Small": {"arcsec_per_slice": 0.35, "pixels_per_spectral_element": 2},
}

GRATING_PARAMS = {
    "BL": {"A_per_pixel": 0.5, "wave_range": (3500, 5600)},
    "BM": {"A_per_pixel": 0.25, "wave_range": (3500, 5600)},
    "BH1": {"A_per_pixel": 0.125, "wave_range": (3500, 4550)},
    "BH2": {"A_per_pixel": 0.125, "wave_range": (4000, 5050)},
    "BH3": {"A_per_pixel": 0.125, "wave_range": (4800, 5850)},
}

# Constants
H = 6.62607015e-27  # erg·s
C = 2.99792458e10  # cm/s
A_GEO_KECK = np.pi / 4.0 * (1000.0) ** 2  # cm² (10m telescope)
PIXEL_SCALE = 0.147  # arcsec/pixel
FWHM_SIGMA_RATIO = 2.35  # Conversion ratio from sigma to FWHM (used by Generic ETC)


def test_configuration(slicer, grating, ccdbin, gratwave, seeing, t_exp, flux_erg, sky_erg):
    """
    Test a single KCWI configuration and return comparison results.
    """
    # Skip invalid combinations
    if ccdbin == "2x2" and slicer == "Small":
        return None

    # Get slicer and grating parameters
    slicer_p = SLICER_PARAMS[slicer]
    grating_p = GRATING_PARAMS[grating]

    # Binning factor (from official KCWI ETC: bin_factor = 0.25 for 2x2, 1.0 for 1x1)
    # With on-chip 2x2 binning, 4 pixels → 1 super-pixel, so read noise reduces by factor of 4
    bin_factor = 0.25 if ccdbin == "2x2" else 1.0

    # Spatial binning
    arcsec_per_slice = slicer_p["arcsec_per_slice"]
    nslices = seeing / arcsec_per_slice
    snr_spatial_bin = seeing * seeing  # arcsec² for point source

    # Spectral binning
    A_per_pixel = grating_p["A_per_pixel"]
    pixels_spectral = slicer_p["pixels_per_spectral_element"]
    snr_spectral_bin = pixels_spectral * A_per_pixel  # Å

    # Wavelength array
    w = np.arange(3000.0, 6000.0, 1.0)
    idx = np.argmin(np.abs(w - gratwave))

    # Photon energy
    E_photon = H * C / (gratwave * 1e-8)

    # Convert flux to photons
    pA_value = flux_erg / E_photon
    sky_photons = sky_erg / E_photon

    # Keck efficiency
    eff_keck = get_params(w, grating)[idx]

    # Skip if efficiency is 0 (wavelength out of range)
    if eff_keck <= 0:
        return None

    # Object counts (Keck)
    obj_cts_keck = eff_keck * A_GEO_KECK * t_exp * pA_value * snr_spatial_bin * snr_spectral_bin

    # Sky counts (Keck)
    sky_cts_keck = eff_keck * A_GEO_KECK * t_exp * sky_photons * snr_spatial_bin * snr_spectral_bin

    # Read noise contribution
    read_noise = 2.7
    pixels_per_arcsec = 1.0 / PIXEL_SCALE
    pixels_spat_bin = pixels_per_arcsec * nslices
    pixels_per_snr_specbin = snr_spectral_bin / A_per_pixel
    rn_contribution = (read_noise**2) * pixels_per_snr_specbin * pixels_spat_bin * bin_factor

    # Keck SNR
    snr_keck = obj_cts_keck / np.sqrt(sky_cts_keck + obj_cts_keck + rn_contribution)

    # Now compute Generic ETC with matching parameters
    instruments, database = load_instruments()
    instrument_name = "KCWI blue"

    # Run Generic ETC with parameters passed directly to Observation
    # This is the correct way to modify parameters dynamically
    acq_time = t_exp / 3600.0

    # For comparison with Keck, we need to match their pixel calculation
    # Keck uses: pixels = (seeing/pixel_scale) × nslices × pixels_spectral
    # To avoid the quadratic addition in Generic ETC exploding with small Slitwidth,
    # we use seeing as Slitwidth (so sqrt(Size_source² + PSF_RMS²) / Slitwidth doesn't explode)

    # Use KCWI's actual efficiency for this grating (to isolate SNR calculation differences)
    # Generic ETC expects total_efficiency = Throughput × QE × Atmosphere
    # We set Throughput = eff_keck and QE = Atmosphere = 1.0
    obs = Observation(
        instruments=instruments,
        instrument=instrument_name,
        exposure_time=t_exp,
        acquisition_time=acq_time,
        SNR_res="per Source per 2λpix",  # New mode: spatially integrated, 2 spectral pixels
        IFS=True,  # KCWI is an IFS!
        test=True,
        # Pass modified parameters directly
        wavelength=gratwave / 10,  # Å to nm
        Slitwidth=arcsec_per_slice,  # Use actual slicer width now!
        dispersion=A_per_pixel / 10,  # Å/pix to nm/pix
        Size_source=seeing,  # Source fills the seeing disk
        # Line_width: use default (1.25nm) - this is a SOURCE property, not instrument bin
        Signal=flux_erg,  # erg/cm²/s/arcsec²/Å
        Sky=sky_erg,  # erg/cm²/s/arcsec²/Å
        # Force KCWI efficiency to match exactly
        Throughput=eff_keck,  # Use KCWI's empirical efficiency
        QE=1.0,  # Rolled into Throughput
        Atmosphere=1.0,  # Rolled into Throughput
    )

    # Debug: check if parameters were applied
    if grating == "BL" and slicer == "Small" and ccdbin == "1x1":
        print(f"DEBUG [{slicer}-{grating}]: dispersion={obs.dispersion:.4f}, expected={A_per_pixel/10:.4f}")
        print(f"DEBUG: Slitwidth={obs.Slitwidth:.4f}, expected={arcsec_per_slice:.4f}")
        print(f"DEBUG: wavelength={obs.wavelength:.1f}, expected={gratwave/10:.1f}")
        print(f"DEBUG: pixels_total_source={obs.pixels_total_source:.2f}, keck_pixels={pixels_per_snr_specbin * pixels_spat_bin:.2f}")
        print(f"DEBUG: Size_source={obs.Size_source:.4f}")
        print(f"DEBUG: Throughput={obs.Throughput:.3f}, QE={obs.QE:.3f}, Atm={obs.Atmosphere:.3f}")
        print(f"DEBUG: Total_eff={obs.Throughput*obs.QE*obs.Atmosphere:.4f}, keck_eff={eff_keck:.4f}")
        print(f"DEBUG: Collecting_area={obs.Collecting_area:.2f} m², keck={A_GEO_KECK/1e4:.2f} m²")
        print(f"DEBUG: Signal_el/pix={obs.Signal_el[obs.i]:.4e}, integrated={obs.Signal_el[obs.i]*obs.pixels_total_source:.4e}")
        print(f"DEBUG: Sky_el/pix={obs.sky[obs.i]:.4e}, integrated={obs.sky[obs.i]*obs.pixels_total_source:.4e}")
        print(f"DEBUG: RN_final={obs.RN_final[obs.i]:.4e}, RN²×pix={(obs.RN_final[obs.i]**2)*obs.pixels_total_source:.4e}")
        print(f"DEBUG: number_pixels_used={obs.number_pixels_used:.2f} (for SNR calculation)")
        print(f"DEBUG: keck_obj_cts={obj_cts_keck:.4e}, keck_sky_cts={sky_cts_keck:.4e}, keck_rn²={rn_contribution:.4e}")
        print(f"DEBUG: source_size_arcsec_after_slit={obs.source_size_arcsec_after_slit:.4f}")
        print(f"DEBUG: Line_width={obs.Line_width:.4f} nm, Bandwidth={obs.Bandwidth:.4f} nm")
        print(f"DEBUG: factor_CU2el={obs.factor_CU2el[obs.i]:.6e}")
        print(f"DEBUG: keck_spatial={snr_spatial_bin:.4f} arcsec², keck_spectral={snr_spectral_bin:.4f} Å")
        # Calculate what Keck uses
        keck_total_arcsec2_A = snr_spatial_bin * snr_spectral_bin
        generic_total = obs.source_size_arcsec_after_slit * min(obs.Line_width*10, obs.Bandwidth*10)
        print(f"DEBUG: keck_total_arcsec²·Å={keck_total_arcsec2_A:.4f}, generic_total={generic_total:.4f}")

    # For "per Source per 2λpix" mode:
    # - Signal/Sky are integrated over pixels_total_source (full spatial × spectral surface)
    # - But RN/Dark use number_pixels_used (limited spectral bin for extraction)
    if obs.SNR_res == "per Source per 2λpix":
        signal_generic = obs.Signal_el[obs.i] * obs.pixels_total_source
        sky_generic = obs.sky[obs.i] * obs.pixels_total_source
        dark_generic = obs.Dark_current_f[obs.i] * obs.number_pixels_used
        rn2_generic = (obs.RN_final[obs.i] ** 2) * obs.number_pixels_used
    else:
        # Standard integration: all over pixels_total_source
        signal_generic = obs.Signal_el[obs.i] * obs.pixels_total_source
        sky_generic = obs.sky[obs.i] * obs.pixels_total_source
        dark_generic = obs.Dark_current_f[obs.i] * obs.pixels_total_source
        rn2_generic = (obs.RN_final[obs.i] ** 2) * obs.pixels_total_source

    # Calculate SNR manually to match KCWI formula
    snr_generic = signal_generic / np.sqrt(signal_generic + sky_generic + dark_generic + rn2_generic)

    return {
        "config": f"{slicer}-{grating}-{ccdbin}",
        "slicer": slicer,
        "grating": grating,
        "binning": ccdbin,
        "gratwave": gratwave,
        "signal_generic": signal_generic,
        "signal_keck": obj_cts_keck,
        "signal_ratio": signal_generic / obj_cts_keck if obj_cts_keck > 0 else 0,
        "sky_generic": sky_generic,
        "sky_keck": sky_cts_keck,
        "sky_ratio": sky_generic / sky_cts_keck if sky_cts_keck > 0 else 0,
        "snr_generic": snr_generic,
        "snr_keck": snr_keck,
        "snr_ratio": snr_generic / snr_keck if snr_keck > 0 else 0,
        "eff_keck": eff_keck,
        "pixels_generic": obs.pixels_total_source,
        "pixels_keck": pixels_per_snr_specbin * pixels_spat_bin,
    }


def run_all_tests():
    """
    Run tests for all valid KCWI configurations.
    """
    print("=" * 100)
    print("KCWI ETC COMPARISON - ALL CONFIGURATIONS")
    print("=" * 100)

    # Test parameters
    t_exp = 3600.0  # seconds
    seeing = 0.75  # arcsec
    flux_erg = 5.6e-19  # erg/cm²/s/arcsec²/Å
    sky_erg = 8e-18  # erg/cm²/s/arcsec²/Å

    results = []

    # Test all combinations
    slicers = ["Large", "Medium", "Small"]
    gratings = ["BL", "BM", "BH1", "BH2", "BH3"]
    binnings = ["1x1"]  # Generic ETC does not model CCD binning, so only test 1x1

    for slicer in slicers:
        for grating in gratings:
            # Use appropriate wavelength for grating
            wave_range = GRATING_PARAMS[grating]["wave_range"]
            gratwave = (wave_range[0] + wave_range[1]) / 2

            for ccdbin in binnings:
                result = test_configuration(
                    slicer, grating, ccdbin, gratwave, seeing, t_exp, flux_erg, sky_erg
                )
                if result:
                    results.append(result)

    # Print results table
    print("\n" + "-" * 100)
    print(f"{'Config':<20} {'λ(Å)':<8} {'Signal Ratio':<14} {'Sky Ratio':<12} {'SNR Ratio':<12} {'SNR Gen':<12} {'SNR Keck':<12}")
    print("-" * 100)

    for r in results:
        print(
            f"{r['config']:<20} {r['gratwave']:<8.0f} {r['signal_ratio']:<14.3f} "
            f"{r['sky_ratio']:<12.3f} {r['snr_ratio']:<12.3f} {r['snr_generic']:<12.2f} {r['snr_keck']:<12.2f}"
        )

    # Summary statistics
    print("\n" + "=" * 100)
    print("SUMMARY")
    print("=" * 100)

    signal_ratios = [r["signal_ratio"] for r in results]
    sky_ratios = [r["sky_ratio"] for r in results]
    snr_ratios = [r["snr_ratio"] for r in results]

    print(f"\nSignal Ratio (Generic/Keck):")
    print(f"  Mean: {np.mean(signal_ratios):.3f}")
    print(f"  Std:  {np.std(signal_ratios):.3f}")
    print(f"  Min:  {np.min(signal_ratios):.3f}")
    print(f"  Max:  {np.max(signal_ratios):.3f}")

    print(f"\nSky Ratio (Generic/Keck):")
    print(f"  Mean: {np.mean(sky_ratios):.3f}")
    print(f"  Std:  {np.std(sky_ratios):.3f}")
    print(f"  Min:  {np.min(sky_ratios):.3f}")
    print(f"  Max:  {np.max(sky_ratios):.3f}")

    print(f"\nSNR Ratio (Generic/Keck):")
    print(f"  Mean: {np.mean(snr_ratios):.3f}")
    print(f"  Std:  {np.std(snr_ratios):.3f}")
    print(f"  Min:  {np.min(snr_ratios):.3f}")
    print(f"  Max:  {np.max(snr_ratios):.3f}")

    # Validation
    print("\n" + "=" * 100)
    print("VALIDATION")
    print("=" * 100)

    all_within_20 = all(0.8 <= r <= 1.2 for r in snr_ratios)
    all_within_10 = all(0.9 <= r <= 1.1 for r in snr_ratios)

    if all_within_10:
        print("✓ All configurations within 10% agreement")
    elif all_within_20:
        print("✓ All configurations within 20% agreement")
    else:
        outliers = [r for r in results if r["snr_ratio"] < 0.8 or r["snr_ratio"] > 1.2]
        print(f"✗ {len(outliers)} configurations outside 20% agreement:")
        for o in outliers:
            print(f"  - {o['config']}: SNR ratio = {o['snr_ratio']:.3f}")

    return results


def test_single_config():
    """
    Test a single configuration matching test_kcwi_comparison.py
    """
    print("=" * 100)
    print("SINGLE CONFIGURATION TEST (matching test_kcwi_comparison.py)")
    print("=" * 100)

    # Exact same parameters as test_kcwi_comparison.py
    t_exp = 3600.0
    seeing = 0.75
    flux_erg = 5.6e-19
    sky_erg = 8e-18

    slicer = "Small"
    grating = "BL"
    gratwave = 4500.0  # Exact wavelength from comparison.py
    ccdbin = "1x1"

    result = test_configuration(
        slicer, grating, ccdbin, gratwave, seeing, t_exp, flux_erg, sky_erg
    )

    if result:
        print("\n" + "-" * 100)
        print("RESULT:")
        print("-" * 100)
        print(f"  Config: {result['config']}")
        print(f"  Wavelength: {result['gratwave']} Å")
        print(f"  Signal ratio (Generic/Keck): {result['signal_ratio']:.4f}")
        print(f"  Sky ratio (Generic/Keck): {result['sky_ratio']:.4f}")
        print(f"  SNR ratio (Generic/Keck): {result['snr_ratio']:.4f}")
        print(f"  SNR Generic: {result['snr_generic']:.4f}")
        print(f"  SNR Keck: {result['snr_keck']:.4f}")
        print("\n" + "-" * 100)
        if 0.8 <= result['snr_ratio'] <= 1.2:
            print("✓ SNR ratio within 20% agreement")
        else:
            print(f"✗ SNR ratio outside 20% agreement: {result['snr_ratio']:.4f}")

    return result


if __name__ == "__main__":
    # First test single config matching test_kcwi_comparison.py
    print("\n" + "=" * 100)
    print("STEP 1: Testing single configuration (Small-BL-1x1 at 4500Å)")
    print("=" * 100)
    single_result = test_single_config()

    # Then run all tests
    print("\n\n" + "=" * 100)
    print("STEP 2: Testing all configurations")
    print("=" * 100)
    results = run_all_tests()
