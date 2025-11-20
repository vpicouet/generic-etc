#!/usr/bin/env python3
"""
Compare Generic ETC with Keck KCWI ETC

Tests the conversions for KCWI blue spectrograph
"""

import sys
import os
import numpy as np
from pathlib import Path

# In notebook, use current working directory
script_dir = str(Path.cwd())
script_dir = "/Users/Vincent/Github/generic-etc/test/cross-validation"
notebooks_dir = str(Path.cwd().parent.parent / "notebooks")
os.chdir(notebooks_dir)

# Setup paths
sys.path.insert(0, ".")
sys.path.insert(0, os.path.join(script_dir, "exposureTimeCalculator"))

from Observation import Observation, load_instruments
from astropy.stats import signal_to_noise_oir_ccd

# Import Keck ETC functions
# sys.path.insert(0, os.path.join(script_dir, "exposureTimeCalculator"))
sys.path.insert(0, "/Users/Vincent/Github/exposureTimeCalculator")
sys.path.insert(0, "/Users/Vincent/Github/exposureTimeCalculator/etc_kcwi.py")

import etc_kcwi
from etc_kcwi import obj_cts, sky_cts, wpA, get_params, sky_mk

# Fix datadir for Keck ETC
etc_kcwi.datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/kcwi/"

print("✓ Imports OK")


def compare_kcwi_blue():
    """
    Compare Generic ETC with Keck ETC for KCWI Blue
    """
    print("=" * 80)
    print("COMPARISON: Generic ETC vs Keck KCWI ETC")
    print("=" * 80)

    # Load Generic ETC instruments
    instruments, database = load_instruments()
    instrument_name = "KCWI blue"

    # Extract Generic ETC parameters
    params = {}
    for i, charact in enumerate(instruments["Charact."]):
        if charact and not isinstance(charact, np.ma.core.MaskedConstant):
            value = instruments[instrument_name][i]
            if not isinstance(value, np.ma.core.MaskedConstant):
                params[charact] = value

    print("\n" + "-" * 80)
    print("GENERIC ETC PARAMETERS (KCWI blue)")
    print("-" * 80)
    for key in [
        "wavelength",
        "Collecting_area",
        "pixel_scale",
        "dispersion",
        "Throughput",
        "QE",
        "Atmosphere",
        "Signal",
        "Sky",
        "Dark_current",
        "Size_source",
        "Slitwidth",
        "Slitlength",
        "Line_width",
        "Bandwidth",
        "RN",
    ]:
        if key in params:
            print(f"  {key:20s} = {params[key]}")

    # Test parameters - use from spreadsheet
    t_exp = params.get('exposure_time', 3600.0)
    readout_time = params.get("Readout", 0.0)
    acq_time = (t_exp + readout_time) / 3600.0

    # For Keck KCWI ETC
    slicer = "Small"  # 0.35" slices
    grating = "BL"  # Blue low resolution
    gratwave = 4500.0  # Å
    seeing = 0.75  # arcsec
    ccdbin = "1x1"

    # Flux for comparison - use same as Generic ETC
    # Convert erg/cm²/s/arcsec²/Å to erg/cm²/s/Å (point source)
    flux_erg = params["Signal"]  # erg/cm²/s/arcsec²/Å

    # For point source: integrate over seeing disk
    # flux_total = flux_erg × seeing²
    flux_total_erg = flux_erg * seeing**2

    print("\n" + "-" * 80)
    print("TEST SETUP")
    print("-" * 80)
    print(f"  Exposure time: {t_exp} s")
    print(f"  Seeing: {seeing} arcsec")
    print(f'  Slicer: {slicer} (0.35" slices)')
    print(f"  Grating: {grating}")
    print(f"  Central wavelength: {gratwave} Å")

    print("\n" + "=" * 80)
    print("1. GENERIC ETC CALCULATION")
    print("=" * 80)

    # Run Generic ETC
    obs = Observation(
        instruments=instruments,
        instrument=instrument_name,
        exposure_time=t_exp,
        acquisition_time=acq_time,
        SNR_res="per pix",
        IFS=False,  # Test as slit mode to match Keck binning
        test=True,
    )

    # Find wavelength index closest to gratwave
    # Generic ETC uses wavelength from database
    wave_generic = params["wavelength"]

    # Integrated values (total electrons in aperture)
    signal_integrated = obs.Signal_el[obs.i] * obs.pixels_total_source
    sky_integrated = obs.sky[obs.i] * obs.pixels_total_source
    dark_integrated = obs.Dark_current_f[obs.i] * obs.pixels_total_source
    rn2_integrated = (obs.RN_final[obs.i] ** 2) * obs.pixels_total_source

    print(f"\nGeneric ETC Results:")
    print(f"  Wavelength: {wave_generic} nm = {wave_generic * 10} Å")
    print(f"\n  === Per pixel values ===")
    print(f"  Signal_el: {obs.Signal_el[obs.i]:.6e} e⁻/pix")
    print(f"  sky: {obs.sky[obs.i]:.6e} e⁻/pix")
    print(f"  Dark_current_f: {obs.Dark_current_f[obs.i]:.6e} e⁻/pix")
    print(f"  RN_final: {obs.RN_final[obs.i]:.6e} e⁻/pix")
    print(f"  SNR per pix: {obs.SNR[obs.i]:.6e}")
    print(f"\n  pixels_total_source: {obs.pixels_total_source:.3f}")
    print(f"  number_pixels_used: {obs.number_pixels_used}")

    print(f"\n  === Integrated values (total in aperture) ===")
    print(f"  Signal integrated: {signal_integrated:.6e} e⁻")
    print(f"  Sky integrated: {sky_integrated:.6e} e⁻")
    print(f"  Dark integrated: {dark_integrated:.6e} e⁻")
    print(f"  RN² integrated: {rn2_integrated:.6e} e⁻²")

    # Calculate integrated SNR
    snr_integrated = signal_integrated / np.sqrt(signal_integrated + sky_integrated + dark_integrated + rn2_integrated)
    print(f"  SNR integrated: {snr_integrated:.6e}")

    # Additional debug info
    print(f"\n  === DEBUG: Integrated values ===")
    print(f"  Signal_el × pixels_total_source = {obs.Signal_el[obs.i] * obs.pixels_total_source:.6e} e⁻")
    print(f"  sky × pixels_total_source = {obs.sky[obs.i] * obs.pixels_total_source:.6e} e⁻")

    # Internal variables
    for attr in ['factor_CU2el', 'factor_CU2el_tot', 'factor_CU2el_sky',
                 'source_size_arcsec_after_slit', 'slit_size_arcsec_after_slit',
                 'effective_area', 'Signal_LU', 'Sky_CU', 'flux_fraction_slit_applied']:
        if hasattr(obs, attr):
            val = getattr(obs, attr)
            if hasattr(val, '__getitem__'):
                try:
                    print(f"  {attr}: {val[obs.i]:.6e}")
                except:
                    print(f"  {attr}: {val}")
            else:
                print(f"  {attr}: {val}")

    # Spatial calculations
    fwhm_sigma_ratio = 2.35
    sigma_x_pix = params['Size_source'] / params['pixel_scale'] / fwhm_sigma_ratio
    sigma_y_pix = params['Line_width'] / params['dispersion'] / fwhm_sigma_ratio
    pixels_calc = 2 * np.pi * sigma_x_pix * sigma_y_pix

    print(f"\n  === DEBUG: Pixel calculations ===")
    print(f"  sigma_x_pix: {sigma_x_pix:.3f} pix (spatial)")
    print(f"  sigma_y_pix: {sigma_y_pix:.3f} pix (spectral)")
    print(f"  2π σx σy: {pixels_calc:.3f} pix")

    print("\n" + "=" * 80)
    print("2. KECK KCWI ETC CALCULATION")
    print("=" * 80)

    # Keck telescope: 10m diameter
    A_geo_keck = np.pi / 4.0 * (1000.0) ** 2  # cm²
    print(f"\nKeck telescope:")
    print(f"  Diameter: 10 m")
    print(f"  Geometric area: {A_geo_keck:.2f} cm²")

    # Generic ETC telescope
    A_geo_generic = params["Collecting_area"] * 1e4  # m² → cm²
    print(f"\nGeneric ETC telescope (KCWI blue):")
    print(f"  Collecting area: {params['Collecting_area']:.4f} m²")
    print(f"  = {A_geo_generic:.2f} cm²")
    print(f"  Ratio Generic/Keck: {A_geo_generic / A_geo_keck:.4f}")

    # Use full wavelength array like Keck ETC does
    w = np.arange(3000.0, 6000.0, 1.0)  # Å
    wave_target = wave_generic * 10  # nm → Å
    idx = np.argmin(np.abs(w - wave_target))  # Find closest wavelength

    # Get photons/Angstrom from flux
    # Keck ETC uses: pA = flux/(2.0*10**-8/w)  for linear f power law
    # flux here is in erg/cm²/s/Å
    h = 6.62607015e-27  # erg·s
    c = 2.99792458e10  # cm/s
    E_photon = h * c / (w[idx] * 1e-8)  # erg

    # Convert erg flux to photon flux
    pA_value = flux_erg / E_photon  # photons/cm²/s/Å/arcsec²

    # Create constant photon flux array (for simplicity)
    pA = np.ones_like(w) * pA_value

    print(f"\nFlux conversion:")
    print(f"  E_photon: {E_photon:.6e} erg")
    print(f"  Flux (erg/cm²/s/Å/arcsec²): {flux_erg:.6e}")
    print(f"  Flux (photons/cm²/s/Å/arcsec²): {pA_value:.6e}")

    # Get Keck efficiency at this wavelength
    eff_keck_array = get_params(w, grating)
    eff_keck = eff_keck_array[idx]
    print(f"\nKeck KCWI efficiency:")
    print(f"  Grating: {grating}")
    print(f"  Efficiency (includes throughput × atmosphere): {eff_keck:.4f}")

    # Generic ETC efficiency
    eff_generic = params["Throughput"] * params["QE"] * params["Atmosphere"]
    print(f"\nGeneric ETC efficiency:")
    print(f"  Throughput: {params['Throughput']:.4f}")
    print(f"  QE: {params['QE']:.4f}")
    print(f"  Atmosphere: {params['Atmosphere']:.4f}")
    print(f"  Total: {eff_generic:.4f}")
    print(f"  Ratio Generic/Keck: {eff_generic / eff_keck:.4f}")

    # Calculate object counts with Keck formula
    # Need to account for spatial and spectral binning

    # KCWI Small slicer: 0.35" per slice
    arcsec_per_slice = 0.35
    nslices = seeing / arcsec_per_slice  # number of slices

    # Spectral: Small slicer BL grating = 0.625 Å/pixel, 2 pixels per spectral element
    A_per_pixel = 0.625  # Å/pixel
    pixels_spectral = 2  # for Small slicer
    A_per_specbin = pixels_spectral * A_per_pixel  # Å

    # Spatial bin
    snr_spatial_bin = seeing * seeing  # arcsec² for point source
    snr_spectral_bin = A_per_specbin  # Å

    print(f"\nKCWI binning:")
    print(f'  Slicer: {slicer} ({arcsec_per_slice}" per slice)')
    print(f"  Number of slices: {nslices:.2f}")
    print(f"  Spatial bin: {snr_spatial_bin:.2f} arcsec²")
    print(f"  Spectral bin: {snr_spectral_bin:.2f} Å")
    print(f"  Dispersion: {A_per_pixel} Å/pixel")

    # Object counts (Keck formula)
    # obj_cts multiplies by area and spectral bin outside
    obj_cts_keck_raw = eff_keck * A_geo_keck * t_exp * pA_value  # per arcsec²/Å
    obj_cts_keck = (
        obj_cts_keck_raw * snr_spatial_bin * snr_spectral_bin
    )  # total electrons

    print(f"\nKeck KCWI object counts:")
    print(f"  Raw (per arcsec²/Å): {obj_cts_keck_raw:.6e} e⁻")
    print(f"  With binning: {obj_cts_keck:.6e} e⁻")

    # Sky counts (Keck formula)
    # Use same sky as Generic ETC (convert erg to photons)
    sky_erg = params['Sky']  # erg/cm²/s/arcsec²/Å
    sky_photons_per_arcsec2_per_A = sky_erg / E_photon  # photons/cm²/s/Å/arcsec²
    airmass = 1.0  # already included in Generic ETC atmosphere

    sky_cts_keck_raw = (
        eff_keck * A_geo_keck * t_exp * sky_photons_per_arcsec2_per_A * airmass
    )
    sky_cts_keck = sky_cts_keck_raw * snr_spatial_bin * snr_spectral_bin

    print(f"\nKeck KCWI sky counts:")
    print(f"  Sky flux (photons/cm²/s/Å/arcsec²): {sky_photons_per_arcsec2_per_A:.6e}")
    print(f"  Airmass: {airmass}")
    print(f"  Raw (per arcsec²/Å): {sky_cts_keck_raw:.6e} e⁻")
    print(f"  With binning: {sky_cts_keck:.6e} e⁻")

    # Read noise
    read_noise_keck = 2.7  # e⁻ from Keck ETC
    pixels_per_arcsec = 1.0 / 0.147  # KCWI pixel scale
    pixels_spat_bin = pixels_per_arcsec * nslices
    pixels_per_snr_specbin = snr_spectral_bin / A_per_pixel
    nframes = 1.0
    bin_factor = 1.0  # for 1x1 binning

    rn_contribution_keck = (
        nframes
        * (read_noise_keck**2)
        * pixels_per_snr_specbin
        * pixels_spat_bin
        * bin_factor
    )

    print(f"\nKeck KCWI read noise:")
    print(f"  Read noise: {read_noise_keck} e⁻")
    print(f"  Pixels (spectral): {pixels_per_snr_specbin:.2f}")
    print(f"  Pixels (spatial): {pixels_spat_bin:.2f}")
    print(f"  RN² contribution: {rn_contribution_keck:.6e} e⁻²")

    # Keck SNR
    snr_keck = obj_cts_keck / np.sqrt(
        sky_cts_keck + obj_cts_keck + rn_contribution_keck
    )

    print(f"\nKeck KCWI SNR:")
    print(f"  SNR: {snr_keck:.6e}")

    print("\n" + "=" * 80)
    print("3. COMPARISON (Integrated values)")
    print("=" * 80)

    print(f"\nSignal (electrons):")
    print(f"  Generic ETC: {signal_integrated:.6e} e⁻")
    print(f"  Keck ETC: {obj_cts_keck:.6e} e⁻")
    print(f"  Ratio (Generic/Keck): {signal_integrated / obj_cts_keck:.4f}")

    print(f"\nSky (electrons):")
    print(f"  Generic ETC: {sky_integrated:.6e} e⁻")
    print(f"  Keck ETC: {sky_cts_keck:.6e} e⁻")
    print(f"  Ratio (Generic/Keck): {sky_integrated / sky_cts_keck:.4f}")

    print(f"\nSNR (integrated):")
    print(f"  Generic ETC: {snr_integrated:.6e}")
    print(f"  Keck ETC: {snr_keck:.6e}")
    print(f"  Ratio (Generic/Keck): {snr_integrated / snr_keck:.4f}")

    print("\n" + "=" * 80)
    print("ANALYSIS")
    print("=" * 80)

    # Check if ratios make sense
    telescope_ratio = A_geo_generic / A_geo_keck
    eff_ratio = eff_generic / eff_keck

    print(f"\nExpected ratios:")
    print(f"  Telescope area (Generic/Keck): {telescope_ratio:.4f}")
    print(f"  Efficiency (Generic/Keck): {eff_ratio:.4f}")
    print(f"  Combined (area × eff): {telescope_ratio * eff_ratio:.4f}")

    print(f"\nObserved ratios:")
    print(f"  Signal ratio: {obs.Signal_el[obs.i] / obj_cts_keck:.4f}")
    print(f"  Sky ratio: {obs.sky[obs.i] / sky_cts_keck:.4f}")

    expected_ratio = telescope_ratio * eff_ratio
    signal_ratio = signal_integrated / obj_cts_keck
    sky_ratio = sky_integrated / sky_cts_keck
    snr_ratio = snr_integrated / snr_keck

    print(f"\nValidation:")
    if abs(signal_ratio / expected_ratio - 1) < 0.2:
        print(f"  ✓ Signal ratio matches expected (within 20%)")
    else:
        print(
            f"  ✗ Signal ratio differs from expected by {abs(signal_ratio / expected_ratio - 1) * 100:.1f}%"
        )

    if abs(sky_ratio / expected_ratio - 1) < 0.2:
        print(f"  ✓ Sky ratio matches expected (within 20%)")
    else:
        print(
            f"  ✗ Sky ratio differs from expected by {abs(sky_ratio / expected_ratio - 1) * 100:.1f}%"
        )

    # SNR validation - should scale as sqrt of expected ratio
    expected_snr_ratio = np.sqrt(expected_ratio)
    print(f"\nSNR Validation:")
    print(f"  Expected SNR ratio (sqrt of {expected_ratio:.4f}): {expected_snr_ratio:.4f}")
    print(f"  Actual SNR ratio: {snr_ratio:.4f}")
    if abs(snr_ratio / expected_snr_ratio - 1) < 0.2:
        print(f"  ✓ SNR ratio matches expected (within 20%)")
    else:
        print(
            f"  ✗ SNR ratio differs from expected by {abs(snr_ratio / expected_snr_ratio - 1) * 100:.1f}%"
        )

    return {
        "generic_signal": signal_integrated,
        "keck_signal": obj_cts_keck,
        "generic_sky": sky_integrated,
        "keck_sky": sky_cts_keck,
        "generic_snr": snr_integrated,
        "keck_snr": snr_keck,
        "telescope_ratio": telescope_ratio,
        "eff_ratio": eff_ratio,
        "snr_ratio": snr_ratio,
    }


if __name__ == "__main__":
    result = compare_kcwi_blue()
