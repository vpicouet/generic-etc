#!/usr/bin/env python3
"""
Compare Generic ETC with Keck NIRC2 ETC
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

from Observation import Observation, load_instruments

print("✓ Imports OK")


def calculate_nirc2_snr(magnitude, strehl, exp_time, coadds, num_dith, camera, img_filter, num_read):
    """
    NIRC2 SNR calculation (from etc_nirc2.py)
    """
    noise_read = 56.0  # e-, drops as sqrt(reads)
    if num_read > 2:
        noise_read = noise_read / np.sqrt(num_read)

    gain = 4.0  # e-/DN

    # Camera-dependent parameters
    if camera == "wide":
        background = {'J': 0.5, 'H': 4.0, 'K': 5.7, 'Kp': 5.6, 'Lp': 18535, 'Ms': 18535}
        zero_point = {'J': 26.9, 'H': 26.96, 'K': 26.18, 'Kp': 26.30, 'Lp': 25.08, 'Ms': 22.87}
        num_pix = {'J': 12.5, 'H': 12.5, 'K': 12.5, 'Kp': 12.5, 'Lp': 28.4, 'Ms': 38.5}
    elif camera == "narrow":
        background = {'J': 0.5, 'H': 4.0, 'K': 5.7, 'Kp': 5.6, 'Lp': 18535, 'Ms': 18535}
        zero_point = {'J': 26.9, 'H': 26.96, 'K': 26.18, 'Kp': 26.30, 'Lp': 25.08, 'Ms': 22.87}
        num_pix = {'J': 78.95, 'H': 50.2, 'K': 95.2, 'Kp': 95.2, 'Lp': 283.5, 'Ms': 490.8}

    # Background scaling for wide camera
    if camera == "wide":
        bg = background[img_filter] * gain * 16
    elif camera == "narrow":
        bg = background[img_filter] * gain

    num_exp = coadds * num_dith

    # Zero point adjusted for Strehl
    mag_zero = zero_point[img_filter] + 2.5 * np.log10(strehl)

    # Signal calculation
    signal = num_exp * exp_time * np.power(10, 0.4 * (mag_zero - magnitude))

    # Noise calculation
    noise = np.sqrt(
        num_exp * np.square(noise_read) * num_pix[img_filter] +
        num_pix[img_filter] * background[img_filter] * num_exp * exp_time +
        signal
    )

    s2n = signal / noise

    return {
        'snr': s2n,
        'signal': signal / gain,  # in DN
        'noise': noise / gain,
        'num_pix': num_pix[img_filter],
        'background_per_frame': bg * exp_time / gain,
        'zero_point': zero_point[img_filter],
        'mag_zero_strehl': mag_zero,
        'read_noise': noise_read
    }


def compare_nirc2():
    """
    Compare Generic ETC with Keck NIRC2 ETC
    """
    print("=" * 80)
    print("COMPARISON: Generic ETC vs Keck NIRC2 ETC")
    print("=" * 80)

    # Load Generic ETC instruments
    instruments, database = load_instruments()
    instrument_name = "NIRC2"

    # Extract Generic ETC parameters
    params = {}
    for i, charact in enumerate(instruments["Charact."]):
        if charact and not isinstance(charact, np.ma.core.MaskedConstant):
            value = instruments[instrument_name][i]
            if not isinstance(value, np.ma.core.MaskedConstant):
                params[charact] = value

    print("\n" + "-" * 80)
    print("GENERIC ETC PARAMETERS (NIRC2)")
    print("-" * 80)
    for key in [
        "wavelength",
        "Collecting_area",
        "pixel_scale",
        "Throughput",
        "QE",
        "Atmosphere",
        "Signal",
        "Sky",
        "Dark_current",
        "Size_source",
        "RN",
        "Bandwidth",
    ]:
        if key in params:
            print(f"  {key:20s} = {params[key]}")

    # Test parameters
    t_exp = params.get('exposure_time', 60.0)
    if t_exp == 0 or np.ma.is_masked(t_exp):
        t_exp = 60.0

    # NIRC2 parameters
    strehl = 0.4  # K-band Strehl
    coadds = 10
    num_dith = 5
    camera = "narrow"  # Based on pixel_scale 0.01
    img_filter = "K"
    num_read = 16

    print("\n" + "-" * 80)
    print("TEST SETUP")
    print("-" * 80)
    print(f"  Exposure time: {t_exp} s")
    print(f"  Coadds: {coadds}")
    print(f"  Dithers: {num_dith}")
    print(f"  Camera: {camera}")
    print(f"  Filter: {img_filter}")
    print(f"  Strehl: {strehl}")
    print(f"  Number of reads: {num_read}")

    print("\n" + "=" * 80)
    print("1. GENERIC ETC CALCULATION")
    print("=" * 80)

    # Total acquisition for Generic ETC
    total_exp = t_exp * coadds * num_dith
    acq_time = total_exp / 3600.0

    # Run Generic ETC
    # For imager mode, need to pass Throughput_FWHM (filter bandwidth in Å)
    bandwidth_A = params.get("Bandwidth", 40000)  # Filter bandwidth in Å

    obs = Observation(
        instruments=instruments,
        instrument=instrument_name,
        exposure_time=t_exp,
        acquisition_time=acq_time,
        SNR_res="per pix",
        IFS=False,
        test=True,
        Throughput_FWHM=bandwidth_A,
        spectrograph=False,  # This is an imager
    )

    # Debug: print intermediate values
    print(f"\n  === DEBUG: Intermediate values ===")
    print(f"  factor_CU2el: {obs.factor_CU2el:.6e}")
    print(f"  Throughput_FWHM: {bandwidth_A}")
    print(f"  spectro: {obs.spectro}")
    if hasattr(obs, 'pixels_total_source'):
        print(f"  pixels_total_source: {obs.pixels_total_source:.6e}")
    else:
        print(f"  pixels_total_source: NOT DEFINED")

    # For imager, pixels_total_source should be the PSF area
    signal_integrated = obs.Signal_el[obs.i] * obs.pixels_total_source
    sky_integrated = obs.sky[obs.i] * obs.pixels_total_source
    dark_integrated = obs.Dark_current_f[obs.i] * obs.pixels_total_source
    rn2_integrated = (obs.RN_final[obs.i] ** 2) * obs.pixels_total_source

    print(f"\nGeneric ETC Results (single exposure):")
    print(f"  Wavelength: {params['wavelength']} nm")
    print(f"\n  === Per pixel values ===")
    print(f"  Signal_el: {obs.Signal_el[obs.i]:.6e} e⁻/pix")
    print(f"  sky: {obs.sky[obs.i]:.6e} e⁻/pix")
    print(f"  Dark_current_f: {obs.Dark_current_f[obs.i]:.6e} e⁻/pix")
    print(f"  RN_final: {obs.RN_final[obs.i]:.6e} e⁻/pix")
    print(f"  SNR per pix: {obs.SNR[obs.i]:.6e}")
    print(f"\n  pixels_total_source: {obs.pixels_total_source:.3f}")

    print(f"\n  === Integrated values (single exp) ===")
    print(f"  Signal integrated: {signal_integrated:.6e} e⁻")
    print(f"  Sky integrated: {sky_integrated:.6e} e⁻")
    print(f"  Dark integrated: {dark_integrated:.6e} e⁻")
    print(f"  RN² integrated: {rn2_integrated:.6e} e⁻²")

    # SNR for stacked images
    n_images = coadds * num_dith
    snr_integrated = (signal_integrated * n_images) / np.sqrt(
        (signal_integrated + sky_integrated + dark_integrated) * n_images + rn2_integrated * n_images
    )
    print(f"\n  === Stacked ({n_images} images) ===")
    print(f"  Total signal: {signal_integrated * n_images:.6e} e⁻")
    print(f"  SNR integrated: {snr_integrated:.6e}")

    print("\n" + "=" * 80)
    print("2. KECK NIRC2 ETC CALCULATION")
    print("=" * 80)

    # Convert flux to magnitude
    flux_erg = params["Signal"]  # erg/cm²/s/arcsec²/Å
    seeing = params.get("Size_source", 0.05) * 2.35  # FWHM

    # Integrate over PSF
    flux_total = flux_erg * (np.pi * (seeing/2)**2)  # erg/cm²/s/Å

    # Convert to erg/cm²/s/cm
    flux_per_cm = flux_total * 1e8

    # Convert to AB magnitude
    lambda_cm = params["wavelength"] * 1e-7  # nm -> cm
    c_cgs = 3e10
    f_nu = flux_per_cm * lambda_cm**2 / c_cgs
    mag_AB = -2.5 * np.log10(f_nu) - 48.6

    print(f"\nFlux conversion:")
    print(f"  Signal flux: {flux_erg:.3e} erg/cm²/s/arcsec²/Å")
    print(f"  PSF FWHM: {seeing:.3f} arcsec")
    print(f"  Integrated flux: {flux_total:.3e} erg/cm²/s/Å")
    print(f"  AB magnitude: {mag_AB:.2f}")

    # Call NIRC2 ETC
    result = calculate_nirc2_snr(
        magnitude=mag_AB,
        strehl=strehl,
        exp_time=t_exp,
        coadds=coadds,
        num_dith=num_dith,
        camera=camera,
        img_filter=img_filter,
        num_read=num_read
    )

    print(f"\nNIRC2 ETC Input:")
    print(f"  Magnitude: {mag_AB:.2f}")
    print(f"  Zero point (K): {result['zero_point']}")
    print(f"  Zero point with Strehl: {result['mag_zero_strehl']:.2f}")

    print(f"\nNIRC2 ETC Result:")
    print(f"  Signal: {result['signal']:.6e} DN")
    print(f"  Noise: {result['noise']:.6e} DN")
    print(f"  SNR: {result['snr']:.6e}")
    print(f"  Aperture pixels: {result['num_pix']}")
    print(f"  Background/frame: {result['background_per_frame']:.2f} DN")
    print(f"  Read noise: {result['read_noise']:.2f} e⁻")

    print("\n" + "=" * 80)
    print("3. COMPARISON")
    print("=" * 80)

    print(f"\nSNR Comparison:")
    print(f"  Generic ETC: {snr_integrated:.6e}")
    print(f"  NIRC2 ETC: {result['snr']:.6e}")
    if result['snr'] > 0:
        print(f"  Ratio (Generic/NIRC2): {snr_integrated / result['snr']:.4f}")

    print(f"\nKey Parameters:")
    print(f"  Generic - Collecting area: {params['Collecting_area']} m²")
    print(f"  NIRC2 - Collecting area: 78.54 m² (Keck)")
    print(f"  Generic - Throughput: {params['Throughput']}")
    print(f"  Generic - QE: {params['QE']}")
    print(f"  Generic - Atmosphere: {params['Atmosphere']}")
    print(f"  Total efficiency: {params['Throughput'] * params['QE'] * params['Atmosphere']:.4f}")
    print(f"  Generic - pixels: {obs.pixels_total_source:.1f}")
    print(f"  NIRC2 - pixels: {result['num_pix']}")

    return {
        "generic_snr": snr_integrated,
        "nirc2_snr": result['snr'],
    }


if __name__ == "__main__":
    result = compare_nirc2()
