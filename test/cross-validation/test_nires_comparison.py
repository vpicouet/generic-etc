#!/usr/bin/env python3
"""
Compare Generic ETC with Keck NIRES ETC
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

from Observation import Observation, load_instruments

# Import NIRES ETC calculation function directly (avoid pysynphot import issue)
import scipy.constants as sc
from astropy.io import ascii
import math

# NIRES datadir
nires_datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/nires/"

def calculate_nires_s2n_pointsource(mag_src, obs_wave, exp_time, coadds, dither, dither_repeat, seeing, num_reads):
    """NIRES SNR calculation for point source (copied from etc_nires.py)"""

    dark_current = 0.01  # e-/s
    throughput = ascii.read(nires_datadir + "nires_throughput_18_edit.dat")
    tpt_wave = np.array(throughput.columns[0])
    tpt_val = np.array(throughput.columns[1])
    tpt_val_interp = np.interp(obs_wave, tpt_wave, tpt_val)
    gain = 3.8  # e-/ADU
    pix_size = 0.123  # arcsec/pix
    collect_area = 76 * 1e4  # 76 m² = 76e4 cm²
    del_lmbd = obs_wave / 2700 / 1e4  # cm
    spatial_cov = (0.55 * seeing) / (math.pi * seeing**2)
    num_pix_slt = 18 * 0.55 / pix_size**2
    num_pix_src = 0.55 * seeing

    # Mauna Kea IR sky background
    IR_sky = ascii.read(nires_datadir + "mk_skybg_zm_10_10_ph.dat")
    sky_wave = np.array(IR_sky.columns[0]) / 1000  # μm
    sky_bkg = 1000 * np.array(IR_sky.columns[1])  # photon/s/arcsec²/cm/cm²
    sky_bkg_interp = np.interp(obs_wave, sky_wave, sky_bkg)

    flux_src = np.power(10, -0.4 * mag_src) / 1e8  # input source flux in erg/s/cm²/cm

    photon_energy = sc.h * 1e7 * (sc.c * 1e2) / (obs_wave * 1e4)  # single photon energy in erg

    flux_src_phot = flux_src / photon_energy  # source flux in photon/s/cm²/cm

    sig_src = flux_src_phot * collect_area * tpt_val_interp * del_lmbd * spatial_cov  # source signal in e-/s

    if num_reads == 2:
        read_noise = 15  # CDS
    elif num_reads == 16:
        read_noise = 5  # MCDS

    if dither == "AB":
        num_dith = 2
    elif dither == "ABBA":
        num_dith = 4

    sig_src_int = sig_src * exp_time  # source flux integrated over time, single frame
    noise_int = (sig_src * exp_time + num_pix_slt * dark_current * exp_time +
                 (num_pix_slt - num_pix_src) * sky_bkg_interp * exp_time +
                 num_pix_slt * read_noise**2)**(1/2)  # noise, single frame

    s2n = sig_src_int / noise_int * math.sqrt(num_dith)  # s/n

    return s2n, sig_src_int, noise_int, tpt_val_interp, sky_bkg_interp

print("✓ Imports OK")


def compare_nires():
    """
    Compare Generic ETC with Keck NIRES ETC
    """
    print("=" * 80)
    print("COMPARISON: Generic ETC vs Keck NIRES ETC")
    print("=" * 80)

    # Load Generic ETC instruments
    instruments, database = load_instruments()
    instrument_name = "NIRES"

    # Extract Generic ETC parameters
    params = {}
    for i, charact in enumerate(instruments["Charact."]):
        if charact and not isinstance(charact, np.ma.core.MaskedConstant):
            value = instruments[instrument_name][i]
            if not isinstance(value, np.ma.core.MaskedConstant):
                params[charact] = value

    print("\n" + "-" * 80)
    print("GENERIC ETC PARAMETERS (NIRES)")
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
        "Spectral_resolution",
    ]:
        if key in params:
            print(f"  {key:20s} = {params[key]}")

    # Test parameters
    t_exp = params.get('exposure_time', 300.0)
    if t_exp == 0 or np.ma.is_masked(t_exp):
        t_exp = 300.0  # Default exposure time

    readout_time = params.get("Readout", 10.0)
    acq_time = (t_exp + readout_time) / 3600.0

    # NIRES parameters for Keck ETC
    obs_wave = params["wavelength"] / 1000  # nm -> μm
    seeing = 0.75  # arcsec

    print("\n" + "-" * 80)
    print("TEST SETUP")
    print("-" * 80)
    print(f"  Exposure time: {t_exp} s")
    print(f"  Seeing: {seeing} arcsec")
    print(f"  Wavelength: {obs_wave} μm = {params['wavelength']} nm")

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
        IFS=False,
        test=True,
    )

    # Integrated values
    signal_integrated = obs.Signal_el[obs.i] * obs.pixels_total_source
    sky_integrated = obs.sky[obs.i] * obs.pixels_total_source
    dark_integrated = obs.Dark_current_f[obs.i] * obs.pixels_total_source
    rn2_integrated = (obs.RN_final[obs.i] ** 2) * obs.pixels_total_source

    print(f"\nGeneric ETC Results:")
    print(f"  Wavelength: {params['wavelength']} nm")
    print(f"\n  === Per pixel values ===")
    print(f"  Signal_el: {obs.Signal_el[obs.i]:.6e} e⁻/pix")
    print(f"  sky: {obs.sky[obs.i]:.6e} e⁻/pix")
    print(f"  Dark_current_f: {obs.Dark_current_f[obs.i]:.6e} e⁻/pix")
    print(f"  RN_final: {obs.RN_final[obs.i]:.6e} e⁻/pix")
    print(f"  SNR per pix: {obs.SNR[obs.i]:.6e}")
    print(f"\n  pixels_total_source: {obs.pixels_total_source:.3f}")

    print(f"\n  === Integrated values ===")
    print(f"  Signal integrated: {signal_integrated:.6e} e⁻")
    print(f"  Sky integrated: {sky_integrated:.6e} e⁻")
    print(f"  Dark integrated: {dark_integrated:.6e} e⁻")
    print(f"  RN² integrated: {rn2_integrated:.6e} e⁻²")

    # Calculate integrated SNR
    snr_integrated = signal_integrated / np.sqrt(signal_integrated + sky_integrated + dark_integrated + rn2_integrated)
    print(f"  SNR integrated: {snr_integrated:.6e}")

    print("\n" + "=" * 80)
    print("2. KECK NIRES ETC CALCULATION")
    print("=" * 80)

    # NIRES ETC parameters
    # Convert signal flux to magnitude
    # flux in erg/cm²/s/arcsec²/Å -> magnitude
    flux_erg = params["Signal"]  # erg/cm²/s/arcsec²/Å

    # For point source, integrate over seeing disk to get total flux
    # flux_total in erg/cm²/s/Å
    flux_total = flux_erg * seeing**2  # erg/cm²/s/Å

    # Convert to erg/cm²/s/cm (NIRES uses this unit)
    # 1 Å = 1e-8 cm, so flux_per_Å * 1e8 = flux_per_cm
    flux_per_cm = flux_total * 1e8  # erg/cm²/s/Å -> erg/cm²/s/cm

    # NIRES ETC uses: flux_src = 10^(-0.4*mag) / 1e8 [erg/s/cm²/cm]
    # So: 10^(-0.4*mag) = flux_per_cm * 1e8
    # mag = -2.5 * log10(flux_per_cm * 1e8)
    mag_src = -2.5 * np.log10(flux_per_cm * 1e8)

    print(f"\nFlux conversion:")
    print(f"  Signal flux: {flux_erg:.3e} erg/cm²/s/arcsec²/Å")
    print(f"  Integrated (seeing²={seeing**2:.2f} arcsec²): {flux_total:.3e} erg/cm²/s/Å")
    print(f"  In erg/cm²/s/cm: {flux_per_cm:.3e}")
    print(f"  Converted magnitude: {mag_src:.2f}")

    # Verify: this is extremely faint, should be ~25-30 mag, not 6!
    # The issue is that 5.6e-19 * 0.56 * 1e8 = 3.15e-11
    # 10^(-0.4*mag) = 3.15e-11 * 1e8 = 3.15e-3
    # mag = -2.5 * log10(3.15e-3) = 6.25
    # But this is WRONG - we need AB magnitude

    # AB magnitude: m_AB = -2.5 log10(f_nu) - 48.6
    # f_nu = f_lambda * lambda^2 / c
    # At 1.5 um = 1.5e-4 cm, c = 3e10 cm/s
    lambda_cm = obs_wave * 1e-4  # μm -> cm
    c_cgs = 3e10  # cm/s
    f_nu = flux_per_cm * lambda_cm**2 / c_cgs  # erg/cm²/s/Hz
    mag_AB = -2.5 * np.log10(f_nu) - 48.6

    print(f"  f_nu: {f_nu:.3e} erg/cm²/s/Hz")
    print(f"  AB magnitude: {mag_AB:.2f}")

    # Use AB magnitude
    mag_src = mag_AB

    print(f"\nNIRES ETC Input:")
    print(f"  Wavelength: {obs_wave} μm")
    print(f"  Exposure time: {t_exp} s")
    print(f"  Seeing: {seeing} arcsec")
    print(f"  Test magnitude: {mag_src}")

    # Call NIRES ETC
    try:
        snr_nires, sig_nires, noise_nires, tpt_nires, sky_nires = calculate_nires_s2n_pointsource(
            mag_src=mag_src,
            obs_wave=obs_wave,
            exp_time=t_exp,
            coadds=1,
            dither="AB",
            dither_repeat=1,
            seeing=seeing,
            num_reads=16  # MCDS for lower read noise
        )
        print(f"\nNIRES ETC Result:")
        print(f"  Signal (single frame): {sig_nires:.6e} e⁻")
        print(f"  Noise (single frame): {noise_nires:.6e} e⁻")
        print(f"  Throughput at {obs_wave} μm: {tpt_nires:.4f}")
        print(f"  Sky background: {sky_nires:.6e} photons/s/arcsec²/cm/cm²")
        print(f"  SNR (with AB dither): {snr_nires:.6e}")
    except Exception as e:
        print(f"\nNIRES ETC Error: {e}")
        import traceback
        traceback.print_exc()
        snr_nires = None
        sig_nires = None

    print("\n" + "=" * 80)
    print("3. COMPARISON")
    print("=" * 80)

    if snr_nires:
        print(f"\nSNR Comparison:")
        print(f"  Generic ETC: {snr_integrated:.6e}")
        print(f"  NIRES ETC: {snr_nires:.6e}")
        print(f"  Ratio (Generic/NIRES): {snr_integrated / snr_nires:.4f}")

    # Key parameters comparison
    print(f"\nKey Parameters:")
    print(f"  Generic - Collecting area: {params['Collecting_area']} m²")
    print(f"  NIRES - Collecting area: 76 m²")
    print(f"  Generic - Throughput: {params['Throughput']}")
    print(f"  Generic - QE: {params['QE']}")
    print(f"  Generic - Atmosphere: {params['Atmosphere']}")
    print(f"  Total efficiency: {params['Throughput'] * params['QE'] * params['Atmosphere']:.4f}")

    return {
        "generic_signal": signal_integrated,
        "generic_sky": sky_integrated,
        "generic_snr": snr_integrated,
        "nires_snr": snr_nires,
    }


if __name__ == "__main__":
    result = compare_nires()
