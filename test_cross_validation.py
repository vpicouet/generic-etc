#!/usr/bin/env python3
"""
Cross-validation script for Generic ETC vs Astropy
Tests various parameters to identify discrepancies
"""

import sys
sys.path.insert(0, 'notebooks')

import numpy as np
from Observation import Observation, load_instruments
from astropy.stats import signal_to_noise_oir_ccd

# Load instruments
instruments, database = load_instruments()

# Test with GALEX FUV (simple case: no read noise)
instrument_name = "GALEX FUV"

# Extract parameters
params = {}
for i, charact in enumerate(instruments['Charact.']):
    if charact and not isinstance(charact, np.ma.core.MaskedConstant):
        value = instruments[instrument_name][i]
        if not isinstance(value, np.ma.core.MaskedConstant):
            params[charact] = value

print("="*70)
print(f"CROSS-VALIDATION: {instrument_name}")
print("="*70)

# Reference case
t_exp = 1000.0  # seconds

# Generic ETC
obs = Observation(
    instruments=instruments,
    instrument=instrument_name,
    exposure_time=t_exp,
    SNR_res="per pix",
    IFS=False,
    test=True
)

print(f"\n{instrument_name} parameters:")
print(f"  Wavelength: {params['wavelength']:.1f} nm")
print(f"  Collecting area: {params['Collecting_area']:.3f} m²")
print(f"  Pixel scale: {params['pixel_scale']:.2f} arcsec/pix")
print(f"  Dispersion: {params['dispersion']:.2f} Å/pix")
print(f"  Dark current: {params['Dark_current']:.4f} e⁻/pix/hour")
print(f"  Signal: {params['Signal']:.2e} erg/cm²/s/arcsec²/Å")
print(f"  Sky: {params['Sky']:.2e} erg/cm²/s/arcsec²/Å")

print(f"\n{'='*70}")
print("GENERIC ETC - Internal values")
print("="*70)
print(f"  pixels_total_source: {obs.pixels_total_source}")
print(f"  elem_size: {obs.elem_size}")
print(f"  number_pixels_used: {obs.number_pixels_used}")
print(f"  flux_fraction: {obs.flux_fraction:.6f}")
print(f"\n  Signal_el (e⁻/pix/exp): {obs.Signal_el[obs.i]:.6e}")
print(f"  sky (e⁻/pix/exp): {obs.sky[obs.i]:.6e}")
print(f"  Dark_current_f (e⁻/pix/exp): {obs.Dark_current_f[obs.i]:.6e}")
print(f"  Sky_noise (e⁻): {obs.Sky_noise[obs.i]:.6e}")
print(f"  SNR: {obs.SNR[obs.i]:.6e}")

# Astropy calculation
wavelength_nm = params['wavelength']
E_photon = 1.986e-8 / wavelength_nm  # erg

signal_flux = params['Signal']  # erg/cm²/s/arcsec²/Å
sky_flux = params['Sky']        # erg/cm²/s/arcsec²/Å

# Convert to photons
signal_photons_rate = signal_flux / E_photon  # photons/cm²/s/arcsec²/Å
sky_photons_rate = sky_flux / E_photon        # photons/cm²/s/arcsec²/Å

# Telescope/instrument
collecting_area_cm2 = params['Collecting_area'] * 1e4  # m² → cm²
pixel_area_arcsec2 = params['pixel_scale']**2          # arcsec²
wavelength_range = params['dispersion']                 # Å/pixel
total_throughput = params['Throughput'] * params['QE'] * params['Atmosphere']

# Calculate electron rates PER PIXEL
source_eps = (signal_photons_rate *
              collecting_area_cm2 *
              pixel_area_arcsec2 *
              wavelength_range *
              total_throughput)

sky_eps = (sky_photons_rate *
           collecting_area_cm2 *
           pixel_area_arcsec2 *
           wavelength_range *
           total_throughput)

dark_eps = params['Dark_current'] / 3600.0  # e⁻/pix/hour → e⁻/s/pix

print(f"\n{'='*70}")
print("ASTROPY - Electron rates (per pixel)")
print("="*70)
print(f"  source_eps: {source_eps:.6e} e⁻/s/pix")
print(f"  sky_eps: {sky_eps:.6e} e⁻/s/pix")
print(f"  dark_eps: {dark_eps:.6e} e⁻/s/pix")

# Total electrons over exposure
n_pixels = 1  # per pixel mode

signal_total = source_eps * t_exp
sky_total = sky_eps * t_exp * n_pixels
dark_total = dark_eps * t_exp * n_pixels

print(f"\n  Total electrons (t_exp = {t_exp}s, npix = {n_pixels}):")
print(f"    Signal: {signal_total:.6e} e⁻")
print(f"    Sky: {sky_total:.6e} e⁻")
print(f"    Dark: {dark_total:.6e} e⁻")

# SNR calculation
snr_astropy = signal_to_noise_oir_ccd(
    t=t_exp,
    source_eps=source_eps,
    sky_eps=sky_eps,
    dark_eps=dark_eps,
    rd=params['RN'],
    npix=n_pixels,
    gain=1.0
)
print(f"  SNR (astropy): {snr_astropy:.6e}")

print(f"\n{'='*70}")
print("COMPARISON: Generic ETC vs Astropy")
print("="*70)

print(f"\nElectrons per pixel per exposure:")
print(f"  Signal: Generic={obs.Signal_el[obs.i]:.6e}, Astropy={signal_total:.6e}")
print(f"         Ratio (Gen/Ast): {obs.Signal_el[obs.i]/signal_total:.4f}")

print(f"  Sky:    Generic={obs.sky[obs.i]:.6e}, Astropy={sky_total:.6e}")
print(f"         Ratio (Gen/Ast): {obs.sky[obs.i]/sky_total:.4f}")

print(f"  Dark:   Generic={obs.Dark_current_f[obs.i]:.6e}, Astropy={dark_total:.6e}")
print(f"         Ratio (Gen/Ast): {obs.Dark_current_f[obs.i]/dark_total:.4f}")

print(f"\nSNR:")
print(f"  Generic: {obs.SNR[obs.i]:.6e}")
print(f"  Astropy: {snr_astropy:.6e}")
print(f"  Ratio (Gen/Ast): {obs.SNR[obs.i]/snr_astropy:.4f}")

# DIAGNOSIS
print(f"\n{'='*70}")
print("DIAGNOSIS")
print("="*70)

# Check if the issue is pixels_total_source division
print(f"\nGeneric ETC uses:")
print(f"  factor_CU2el = ... / pixels_total_source")
print(f"  pixels_total_source = {obs.pixels_total_source}")

# Recalculate what Generic ETC would give WITHOUT pixel division
Signal_el_corrected = obs.Signal_el[obs.i] * obs.pixels_total_source
sky_corrected = obs.sky[obs.i] * obs.pixels_total_source

print(f"\nIf we multiply back by pixels_total_source:")
print(f"  Signal_el_corrected: {Signal_el_corrected:.6e} e⁻")
print(f"  sky_corrected: {sky_corrected:.6e} e⁻")
print(f"  Ratio Signal (corrected/astropy): {Signal_el_corrected/signal_total:.4f}")
print(f"  Ratio Sky (corrected/astropy): {sky_corrected/sky_total:.4f}")

print(f"\n{'='*70}")
print("CONCLUSION")
print("="*70)
if abs(Signal_el_corrected/signal_total - 1) < 0.1:
    print("✓ Signal matches after multiplying by pixels_total_source!")
    print("  → Generic ETC divides by pixels_total_source")
    print("  → Astropy expects values per pixel")
else:
    print("✗ Signal still doesn't match - need further investigation")

if abs(sky_corrected/sky_total - 1) < 0.1:
    print("✓ Sky matches after multiplying by pixels_total_source!")
else:
    print("✗ Sky still doesn't match - need further investigation")
