#!/usr/bin/env python3
"""
Debug: Compare sky calculation between Generic ETC and KCWI ETC
"""

import sys
import os
import numpy as np

# Setup paths
sys.path.insert(0, "/Users/Vincent/Github/generic-etc/notebooks")
sys.path.insert(0, "/Users/Vincent/Github/exposureTimeCalculator")

import etc_kcwi
etc_kcwi.datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/kcwi/"

from etc_kcwi import wpA, obj_cts, sky_cts, get_params

# Test parameters
slicer = "Small"
grating = "BL"
gratwave = 4550
seeing = 0.75
t_exp = 3600
flux_erg = 5.6e-19
sky_erg = 8e-18  # Nominal
read_noise = 2.7

print("="*80)
print("KCWI ETC Internal Calculation Breakdown")
print("="*80)

# KCWI's calculation (from their code)
w = np.arange(3000.0, 6000.0, 1.0)

# Spatial and spectral bins
seeing2 = seeing  # Small slicer
snr_spatial_bin = seeing * seeing2
pixels_spectral = 2
arcsec_per_slice = 0.35
nslices = seeing / arcsec_per_slice

A_per_pixel_kcwi = 0.625
snr_spectral_bin = pixels_spectral * A_per_pixel_kcwi

pixels_per_arcsec = 1.0 / 0.147
pixels_spat_bin = pixels_per_arcsec * nslices
pixels_per_snr_specbin = snr_spectral_bin / A_per_pixel_kcwi

print(f"\nSpatial configuration:")
print(f"  seeing = {seeing}\"")
print(f"  seeing2 = {seeing2}\" (Small slicer uses seeing)")
print(f"  snr_spatial_bin = {seeing} × {seeing2} = {snr_spatial_bin} arcsec²")
print(f"  arcsec_per_slice = {arcsec_per_slice}\"")
print(f"  nslices = {seeing} / {arcsec_per_slice} = {nslices:.3f}")
print(f"  pixels_spat_bin = {pixels_per_arcsec:.3f} pix/\" × {nslices:.3f} = {pixels_spat_bin:.3f} pixels")

print(f"\nSpectral configuration:")
print(f"  pixels_spectral = {pixels_spectral}")
print(f"  A_per_pixel = {A_per_pixel_kcwi} Å/pix")
print(f"  snr_spectral_bin = {pixels_spectral} × {A_per_pixel_kcwi} = {snr_spectral_bin} Å")
print(f"  pixels_per_snr_specbin = {snr_spectral_bin} / {A_per_pixel_kcwi} = {pixels_per_snr_specbin:.1f} pixels")

# Get photons/Angstrom
_, pA = wpA(flux_erg)

# Calculate counts
c_o_array = obj_cts(w, pA, grating, t_exp)
c_s_array = sky_cts(w, grating, t_exp)

idx = np.argmin(np.abs(w - gratwave))

c_o_per_A = c_o_array[idx]
c_s_per_A = c_s_array[idx]

c_o = c_o_per_A * snr_spatial_bin * snr_spectral_bin
c_s = c_s_per_A * snr_spatial_bin * snr_spectral_bin

bin_factor = 1.0
nframes = 1.0
c_r = nframes * (read_noise**2) * pixels_per_snr_specbin * pixels_spat_bin * bin_factor

snr_kcwi = c_o / np.sqrt(c_s + c_o + c_r)

print(f"\nCounts at {gratwave}Å:")
print(f"  c_o (per Å per arcsec²) = {c_o_per_A:.3f} e-")
print(f"  c_o × spatial_bin × spectral_bin = {c_o_per_A:.3f} × {snr_spatial_bin} × {snr_spectral_bin} = {c_o:.3f} e-")
print(f"  c_s (per Å per arcsec²) = {c_s_per_A:.3f} e-")
print(f"  c_s × spatial_bin × spectral_bin = {c_s_per_A:.3f} × {snr_spatial_bin} × {snr_spectral_bin} = {c_s:.3f} e-")
print(f"  c_r = {read_noise}² × {pixels_per_snr_specbin} × {pixels_spat_bin:.3f} × {bin_factor} = {c_r:.3f} e²")

print(f"\nSNR calculation:")
print(f"  Noise² = signal + sky + rn² = {c_o:.3f} + {c_s:.3f} + {c_r:.3f} = {c_o + c_s + c_r:.3f}")
print(f"  SNR = {c_o:.3f} / sqrt({c_o + c_s + c_r:.3f}) = {snr_kcwi:.3f}")

# Now test with different sky values
print("\n" + "="*80)
print("Sky variation test")
print("="*80)

for sky_factor in [0.5, 1.0, 2.0]:
    sky_test = sky_erg * sky_factor

    # KCWI doesn't change anything based on sky_erg input!!
    # It uses its own sky model from sky_mk()
    c_s_test_array = sky_cts(w, grating, t_exp)
    c_s_test = c_s_test_array[idx] * snr_spatial_bin * snr_spectral_bin

    snr_test = c_o / np.sqrt(c_s_test + c_o + c_r)

    print(f"\nsky_erg = {sky_test:.1e}:")
    print(f"  KCWI sky counts = {c_s_test:.3f} e- (SAME AS BEFORE!)")
    print(f"  SNR = {snr_test:.3f}")

print("\n" + "="*80)
print("KEY FINDING:")
print("="*80)
print("KCWI ETC uses its OWN sky model (from mk_sky.dat file)!")
print("It does NOT use the input sky_erg parameter!")
print("This is why varying sky_erg doesn't change KCWI SNR!")
