#!/usr/bin/env python3
"""
Direct comparison of electron counts between Generic and KCWI
"""

import sys
import os
import numpy as np

sys.path.insert(0, "/Users/Vincent/Github/generic-etc/notebooks")
sys.path.insert(0, "/Users/Vincent/Github/exposureTimeCalculator")

from Observation import load_instruments, Observation
import etc_kcwi
from etc_kcwi import calculate_snr as kcwi_calculate_snr, get_params, wpA, obj_cts, sky_cts

etc_kcwi.datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/kcwi/"

# Parameters
slicer = "Small"
grating = "BL"
gratwave = 4550
seeing = 0.75
t_exp = 3600
flux_erg = 5.6e-19
read_noise = 2.7

print("="*80)
print("KCWI Counts Calculation")
print("="*80)

w = np.arange(3000.0, 6000.0, 1.0)
idx = np.argmin(np.abs(w - gratwave))

# Spatial/spectral bins
snr_spatial_bin = seeing * seeing
pixels_spectral = 2
A_per_pixel_kcwi = 0.625
snr_spectral_bin = pixels_spectral * A_per_pixel_kcwi
arcsec_per_slice = 0.35
nslices = seeing / arcsec_per_slice
pixels_per_arcsec = 1.0 / 0.147
pixels_spat_bin = pixels_per_arcsec * nslices

# Get counts
_, pA = wpA(flux_erg)
c_o_array = obj_cts(w, pA, grating, t_exp)
c_s_array = sky_cts(w, grating, t_exp)

c_o_per_unit = c_o_array[idx]
c_s_per_unit = c_s_array[idx]

c_o_total = c_o_per_unit * snr_spatial_bin * snr_spectral_bin
c_s_total = c_s_per_unit * snr_spatial_bin * snr_spectral_bin
c_r_total = (read_noise**2) * pixels_spectral * pixels_spat_bin

snr_kcwi = c_o_total / np.sqrt(c_o_total + c_s_total + c_r_total)

print(f"\nObject counts:")
print(f"  c_o (per arcsec² per Å) = {c_o_per_unit:.3f} e-")
print(f"  × spatial bin ({snr_spatial_bin:.4f} arcsec²) = {c_o_per_unit * snr_spatial_bin:.3f}")
print(f"  × spectral bin ({snr_spectral_bin:.2f} Å) = {c_o_total:.3f} e-")

print(f"\nSky counts:")
print(f"  c_s (per arcsec² per Å) = {c_s_per_unit:.3f} e-")
print(f"  × spatial bin ({snr_spatial_bin:.4f} arcsec²) = {c_s_per_unit * snr_spatial_bin:.3f}")
print(f"  × spectral bin ({snr_spectral_bin:.2f} Å) = {c_s_total:.3f} e-")

print(f"\nRead noise²:")
print(f"  {read_noise}² × {pixels_spectral} × {pixels_spat_bin:.3f} = {c_r_total:.3f} e²")

print(f"\nSNR = {c_o_total:.3f} / sqrt({c_o_total:.3f} + {c_s_total:.3f} + {c_r_total:.3f})")
print(f"    = {c_o_total:.3f} / {np.sqrt(c_o_total + c_s_total + c_r_total):.3f}")
print(f"    = {snr_kcwi:.3f}")

# Now Generic ETC
print("\n" + "="*80)
print("Generic ETC Counts Calculation")
print("="*80)

# Get KCWI sky in erg
H = 6.62607015e-27
C = 2.99792458e10
A_GEO_KECK = np.pi / 4.0 * (1000.0) ** 2
E_photon = H * C / (gratwave * 1e-8)
eff_keck = get_params(w, grating)
eff_keck_at_wave = eff_keck[idx]

sky_photons_per_s = c_s_per_unit / (eff_keck_at_wave * A_GEO_KECK)
sky_erg = sky_photons_per_s * E_photon

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
    dispersion=A_per_pixel_kcwi / 10,
    Size_source=seeing,
    Signal=flux_erg,
    Sky=sky_erg,
    RN=read_noise,
    Dark_current=0.0,
    Throughput=eff_keck_at_wave,
    QE=1.0,
    Atmosphere=1.0,
)

gen_signal_total = obs.Signal_el[obs.i] * obs.pixels_total_source
gen_sky_total = obs.sky[obs.i] * obs.pixels_total_source
gen_rn2_total = obs.RN_final[obs.i]**2 * obs.number_pixels_used
snr_generic = obs.SNR[obs.i]

print(f"\nObject counts:")
print(f"  Signal_el[i] = {obs.Signal_el[obs.i]:.3f} e-/pix")
print(f"  × pixels_total_source ({obs.pixels_total_source:.3f}) = {gen_signal_total:.3f} e-")

print(f"\nSky counts:")
print(f"  sky[i] = {obs.sky[obs.i]:.3f} e-/pix")
print(f"  × pixels_total_source ({obs.pixels_total_source:.3f}) = {gen_sky_total:.3f} e-")

print(f"\nRead noise²:")
print(f"  {obs.RN_final[obs.i]}² × {obs.number_pixels_used:.3f} = {gen_rn2_total:.3f} e²")

print(f"\nSNR = {gen_signal_total:.3f} / sqrt({gen_signal_total:.3f} + {gen_sky_total:.3f} + {gen_rn2_total:.3f})")
print(f"    = {gen_signal_total:.3f} / {np.sqrt(gen_signal_total + gen_sky_total + gen_rn2_total):.3f}")
print(f"    = {snr_generic:.3f}")

# Compare
print("\n" + "="*80)
print("COMPARISON")
print("="*80)
print(f"\nObject signal:")
print(f"  KCWI: {c_o_total:.3f} e-")
print(f"  Generic: {gen_signal_total:.3f} e-")
print(f"  Ratio: {gen_signal_total/c_o_total:.4f}")

print(f"\nSky:")
print(f"  KCWI: {c_s_total:.3f} e-")
print(f"  Generic: {gen_sky_total:.3f} e-")
print(f"  Ratio: {gen_sky_total/c_s_total:.4f}")

print(f"\nRead noise²:")
print(f"  KCWI: {c_r_total:.3f} e²")
print(f"  Generic: {gen_rn2_total:.3f} e²")
print(f"  Ratio: {gen_rn2_total/c_r_total:.4f}")

print(f"\nSNR:")
print(f"  KCWI: {snr_kcwi:.3f}")
print(f"  Generic: {snr_generic:.3f}")
print(f"  Ratio: {snr_generic/snr_kcwi:.4f} ({100*(snr_generic/snr_kcwi-1):.1f}%)")
