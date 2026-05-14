#!/usr/bin/env python3
"""
Fast version: Pre-compute all data, then plot.
Shows SNR ratios as percentage differences from perfect agreement.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import importlib

# Setup paths
script_dir = "/Users/Vincent/Github/generic-etc/test/cross-validation"
notebooks_dir = str(Path(script_dir).parent.parent / "notebooks")
os.chdir(notebooks_dir)

sys.path.insert(0, ".")
sys.path.insert(0, "/Users/Vincent/Github/exposureTimeCalculator")

from Observation import load_instruments, Observation
import etc_kcwi
from etc_kcwi import get_params
import importlib

# Force reload etc_kcwi to get latest version
importlib.reload(etc_kcwi)

etc_kcwi.datadir = "/Users/Vincent/Github/exposureTimeCalculator/datafiles/kcwi/"

# Constants
SLICER_PARAMS = {
    "Large": {"arcsec_per_slice": 1.35, "pixels_per_spectral_element": 8},
    "Medium": {"arcsec_per_slice": 0.69, "pixels_per_spectral_element": 4},
    "Small": {"arcsec_per_slice": 0.35, "pixels_per_spectral_element": 2},
}

GRATING_PARAMS = {
    "BL": {"A_per_pixel": 0.625, "wave_range": (3500, 5600)},  # Fixed: was 0.5, should be 0.625
    "BM": {"A_per_pixel": 0.28, "wave_range": (3500, 5600)},
    "BH2": {"A_per_pixel": 0.125, "wave_range": (4050, 4860)},
}

PIXEL_SCALE = 0.147
H = 6.62607015e-27
C = 2.99792458e10
A_GEO_KECK = np.pi / 4.0 * (1000.0) ** 2

plt.style.use('default')
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9


def calculate_snr_ratio(slicer, grating, gratwave, seeing, t_exp, flux_erg, sky_erg_input,
                       read_noise=2.7, dark_current=0.0, verbose=False, use_kcwi_sky=True):
    """Calculate SNR ratio (Generic/KCWI) as percentage difference.

    Parameters:
    -----------
    use_kcwi_sky : bool
        If True, use KCWI's sky model (for fair comparison)
        If False, use sky_erg_input (for testing Generic ETC's sky handling)

    Note: KCWI ETC has NO dark current, so varying dark_current only affects Generic ETC.
    """

    slicer_p = SLICER_PARAMS[slicer]
    grating_p = GRATING_PARAMS[grating]

    arcsec_per_slice = slicer_p["arcsec_per_slice"]
    pixels_per_spectral = slicer_p["pixels_per_spectral_element"]
    A_per_pixel = grating_p["A_per_pixel"]

    # KCWI ETC - use their calculation directly
    try:
        from etc_kcwi import calculate_snr as kcwi_calculate_snr, get_params, sky_cts, obj_cts, wpA

        # Call KCWI's SNR calculation
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

        # Get SNR at wavelength
        idx = np.argmin(np.abs(w_kcwi - gratwave))
        snr_keck = snr_kcwi[idx]

        # Get efficiency for Generic ETC matching
        eff_keck = get_params(w_kcwi, grating)
        eff_keck_at_wave = eff_keck[idx]

        # GET KCWI'S SKY VALUE (not the input!)
        c_s_array = sky_cts(w_kcwi, grating, 1.0)  # Per second
        sky_per_s_per_A_per_arcsec2 = c_s_array[idx]

        # Convert to erg for Generic ETC
        E_photon = H * C / (gratwave * 1e-8)
        sky_photons_per_s = sky_per_s_per_A_per_arcsec2 / (eff_keck_at_wave * A_GEO_KECK)
        sky_erg_kcwi = sky_photons_per_s * E_photon

        # Choose which sky to use for Generic ETC
        if use_kcwi_sky:
            sky_erg = sky_erg_kcwi
        else:
            sky_erg = sky_erg_input

        # CALCULATE KCWI COUNTS MANUALLY FOR COMPARISON
        if verbose:
            w, pA = wpA(flux_erg)

            # Spatial bin
            if slicer == 'Small':
                spatial_bin = seeing * seeing
            elif slicer == 'Medium':
                spatial_bin = seeing * max(0.69, seeing)
            elif slicer == 'Large':
                spatial_bin = seeing * 1.38

            # Spectral bin
            spectral_bin = pixels_per_spectral * A_per_pixel

            # Object counts
            c_o_array = obj_cts(w, pA, grating, t_exp)
            c_o = c_o_array[idx] * spatial_bin * spectral_bin

            # Sky counts
            c_s_array_full = sky_cts(w, grating, t_exp)
            c_s = c_s_array_full[idx] * spatial_bin * spectral_bin

            # Read noise
            pixels_per_arcsec = 1.0 / PIXEL_SCALE
            nslices = seeing / arcsec_per_slice
            pixels_spat_bin = pixels_per_arcsec * nslices
            pixels_per_snr_specbin = spectral_bin / A_per_pixel
            c_r = (read_noise**2) * pixels_per_snr_specbin * pixels_spat_bin

            # Manual SNR verification
            snr_manual_kcwi = c_o / np.sqrt(c_s + c_o + c_r)

            print("\n" + "="*60)
            print(f"KCWI CALCULATION ({slicer}-{grating}, gratwave={gratwave})")
            print("="*60)
            print(f"Spatial bin: {spatial_bin:.4f} arcsec²")
            print(f"Spectral bin: {spectral_bin:.4f} Å")
            print(f"Total area: {spatial_bin * spectral_bin:.4f} arcsec²·Å")
            print(f"Signal electrons: {c_o:.2f} e⁻")
            print(f"Sky electrons: {c_s:.2f} e⁻")
            print(f"Read noise²: {c_r:.2f} e²")
            print(f"pixels_spat_bin: {pixels_spat_bin:.2f}")
            print(f"pixels_per_snr_specbin: {pixels_per_snr_specbin:.2f}")
            print(f"SNR (from KCWI func): {snr_keck:.4f}")
            print(f"SNR (manual calc): {snr_manual_kcwi:.4f}")
            print(f"Efficiency: {eff_keck_at_wave:.4f}")
            print(f"Sky (erg): {sky_erg:.3e} erg/cm²/s/Å/arcsec²")

    except Exception as e:
        print(f"KCWI ETC error: {e}")
        import traceback
        traceback.print_exc()
        return None, None, None, None, None

    # Generic ETC
    instruments, _ = load_instruments()

    try:
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
            Sky=sky_erg,  # Either KCWI's sky or input sky
            RN=read_noise,
            Dark_current=dark_current,  # Use the parameter value
            Throughput=eff_keck_at_wave,
            QE=1.0,
            Atmosphere=1.0,
        )

        # Use the SNR calculated by Generic ETC
        snr_generic = obs.SNR[obs.i]

        # Get detailed values for CSV
        signal_total = obs.Signal_el[obs.i] * obs.pixels_total_source
        sky_total = obs.sky[obs.i] * obs.pixels_total_source
        rn2_total = obs.RN_final[obs.i]**2 * obs.number_pixels_used

        if verbose:
            print("\n" + "="*60)
            print("GENERIC ETC CALCULATION")
            print("="*60)
            print(f"Spatial bin (arcsec²): {obs.source_size_arcsec_after_slit:.4f}")
            print(f"Spectral bin (Å): {obs.dispersion * 10 * pixels_per_spectral:.4f}")
            print(f"Pixels total: {obs.pixels_total_source:.2f}")
            print(f"Signal electrons: {signal_total:.2f} e⁻")
            print(f"Sky electrons: {sky_total:.2f} e⁻")
            print(f"Read noise²: {rn2_total:.2f} e²")
            print(f"SNR: {snr_generic:.4f}")
            print(f"Signal per pixel: {obs.Signal_el[obs.i]:.4f} e⁻/pix")
            print(f"Sky per pixel: {obs.sky[obs.i]:.4f} e⁻/pix")
            print(f"factor_CU2el: {obs.factor_CU2el[obs.i]:.3e}")
            print(f"number_pixels_used: {obs.number_pixels_used:.2f}")
            print(f"factor (SNR scaling): {obs.factor[obs.i]:.4f}")

            # Manual SNR calculation
            snr_manual = signal_total / np.sqrt(signal_total + sky_total + rn2_total)
            print(f"\nManual SNR (no factor): {snr_manual:.4f}")
            print(f"Ratio (Manual/Generic): {snr_manual/snr_generic:.4f}")
            print("\n" + "="*60)
            print("COMPARISON")
            print("="*60)
            print(f"Signal ratio (Generic/KCWI): {signal_total/c_o:.4f}" if verbose else "")
            print(f"Sky ratio (Generic/KCWI): {sky_total/c_s:.4f}" if verbose else "")
            print(f"SNR ratio (Generic/KCWI): {snr_generic/snr_keck:.4f}")
            print(f"SNR difference: {100*(snr_generic/snr_keck - 1):.2f}%")

        # Return percentage difference + details
        pct_diff = 100.0 * (snr_generic / snr_keck - 1.0)
        return pct_diff, snr_generic, snr_keck, signal_total, sky_total

    except Exception as e:
        print(f"Generic ETC error: {e}")
        return None, None, None, None, None


print("="*80)
print("PRE-COMPUTING DATA - FULL TEST MODE")
print("="*80)

# FULL TEST: All 6 configs
configs = [
    ("Small", "BL", 4550),
    ("Small", "BM", 4550),
    ("Medium", "BL", 4550),
    ("Medium", "BM", 4550),
    # ("Large", "BL", 4550),
    # ("Large", "BM", 4550),
]

# Nominal parameters
seeing_nominal = 0.75
t_exp_nominal = 3600.0
flux_nominal = 5.6e-19
sky_nominal = 8e-18  # IGNORED - KCWI uses its own sky model
rn_nominal = 2.7
dark_nominal = 0.0  # ALWAYS 0 (KCWI doesn't have dark)

# Pre-compute all datasets
data = {}
csv_rows = []  # For CSV export

# TEST MODE: 3 points per parameter for quick validation
# For full publication-quality plots, increase to 12 points:
# t_exp_values = np.logspace(2, 4, 12)
# flux_factors = np.logspace(-0.5, 0.5, 12)
# etc.

# 1. Exposure time
print("  Computing exposure time variations...")
t_exp_values = np.array([1000, 3600, 10000])  # 3 points for quick test
data['t_exp'] = {'x': t_exp_values, 'y': {}}
for slicer, grating, gratwave in configs:
    key = f"{slicer[0]}-{grating}"
    results = []
    for t_exp in t_exp_values:
        pct_diff, snr_gen, snr_keck, sig, sky = calculate_snr_ratio(
            slicer, grating, gratwave, seeing_nominal, t_exp,
            flux_nominal, sky_nominal, rn_nominal, dark_nominal)
        results.append(pct_diff)
        if pct_diff is not None:
            csv_rows.append({
                'config': key, 'param': 't_exp', 'value': t_exp,
                'pct_diff': pct_diff, 'snr_gen': snr_gen, 'snr_keck': snr_keck,
                'signal_total': sig, 'sky_total': sky
            })
    data['t_exp']['y'][key] = results

# 2. Source flux
print("  Computing source flux variations...")
flux_factors = np.array([0.5, 1.0, 2.0])  # QUICK: 3 points | FULL: np.logspace(-0.5, 0.5, 12)
data['flux'] = {'x': flux_factors * flux_nominal * 1e19, 'y': {}}
for slicer, grating, gratwave in configs:
    key = f"{slicer[0]}-{grating}"
    results = []
    for factor in flux_factors:
        pct_diff, snr_gen, snr_keck, sig, sky = calculate_snr_ratio(
            slicer, grating, gratwave, seeing_nominal, t_exp_nominal,
            flux_nominal * factor, sky_nominal, rn_nominal, dark_nominal)
        results.append(pct_diff)
        if pct_diff is not None:
            csv_rows.append({
                'config': key, 'param': 'flux', 'value': flux_nominal * factor,
                'pct_diff': pct_diff, 'snr_gen': snr_gen, 'snr_keck': snr_keck,
                'signal_total': sig, 'sky_total': sky
            })
    data['flux']['y'][key] = results

# 3. Sky background - Test Generic ETC's sky handling (KCWI always uses own sky)
print("  Computing sky background variations...")
sky_factors = np.array([0.5, 1.0, 2.0])  # QUICK: 3 points | FULL: np.logspace(-0.5, 0.5, 12)
data['sky'] = {'x': sky_factors * sky_nominal * 1e18, 'y': {}}
for slicer, grating, gratwave in configs:
    key = f"{slicer[0]}-{grating}"
    results = []
    for factor in sky_factors:
        pct_diff, snr_gen, snr_keck, sig, sky = calculate_snr_ratio(
            slicer, grating, gratwave, seeing_nominal, t_exp_nominal,
            flux_nominal, sky_nominal * factor, rn_nominal, dark_nominal,
            use_kcwi_sky=False)  # Use input sky to test Generic ETC
        results.append(pct_diff)
        if pct_diff is not None:
            csv_rows.append({
                'config': key, 'param': 'sky', 'value': sky_nominal * factor,
                'pct_diff': pct_diff, 'snr_gen': snr_gen, 'snr_keck': snr_keck,
                'signal_total': sig, 'sky_total': sky
            })
    data['sky']['y'][key] = results

# 4. Read noise
print("  Computing read noise variations...")
rn_values = np.array([1.0, 2.7, 5.0])  # QUICK: 3 points | FULL: np.linspace(0.5, 10.0, 12)
data['rn'] = {'x': rn_values, 'y': {}}
for slicer, grating, gratwave in configs:
    key = f"{slicer[0]}-{grating}"
    results = []
    for rn in rn_values:
        pct_diff, snr_gen, snr_keck, sig, sky = calculate_snr_ratio(
            slicer, grating, gratwave, seeing_nominal, t_exp_nominal,
            flux_nominal, sky_nominal, rn, dark_nominal)
        results.append(pct_diff)
        if pct_diff is not None:
            csv_rows.append({
                'config': key, 'param': 'rn', 'value': rn,
                'pct_diff': pct_diff, 'snr_gen': snr_gen, 'snr_keck': snr_keck,
                'signal_total': sig, 'sky_total': sky
            })
    data['rn']['y'][key] = results

# 5. Dark current - Test Generic ETC's dark handling (KCWI has NO dark, always 0)
print("  Computing dark current variations...")
dark_values = np.array([0.0, 5.0, 10.0])  # QUICK: 3 points | FULL: np.linspace(0.1, 20.0, 12)
data['dark'] = {'x': dark_values, 'y': {}}
for slicer, grating, gratwave in configs:
    key = f"{slicer[0]}-{grating}"
    results = []
    for dark in dark_values:
        pct_diff, snr_gen, snr_keck, sig, sky = calculate_snr_ratio(
            slicer, grating, gratwave, seeing_nominal, t_exp_nominal,
            flux_nominal, sky_nominal, rn_nominal, dark)  # Use variable dark value!
        results.append(pct_diff)
        if pct_diff is not None:
            csv_rows.append({
                'config': key, 'param': 'dark', 'value': dark,
                'pct_diff': pct_diff, 'snr_gen': snr_gen, 'snr_keck': snr_keck,
                'signal_total': sig, 'sky_total': sky
            })
    data['dark']['y'][key] = results

# 6. Seeing - Both ETCs vary
print("  Computing seeing variations...")
seeing_values = np.array([0.5, 0.75, 1.2])  # QUICK: 3 points | FULL: np.linspace(0.4, 2.0, 12)
data['seeing'] = {'x': seeing_values, 'y': {}}
for slicer, grating, gratwave in configs:
    key = f"{slicer[0]}-{grating}"
    results = []
    for seeing in seeing_values:
        pct_diff, snr_gen, snr_keck, sig, sky = calculate_snr_ratio(
            slicer, grating, gratwave, seeing, t_exp_nominal,
            flux_nominal, sky_nominal, rn_nominal, dark_nominal)
        results.append(pct_diff)
        if pct_diff is not None:
            csv_rows.append({
                'config': key, 'param': 'seeing', 'value': seeing,
                'pct_diff': pct_diff, 'snr_gen': snr_gen, 'snr_keck': snr_keck,
                'signal_total': sig, 'sky_total': sky
            })
    data['seeing']['y'][key] = results

# Save CSV
import pandas as pd
df = pd.DataFrame(csv_rows)
csv_path = script_dir + '/kcwi_snr_comparison.csv'
df.to_csv(csv_path, index=False)

print("\n" + "="*80)
print("PLOTTING (fast!)...")
print("="*80)

# Find global y-limits
all_values = []
for dataset in data.values():
    for values in dataset['y'].values():
        all_values.extend([v for v in values if v is not None])

if len(all_values) == 0:
    print("\n⚠️  WARNING: No valid data computed! All calculations returned None.")
    print("Check the error messages above for details.")
    ylim = [-10, 10]  # Default range for empty plots
else:
    ymin, ymax = np.min(all_values), np.max(all_values)
    ymargin = (ymax - ymin) * 0.1
    ylim = [ymin - ymargin, ymax + ymargin]

# Colors
colors = plt.cm.tab10(np.linspace(0, 1, len(configs)))



#%%
# Create figure
fig = plt.figure(figsize=(18, 10))
gs = GridSpec(2, 3, figure=fig, hspace=0.3, wspace=0.3)

# ROW 1: Validation plots with shared y-axis (excellent agreement, <1%)
# Plot 1: Exposure time
ax1 = fig.add_subplot(gs[0, 0])
for i, key in enumerate(data['t_exp']['y'].keys()):
    ax1.semilogx(data['t_exp']['x'], data['t_exp']['y'][key],
                color=colors[i], label=key, linewidth=2.5, alpha=0.9, marker='o', markersize=6)
ax1.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax1.axvline(t_exp_nominal, color='red', linestyle=':', linewidth=2, alpha=0.5, label='Nominal')
ax1.set_xlabel('Exposure Time (s)')
ax1.set_ylabel('SNR Difference (%)')
ax1.set_title('SNR Ratio vs Exposure Time')
ax1.legend(loc='best', ncol=2, framealpha=0.9, fontsize=8)
ax1.grid(True, alpha=0.3, which='both')
# ax1.set_ylim(ylim)

# Plot 2: Source brightness
ax2 = fig.add_subplot(gs[0, 1], sharey=ax1)
for i, key in enumerate(data['flux']['y'].keys()):
    ax2.semilogx(data['flux']['x'], data['flux']['y'][key],
                color=colors[i], label=key, linewidth=2.5, alpha=0.9, marker='o', markersize=6)
ax2.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax2.axvline(flux_nominal * 1e19, color='red', linestyle=':', linewidth=2, alpha=0.5)
ax2.set_xlabel('Source Flux (×10⁻¹⁹ erg/cm²/s/arcsec²/Å)')
ax2.set_title('SNR Ratio vs Source Brightness')
ax2.legend(loc='best', ncol=2, framealpha=0.9, fontsize=8)
ax2.grid(True, alpha=0.3, which='both')
plt.setp(ax2.get_yticklabels(), visible=False)

# Plot 3: Read noise
ax3 = fig.add_subplot(gs[0, 2], sharey=ax1)
for i, key in enumerate(data['rn']['y'].keys()):
    ax3.plot(data['rn']['x'], data['rn']['y'][key],
            color=colors[i], label=key, linewidth=2.5, alpha=0.9, marker='o', markersize=6)
ax3.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax3.axvline(rn_nominal, color='red', linestyle=':', linewidth=2, alpha=0.5)
ax3.set_xlabel('Read Noise (e⁻/pix)')
ax3.set_title('SNR Ratio vs Read Noise')
ax3.legend(loc='best', ncol=2, framealpha=0.9, fontsize=8)
ax3.grid(True, alpha=0.3)
plt.setp(ax3.get_yticklabels(), visible=False)

# ROW 2: Sky and Dark tests (KCWI fixed, Generic varies) with separate y-axis
# Plot 4: Sky background
ax4 = fig.add_subplot(gs[1, 0])
for i, key in enumerate(data['sky']['y'].keys()):
    ax4.semilogx(data['sky']['x'], data['sky']['y'][key],
                color=colors[i], label=key, linewidth=2.5, alpha=0.9, marker='o', markersize=6)
ax4.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax4.axvline(sky_nominal * 1e18, color='red', linestyle=':', linewidth=2, alpha=0.5)
ax4.set_xlabel('Sky Background (×10⁻¹⁸ erg/cm²/s/arcsec²/Å)')
ax4.set_ylabel('SNR Difference (%)')
ax4.set_title('SNR vs Sky (Generic varies, KCWI fixed)')
ax4.legend(loc='best', ncol=2, framealpha=0.9, fontsize=8)
ax4.grid(True, alpha=0.3, which='both')

# Plot 5: Dark current
ax5 = fig.add_subplot(gs[1, 1], sharey=ax4)
for i, key in enumerate(data['dark']['y'].keys()):
    ax5.plot(data['dark']['x'], data['dark']['y'][key],
            color=colors[i], label=key, linewidth=2.5, alpha=0.9, marker='o', markersize=6)
ax5.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax5.axvline(0.0, color='red', linestyle=':', linewidth=2, alpha=0.5)
ax5.set_xlabel('Dark Current (e⁻/pix/hour)')
ax5.set_title('SNR vs Dark (Generic varies, KCWI=0)')
ax5.legend(loc='best', ncol=2, framealpha=0.9, fontsize=8)
ax5.grid(True, alpha=0.3)
plt.setp(ax5.get_yticklabels(), visible=False)

# Plot 6: Seeing - Both ETCs vary
ax6 = fig.add_subplot(gs[1, 2], sharey=ax4)
for i, key in enumerate(data['seeing']['y'].keys()):
    ax6.plot(data['seeing']['x'], data['seeing']['y'][key],
            color=colors[i], label=key, linewidth=2.5, alpha=0.9, marker='o', markersize=6)
ax6.axhline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.7)
ax6.axvline(seeing_nominal, color='red', linestyle=':', linewidth=2, alpha=0.5)
ax6.set_xlabel('Seeing (arcsec)')
ax6.set_title('SNR vs Seeing (Both vary)')
ax6.legend(loc='best', ncol=2, framealpha=0.9, fontsize=8)
ax6.grid(True, alpha=0.3)
plt.setp(ax6.get_yticklabels(), visible=False)

fig.suptitle('Generic ETC vs KCWI ETC: SNR Calculation Validation',
            fontsize=16,  y=0.95)

output_path = '/Users/Vincent/Github/generic-etc/test/cross-validation/kcwi_snr_ratios_comprehensive.png'
plt.tight_layout()
plt.savefig(output_path, dpi=200, bbox_inches='tight')
plt.show()

# %%

# %%
