#!/usr/bin/env python3
"""
Test des conversions Generic ETC : flux → électrons

Vérifie que toutes les conversions physiques sont correctes
"""

import sys
import os

# Store original directory
_original_dir = os.path.dirname(os.path.abspath(__file__))

# If running from project root, need to be in notebooks for data access
if os.path.exists(os.path.join(_original_dir, 'notebooks')):
    os.chdir(os.path.join(_original_dir, 'notebooks'))
    sys.path.insert(0, '.')
else:
    sys.path.insert(0, 'notebooks')

import numpy as np
from Observation import Observation, load_instruments

def test_conversions_galex_fuv():
    """
    Test détaillé des conversions pour GALEX FUV
    """
    print("="*80)
    print("TEST DES CONVERSIONS - GALEX FUV")
    print("="*80)

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

    print("\nParamètres GALEX FUV:")
    for key in ['wavelength', 'Collecting_area', 'pixel_scale', 'dispersion',
                'Throughput', 'QE', 'Atmosphere', 'Signal', 'Sky', 'Dark_current',
                'Size_source', 'Slitwidth', 'Slitlength', 'Line_width', 'Bandwidth', 'RN']:
        if key in params:
            print(f"  {key:20s} = {params[key]}")

    # Test parameters
    t_exp = 1000.0
    readout_time = params.get('Readout', 0.0)
    acq_time = (t_exp + readout_time) / 3600.0

    # Create observation
    obs = Observation(
        instruments=instruments,
        instrument=instrument_name,
        exposure_time=t_exp,
        acquisition_time=acq_time,
        SNR_res="per pix",
        IFS=False,
        test=True
    )

    print("\n" + "="*80)
    print("1. ÉNERGIE DU PHOTON")
    print("="*80)

    # Constants
    h = 6.62606896e-27  # erg s
    c_light = 299792.458  # km/s
    c = c_light * 1e5  # cm/s
    wavelength_nm = params['wavelength']
    wavelength_cm = wavelength_nm / 1e7  # nm → cm

    photon_energy_manual = h * c / wavelength_cm  # erg

    print(f"\nCalcul manuel:")
    print(f"  h = {h:.6e} erg·s")
    print(f"  c = {c:.6e} cm/s")
    print(f"  λ = {wavelength_nm} nm = {wavelength_cm:.6e} cm")
    print(f"  E_photon = h × c / λ = {photon_energy_manual:.6e} erg")

    # Extract from Generic ETC (need to access internal calculation)
    # We'll recalculate to compare
    wavelength_cm_etc = obs.wavelength / 1e7
    photon_energy_etc = h * c / wavelength_cm_etc

    print(f"\nGeneric ETC:")
    print(f"  E_photon = {photon_energy_etc:.6e} erg")

    ratio = photon_energy_manual / photon_energy_etc
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    if abs(ratio - 1.0) < 1e-6:
        print("✓ Énergie photon correcte")
    else:
        print("✗ ERREUR: Énergie photon différente!")

    print("\n" + "="*80)
    print("2. EFFECTIVE AREA (conversion flux → électrons)")
    print("="*80)

    # Manual calculation
    collecting_area_cm2 = params['Collecting_area'] * 1e4  # m² → cm²
    effective_area_manual = (collecting_area_cm2 * params['Throughput'] *
                            params['QE'] * params['Atmosphere'] / photon_energy_manual)

    print(f"\nCalcul manuel:")
    print(f"  Collecting_area = {params['Collecting_area']} m² = {collecting_area_cm2} cm²")
    print(f"  Throughput = {params['Throughput']}")
    print(f"  QE = {params['QE']}")
    print(f"  Atmosphere = {params['Atmosphere']}")
    print(f"  effective_area = A × T × QE × Atm / E_photon")
    print(f"                 = {effective_area_manual:.6e} cm²·photon/erg")

    # Generic ETC should have this calculated
    # Access obs attributes if possible
    print(f"\nGeneric ETC:")
    print(f"  effective_area = {obs.effective_area[obs.i]:.6e} cm²·photon/erg")

    ratio = effective_area_manual / obs.effective_area[obs.i]
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    if abs(ratio - 1.0) < 1e-6:
        print("✓ Effective area correcte")
    else:
        print("✗ ERREUR: Effective area différente!")

    print("\n" + "="*80)
    print("3. AIRES SPATIALES (source vs slit)")
    print("="*80)

    fwhm_sigma_ratio = 2.35

    # Source size after slit
    source_size_spatial = min(params['Size_source'], params['Slitwidth'])
    source_size_spectral = min(params['Size_source'] * fwhm_sigma_ratio, params['Slitlength'])
    source_size_arcsec_manual = source_size_spatial * source_size_spectral

    # Slit size
    slit_size_arcsec_manual = params['Slitwidth'] * params['Slitlength']

    print(f"\nCalcul manuel:")
    print(f"  Size_source = {params['Size_source']} arcsec")
    print(f"  Slitwidth = {params['Slitwidth']} arcsec")
    print(f"  Slitlength = {params['Slitlength']} arcsec")
    print(f"  fwhm_sigma_ratio = {fwhm_sigma_ratio}")
    print(f"\n  Source spatial = min({params['Size_source']}, {params['Slitwidth']}) = {source_size_spatial} arcsec")
    print(f"  Source spectral = min({params['Size_source'] * fwhm_sigma_ratio}, {params['Slitlength']}) = {source_size_spectral} arcsec")
    print(f"  source_size_arcsec = {source_size_arcsec_manual} arcsec²")
    print(f"\n  slit_size_arcsec = {params['Slitwidth']} × {params['Slitlength']} = {slit_size_arcsec_manual} arcsec²")

    print(f"\nGeneric ETC:")
    # Need to access internal variables - recalculate
    source_size_arcsec_etc = (np.minimum(params['Size_source'] * fwhm_sigma_ratio, params['Slitlength']) *
                              np.minimum(params['Size_source'], params['Slitwidth']))
    slit_size_arcsec_etc = params['Slitwidth'] * params['Slitlength']
    print(f"  source_size_arcsec = {source_size_arcsec_etc} arcsec²")
    print(f"  slit_size_arcsec = {slit_size_arcsec_etc} arcsec²")

    print(f"\nRatio slit/source: {slit_size_arcsec_manual / source_size_arcsec_manual:.2f}")
    print(f"→ Le sky devrait être ~{slit_size_arcsec_manual / source_size_arcsec_manual:.1f}× plus grand que le signal")

    print("\n" + "="*80)
    print("4. PIXELS_TOTAL_SOURCE")
    print("="*80)

    sigma_x_pix = params['Size_source'] / params['pixel_scale'] / fwhm_sigma_ratio
    sigma_y_pix = params['Line_width'] / params['dispersion'] / fwhm_sigma_ratio
    pixels_total_source_manual = 2 * np.pi * sigma_x_pix * sigma_y_pix

    print(f"\nCalcul manuel:")
    print(f"  sigma_x = Size_source / pixel_scale / 2.35")
    print(f"          = {params['Size_source']} / {params['pixel_scale']} / {fwhm_sigma_ratio}")
    print(f"          = {sigma_x_pix:.3f} pix")
    print(f"  sigma_y = Line_width / dispersion / 2.35")
    print(f"          = {params['Line_width']} / {params['dispersion']} / {fwhm_sigma_ratio}")
    print(f"          = {sigma_y_pix:.3f} pix")
    print(f"  pixels_total_source = 2π × sigma_x × sigma_y")
    print(f"                      = {pixels_total_source_manual:.3f} pix")

    print(f"\nGeneric ETC:")
    print(f"  pixels_total_source = {obs.pixels_total_source:.3f} pix")

    ratio = pixels_total_source_manual / obs.pixels_total_source
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    if abs(ratio - 1.0) < 1e-6:
        print("✓ pixels_total_source correct")
    else:
        print("✗ ERREUR: pixels_total_source différent!")

    print("\n" + "="*80)
    print("5. FACTOR_CU2EL (signal)")
    print("="*80)

    arcsec2str = ((np.pi / 180.0) / 3600.0) ** 2  # arcsec² → sr
    spectral_width_signal = min(params['Line_width'], params['Bandwidth'])

    factor_CU2el_manual = (effective_area_manual * arcsec2str *
                          spectral_width_signal *
                          source_size_arcsec_manual / pixels_total_source_manual)

    print(f"\nCalcul manuel:")
    print(f"  arcsec2str = {arcsec2str:.6e} sr/arcsec²")
    print(f"  spectral_width = min(Line_width, Bandwidth)")
    print(f"                 = min({params['Line_width']}, {params['Bandwidth']})")
    print(f"                 = {spectral_width_signal} Å")
    print(f"  factor_CU2el = effective_area × arcsec2str × spectral_width")
    print(f"                 × source_size_arcsec / pixels_total_source")
    print(f"               = {factor_CU2el_manual:.6e}")

    print(f"\nGeneric ETC:")
    print(f"  factor_CU2el = {obs.factor_CU2el[obs.i]:.6e}")

    ratio = factor_CU2el_manual / obs.factor_CU2el[obs.i]
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    if abs(ratio - 1.0) < 1e-3:
        print("✓ factor_CU2el correct (within 0.1%)")
    else:
        print(f"⚠ ATTENTION: factor_CU2el différent de {(ratio-1)*100:.2f}%")

    print("\n" + "="*80)
    print("6. FACTOR_CU2EL_SKY")
    print("="*80)

    spectral_width_sky = max(min(params['Line_width'], params['Bandwidth']), params['dispersion'])

    factor_CU2el_sky_manual = (effective_area_manual * arcsec2str *
                               spectral_width_sky *
                               slit_size_arcsec_manual / pixels_total_source_manual)

    print(f"\nCalcul manuel:")
    print(f"  spectral_width_sky = max(min(Line_width, Bandwidth), dispersion)")
    print(f"                     = max(min({params['Line_width']}, {params['Bandwidth']}), {params['dispersion']})")
    print(f"                     = {spectral_width_sky} Å")
    print(f"  factor_CU2el_sky = effective_area × arcsec2str × spectral_width_sky")
    print(f"                     × slit_size_arcsec / pixels_total_source")
    print(f"                   = {factor_CU2el_sky_manual:.6e}")

    print(f"\nGeneric ETC:")
    print(f"  factor_CU2el_sky = {obs.factor_CU2el_sky[obs.i]:.6e}")

    ratio = factor_CU2el_sky_manual / obs.factor_CU2el_sky[obs.i]
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    print(f"\nRatio factor_sky/factor_signal: {obs.factor_CU2el_sky[obs.i] / obs.factor_CU2el[obs.i]:.2f}")
    print(f"  (attendu: ~{slit_size_arcsec_manual / source_size_arcsec_manual:.2f} si même largeur spectrale)")

    if abs(ratio - 1.0) < 1e-3:
        print("✓ factor_CU2el_sky correct (within 0.1%)")
    else:
        print(f"⚠ ATTENTION: factor_CU2el_sky différent de {(ratio-1)*100:.2f}%")

    print("\n" + "="*80)
    print("7. SIGNAL_EL (électrons)")
    print("="*80)

    Signal_el_manual = params['Signal'] * factor_CU2el_manual * t_exp

    print(f"\nCalcul manuel:")
    print(f"  Signal flux = {params['Signal']:.6e} erg/cm²/s/arcsec²/Å")
    print(f"  factor_CU2el = {factor_CU2el_manual:.6e}")
    print(f"  exposure_time = {t_exp} s")
    print(f"  Signal_el = Signal × factor_CU2el × t_exp")
    print(f"            = {Signal_el_manual:.6e} e⁻")

    print(f"\nGeneric ETC:")
    print(f"  Signal_el = {obs.Signal_el[obs.i]:.6e} e⁻")

    ratio = Signal_el_manual / obs.Signal_el[obs.i]
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    if abs(ratio - 1.0) < 1e-3:
        print("✓ Signal_el correct (within 0.1%)")
    else:
        print(f"⚠ ATTENTION: Signal_el différent de {(ratio-1)*100:.2f}%")

    print("\n" + "="*80)
    print("8. SKY (électrons)")
    print("="*80)

    sky_manual = params['Sky'] * factor_CU2el_sky_manual * t_exp

    print(f"\nCalcul manuel:")
    print(f"  Sky flux = {params['Sky']:.6e} erg/cm²/s/arcsec²/Å")
    print(f"  factor_CU2el_sky = {factor_CU2el_sky_manual:.6e}")
    print(f"  exposure_time = {t_exp} s")
    print(f"  sky = Sky × factor_CU2el_sky × t_exp")
    print(f"      = {sky_manual:.6e} e⁻")

    print(f"\nGeneric ETC:")
    print(f"  sky = {obs.sky[obs.i]:.6e} e⁻")

    ratio = sky_manual / obs.sky[obs.i]
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    print(f"\nRatio sky/Signal_el: {obs.sky[obs.i] / obs.Signal_el[obs.i]:.2f}")
    print(f"  (flux_sky/flux_signal = {params['Sky']/params['Signal']:.2f})")
    print(f"  (aire_slit/aire_source = {slit_size_arcsec_manual/source_size_arcsec_manual:.2f})")
    print(f"  → attendu: ~{(params['Sky']/params['Signal']) * (slit_size_arcsec_manual/source_size_arcsec_manual):.2f}")

    if abs(ratio - 1.0) < 1e-3:
        print("✓ sky correct (within 0.1%)")
    else:
        print(f"⚠ ATTENTION: sky différent de {(ratio-1)*100:.2f}%")

    print("\n" + "="*80)
    print("9. DARK_CURRENT_F (électrons)")
    print("="*80)

    Dark_current_f_manual = params['Dark_current'] * t_exp / 3600.0 * obs.number_pixels_used

    print(f"\nCalcul manuel:")
    print(f"  Dark_current = {params['Dark_current']} e⁻/pix/hour")
    print(f"  exposure_time = {t_exp} s = {t_exp/3600:.4f} hour")
    print(f"  number_pixels_used = {obs.number_pixels_used}")
    print(f"  Dark_current_f = Dark_current × t_exp / 3600 × n_pix")
    print(f"                 = {Dark_current_f_manual:.6e} e⁻")

    print(f"\nGeneric ETC:")
    print(f"  Dark_current_f = {obs.Dark_current_f[obs.i]:.6e} e⁻")

    ratio = Dark_current_f_manual / obs.Dark_current_f[obs.i]
    print(f"\nRatio (manuel/ETC): {ratio:.10f}")

    if abs(ratio - 1.0) < 1e-6:
        print("✓ Dark_current_f correct")
    else:
        print("✗ ERREUR: Dark_current_f différent!")

    print("\n" + "="*80)
    print("10. READ NOISE (électrons)")
    print("="*80)

    print(f"\nParamètre:")
    print(f"  RN = {params['RN']} e⁻")

    print(f"\nGeneric ETC:")
    print(f"  RN_final = {obs.RN_final[obs.i]:.6e} e⁻")
    print(f"  EMCCD = {obs.EMCCD}")

    if not obs.EMCCD:
        if abs(obs.RN_final[obs.i] - params['RN']) < 1e-6:
            print("✓ RN_final = RN (mode normal)")
        else:
            print(f"✗ ERREUR: RN_final ≠ RN!")
    else:
        print(f"  EM_gain = {obs.EM_gain}")
        print(f"  RN_EMCCD = {obs.RN_EMCCD}")
        rn_expected = obs.RN_EMCCD / obs.EM_gain
        if abs(obs.RN_final[obs.i] - rn_expected) < 1e-6:
            print(f"✓ RN_final = RN_EMCCD / EM_gain (mode EMCCD)")
        else:
            print(f"✗ ERREUR: RN_final incorrect en mode EMCCD!")

    print("\n" + "="*80)
    print("RÉSUMÉ")
    print("="*80)

    print(f"\nValeurs finales:")
    print(f"  Signal_el = {obs.Signal_el[obs.i]:.6e} e⁻")
    print(f"  sky = {obs.sky[obs.i]:.6e} e⁻")
    print(f"  Dark_current_f = {obs.Dark_current_f[obs.i]:.6e} e⁻")
    print(f"  RN_final = {obs.RN_final[obs.i]:.6e} e⁻")
    print(f"  SNR = {obs.SNR[obs.i]:.6e}")

    print(f"\nRapports:")
    print(f"  sky / Signal_el = {obs.sky[obs.i] / obs.Signal_el[obs.i]:.2f}")
    print(f"  Dark / Signal_el = {obs.Dark_current_f[obs.i] / obs.Signal_el[obs.i]:.4f}")
    print(f"  RN / Signal_el = {obs.RN_final[obs.i] / obs.Signal_el[obs.i]:.4f}")

    return {
        'instrument': instrument_name,
        'Signal_el': obs.Signal_el[obs.i],
        'sky': obs.sky[obs.i],
        'Dark_current_f': obs.Dark_current_f[obs.i],
        'RN_final': obs.RN_final[obs.i],
        'SNR': obs.SNR[obs.i],
        'pixels_total_source': obs.pixels_total_source
    }

if __name__ == "__main__":
    result = test_conversions_galex_fuv()
