#!/usr/bin/env python3
"""
Test de l'hypothèse : le problème vient de la différence d'aire pour le sky
Si on renormalise le Sky d'entrée, tout devrait matcher !
"""

import sys
sys.path.insert(0, 'notebooks')

def compute_manual_calculation():
    """
    Calcul manuel étape par étape pour comprendre exactement ce qui se passe
    """
    import numpy as np

    print("="*80)
    print("CALCUL MANUEL - GALEX FUV")
    print("="*80)

    # Paramètres GALEX FUV
    wavelength = 152.8  # nm
    Collecting_area = 0.196  # m²
    pixel_scale = 1.5  # arcsec/pix
    dispersion = 2.5  # Å/pix
    Throughput = 0.08
    QE = 0.12
    Atmosphere = 1.0
    Signal = 5.6e-19  # erg/cm²/s/arcsec²/Å
    Sky = 1.53e-19  # erg/cm²/s/arcsec²/Å
    Dark_current = 0.049  # e⁻/pix/hour
    Size_source = 4.0  # arcsec
    Slitwidth = 4.0  # arcsec
    Slitlength = 20.0  # arcsec
    Line_width = 2.0  # Å
    Bandwidth = 250.0  # Å
    IFS = False
    fwhm_sigma_ratio = 2.35

    t_exp = 1000.0  # s

    print(f"\nParamètres:")
    print(f"  Size_source: {Size_source} arcsec")
    print(f"  Slitwidth: {Slitwidth} arcsec")
    print(f"  Slitlength: {Slitlength} arcsec")
    print(f"  pixel_scale: {pixel_scale} arcsec/pix")
    print(f"  dispersion: {dispersion} Å/pix")
    print(f"  fwhm_sigma_ratio: {fwhm_sigma_ratio}")

    # Calcul des aires (ligne 2472-2473 de Observation.py)
    source_size_arcsec_after_slit = (np.minimum(Size_source*fwhm_sigma_ratio, Slitlength) *
                                     np.minimum(Size_source*fwhm_sigma_ratio, Slitwidth))

    slit_size_arcsec_after_slit = (np.minimum(Size_source*fwhm_sigma_ratio, Slitlength) *
                                   Slitwidth)

    print(f"\nAires calculées (lignes 2472-2473):")
    print(f"  source_size_arcsec_after_slit: {source_size_arcsec_after_slit:.2f} arcsec²")
    print(f"  slit_size_arcsec_after_slit: {slit_size_arcsec_after_slit:.2f} arcsec²")
    print(f"  pixel_scale²: {pixel_scale**2:.2f} arcsec²")

    # Calcul pixels_total_source (ligne 2399)
    # Simplifié pour IFS=False
    PSF_lambda_pix = 10*wavelength / 1300 / dispersion  # Spectral_resolution ~ 1300
    source_spectral_pixels = np.maximum(1, np.sqrt(PSF_lambda_pix**2 +
                                                   (np.minimum(Line_width, Bandwidth) / dispersion)**2))
    source_spatial_pixels = np.maximum(1, np.minimum(np.sqrt(Size_source**2) * fwhm_sigma_ratio / pixel_scale,
                                                     Slitlength / pixel_scale))
    source_size = (np.maximum(np.minimum(Size_source * fwhm_sigma_ratio, Slitlength) / pixel_scale, 1) *
                  np.sqrt(source_spatial_pixels**2 + source_spectral_pixels**2))
    pixels_total_source = source_size  # IFS=False, donc facteur 1

    print(f"\nPixels (lignes 2395-2399):")
    print(f"  PSF_lambda_pix: {PSF_lambda_pix:.2f}")
    print(f"  source_spectral_pixels: {source_spectral_pixels:.2f}")
    print(f"  source_spatial_pixels: {source_spatial_pixels:.2f}")
    print(f"  source_size: {source_size:.2f}")
    print(f"  pixels_total_source: {pixels_total_source:.2f}")

    # Calcul effective_area
    collecting_area_cm2 = Collecting_area * 1e4
    total_throughput = Throughput * QE * Atmosphere
    effective_area = collecting_area_cm2 * total_throughput
    arcsec2str = (np.pi/180/3600)**2

    print(f"\nEffective area:")
    print(f"  collecting_area: {collecting_area_cm2:.1f} cm²")
    print(f"  total_throughput: {total_throughput:.4f}")
    print(f"  effective_area: {effective_area:.2f} cm²")
    print(f"  arcsec2str: {arcsec2str:.4e}")

    # Formula test=True (ligne 2481, 2483)
    factor_CU2el_tot = (effective_area * arcsec2str *
                        np.minimum(Line_width, Bandwidth) *
                        source_size_arcsec_after_slit / pixels_total_source)

    factor_CU2el_sky_tot = (effective_area * arcsec2str *
                            np.maximum(np.minimum(Line_width, Bandwidth), dispersion) *
                            slit_size_arcsec_after_slit / pixels_total_source)

    print(f"\nFacteurs (test=True, lignes 2481-2483):")
    print(f"  factor_CU2el_tot: {factor_CU2el_tot:.6e}")
    print(f"  factor_CU2el_sky_tot: {factor_CU2el_sky_tot:.6e}")
    print(f"  Ratio sky/signal: {factor_CU2el_sky_tot/factor_CU2el_tot:.2f}")

    # Calcul Signal et Sky en électrons (lignes 2539, 2545)
    # Conversion ergs -> LU
    E_photon = 1.986e-8 / wavelength
    Signal_LU = Signal / E_photon * (np.pi/180/3600)**2
    Sky_CU = Sky / E_photon * (np.pi/180/3600)**2

    flux_fraction_slit_applied = 1.0  # Simplifié

    Signal_el = Signal_LU * factor_CU2el_tot * t_exp * flux_fraction_slit_applied
    sky = Sky_CU * factor_CU2el_sky_tot * t_exp
    Dark_current_f = Dark_current * t_exp / 3600

    print(f"\nÉlectrons (lignes 2539, 2545):")
    print(f"  Signal_LU: {Signal_LU:.6e}")
    print(f"  Sky_CU: {Sky_CU:.6e}")
    print(f"  Signal_el: {Signal_el:.6e} e⁻")
    print(f"  sky: {sky:.6e} e⁻")
    print(f"  Dark_current_f: {Dark_current_f:.6e} e⁻")

    # Calcul Astropy (per pixel)
    signal_photons = Signal / E_photon
    sky_photons = Sky / E_photon

    source_eps_astropy = (signal_photons * collecting_area_cm2 * pixel_scale**2 *
                         dispersion * total_throughput)
    sky_eps_astropy = (sky_photons * collecting_area_cm2 * pixel_scale**2 *
                      dispersion * total_throughput)
    dark_eps_astropy = Dark_current / 3600.0

    signal_total_astropy = source_eps_astropy * t_exp
    sky_total_astropy = sky_eps_astropy * t_exp
    dark_total_astropy = dark_eps_astropy * t_exp

    print(f"\n{'='*80}")
    print(f"ASTROPY (per pixel)")
    print(f"{'='*80}")
    print(f"  source_eps: {source_eps_astropy:.6e} e⁻/s/pix")
    print(f"  sky_eps: {sky_eps_astropy:.6e} e⁻/s/pix")
    print(f"  signal_total: {signal_total_astropy:.6e} e⁻")
    print(f"  sky_total: {sky_total_astropy:.6e} e⁻")
    print(f"  dark_total: {dark_total_astropy:.6e} e⁻")

    # Comparaison brute
    print(f"\n{'='*80}")
    print(f"COMPARAISON BRUTE")
    print(f"{'='*80}")
    print(f"  Signal: Generic={Signal_el:.6e}, Astropy={signal_total_astropy:.6e}")
    print(f"  Ratio: {Signal_el/signal_total_astropy:.4f}")
    print(f"  Sky: Generic={sky:.6e}, Astropy={sky_total_astropy:.6e}")
    print(f"  Ratio: {sky/sky_total_astropy:.4f}")
    print(f"  Dark: Generic={Dark_current_f:.6e}, Astropy={dark_total_astropy:.6e}")
    print(f"  Ratio: {Dark_current_f/dark_total_astropy:.4f}")

    # Correction par pixels_total_source
    Signal_el_corrected = Signal_el * pixels_total_source
    sky_corrected = sky * pixels_total_source

    print(f"\n{'='*80}")
    print(f"AVEC CORRECTION × pixels_total_source ({pixels_total_source:.2f})")
    print(f"{'='*80}")
    print(f"  Signal: {Signal_el_corrected:.6e} / {signal_total_astropy:.6e} = {Signal_el_corrected/signal_total_astropy:.4f}")
    print(f"  Sky: {sky_corrected:.6e} / {sky_total_astropy:.6e} = {sky_corrected/sky_total_astropy:.4f}")

    # Renormalisation des inputs
    area_ratio_signal = source_size_arcsec_after_slit / (pixel_scale**2)
    area_ratio_sky = slit_size_arcsec_after_slit / (pixel_scale**2)

    spectral_ratio_signal = np.minimum(Line_width, Bandwidth) / dispersion
    spectral_ratio_sky = np.maximum(np.minimum(Line_width, Bandwidth), dispersion) / dispersion

    print(f"\n{'='*80}")
    print(f"RATIOS POUR RENORMALISATION")
    print(f"{'='*80}")
    print(f"  Spatial:")
    print(f"    source_size / pixel²: {area_ratio_signal:.2f}")
    print(f"    slit_size / pixel²: {area_ratio_sky:.2f}")
    print(f"  Spectral:")
    print(f"    min(Line_width, BW) / dispersion: {spectral_ratio_signal:.2f}")
    print(f"    max(min(LW, BW), disp) / dispersion: {spectral_ratio_sky:.2f}")
    print(f"  Total:")
    print(f"    Signal: {area_ratio_signal * spectral_ratio_signal:.2f}")
    print(f"    Sky: {area_ratio_sky * spectral_ratio_sky:.2f}")

    # Renormaliser
    Signal_renorm = Signal / (area_ratio_signal * spectral_ratio_signal) * pixels_total_source
    Sky_renorm = Sky / (area_ratio_sky * spectral_ratio_sky) * pixels_total_source

    # Recalculer
    Signal_LU_renorm = Signal_renorm / E_photon * (np.pi/180/3600)**2
    Sky_CU_renorm = Sky_renorm / E_photon * (np.pi/180/3600)**2

    Signal_el_renorm = Signal_LU_renorm * factor_CU2el_tot * t_exp * flux_fraction_slit_applied
    sky_renorm = Sky_CU_renorm * factor_CU2el_sky_tot * t_exp

    Signal_el_renorm_final = Signal_el_renorm * pixels_total_source
    sky_renorm_final = sky_renorm * pixels_total_source

    print(f"\n{'='*80}")
    print(f"AVEC RENORMALISATION DES INPUTS")
    print(f"{'='*80}")
    print(f"  Signal_renorm: {Signal_renorm:.6e} (÷ {area_ratio_signal * spectral_ratio_signal:.2f} × {pixels_total_source:.2f})")
    print(f"  Sky_renorm: {Sky_renorm:.6e} (÷ {area_ratio_sky * spectral_ratio_sky:.2f} × {pixels_total_source:.2f})")
    print(f"\n  Signal_el final: {Signal_el_renorm_final:.6e}")
    print(f"  Astropy: {signal_total_astropy:.6e}")
    print(f"  Ratio: {Signal_el_renorm_final/signal_total_astropy:.4f}")
    print(f"\n  sky final: {sky_renorm_final:.6e}")
    print(f"  Astropy: {sky_total_astropy:.6e}")
    print(f"  Ratio: {sky_renorm_final/sky_total_astropy:.4f}")

    # CONCLUSION
    print(f"\n{'='*80}")
    print(f"CONCLUSION")
    print(f"{'='*80}")

    if abs(Signal_el_renorm_final/signal_total_astropy - 1) < 0.15:
        print(f"✓ Signal MATCH après renormalisation!")
    else:
        print(f"✗ Signal ne match pas: {Signal_el_renorm_final/signal_total_astropy:.4f}")

    if abs(sky_renorm_final/sky_total_astropy - 1) < 0.15:
        print(f"✓ Sky MATCH après renormalisation!")
    else:
        print(f"✗ Sky ne match pas: {sky_renorm_final/sky_total_astropy:.4f}")

    return {
        'Signal_el': Signal_el,
        'sky': sky,
        'pixels_total_source': pixels_total_source,
        'signal_astropy': signal_total_astropy,
        'sky_astropy': sky_total_astropy,
        'area_ratio_signal': area_ratio_signal,
        'area_ratio_sky': area_ratio_sky,
    }

def test_sky_renormalization():
    """
    Si Generic ETC intègre le sky sur slit_size_arcsec_after_slit,
    et Astropy intègre sur pixel_scale²,
    alors il faut renormaliser le Sky d'entrée pour comparer !
    """

    try:
        import numpy as np
        from Observation import Observation, load_instruments
        from astropy.stats import signal_to_noise_oir_ccd

        instruments, database = load_instruments()
        instrument_name = "GALEX FUV"

        params = {}
        for i, charact in enumerate(instruments['Charact.']):
            if charact and not isinstance(charact, np.ma.core.MaskedConstant):
                value = instruments[instrument_name][i]
                if not isinstance(value, np.ma.core.MaskedConstant):
                    params[charact] = value

        print("="*80)
        print("TEST: Sky renormalization hypothesis")
        print("="*80)

        # Paramètres
        t_exp = 1000.0
        wavelength_nm = params['wavelength']
        E_photon = 1.986e-8 / wavelength_nm

        collecting_area_cm2 = params['Collecting_area'] * 1e4
        pixel_area_arcsec2 = params['pixel_scale']**2
        wavelength_range = params['dispersion']
        total_throughput = params['Throughput'] * params['QE'] * params['Atmosphere']

        # TEST 1: Sans renormalisation (ce qu'on fait actuellement)
        print("\n" + "="*80)
        print("TEST 1: Sans renormalisation")
        print("="*80)

        obs = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=t_exp, SNR_res="per pix", IFS=False, test=True
        )

        signal_flux = params['Signal']
        sky_flux = params['Sky']

        signal_photons = signal_flux / E_photon
        sky_photons = sky_flux / E_photon

        source_eps = (signal_photons * collecting_area_cm2 * pixel_area_arcsec2 *
                     wavelength_range * total_throughput)
        sky_eps = (sky_photons * collecting_area_cm2 * pixel_area_arcsec2 *
                  wavelength_range * total_throughput)
        dark_eps = params['Dark_current'] / 3600.0

        snr_astropy = signal_to_noise_oir_ccd(
            t=t_exp, source_eps=source_eps, sky_eps=sky_eps,
            dark_eps=dark_eps, rd=params['RN'], npix=1, gain=1.0
        )

        print(f"\nGeneric ETC:")
        print(f"  Signal_el: {obs.Signal_el[obs.i]:.6e} e⁻")
        print(f"  sky: {obs.sky[obs.i]:.6e} e⁻")
        print(f"  SNR: {obs.SNR[obs.i]:.6e}")
        print(f"  pixels_total_source: {obs.pixels_total_source:.2f}")

        print(f"\nAstropy:")
        print(f"  signal_total: {source_eps * t_exp:.6e} e⁻")
        print(f"  sky_total: {sky_eps * t_exp:.6e} e⁻")
        print(f"  SNR: {snr_astropy:.6e}")

        # Correction simple par pixels_total_source
        signal_corrected = obs.Signal_el[obs.i] * obs.pixels_total_source
        sky_corrected = obs.sky[obs.i] * obs.pixels_total_source

        print(f"\nRatios (Generic × npix / Astropy):")
        print(f"  Signal: {signal_corrected / (source_eps * t_exp):.4f}")
        print(f"  Sky: {sky_corrected / (sky_eps * t_exp):.4f}")

        # TEST 2: Avec renormalisation du Sky d'entrée
        print("\n" + "="*80)
        print("TEST 2: Avec renormalisation du Sky d'entrée")
        print("="*80)

        # Calculer le ratio des aires
        print(f"\nAires utilisées:")
        print(f"  pixel_scale²: {pixel_area_arcsec2:.2f} arcsec²")
        print(f"  slit_size_arcsec_after_slit: {obs.slit_size_arcsec_after_slit:.2f} arcsec²")
        print(f"  source_size_arcsec_after_slit: {obs.source_size_arcsec_after_slit:.2f} arcsec²")

        # Generic ETC utilise slit_size pour le sky
        area_ratio_sky = obs.slit_size_arcsec_after_slit / pixel_area_arcsec2
        area_ratio_signal = obs.source_size_arcsec_after_slit / pixel_area_arcsec2

        print(f"\nRatios d'aire:")
        print(f"  slit_size / pixel² = {area_ratio_sky:.2f}")
        print(f"  source_size / pixel² = {area_ratio_signal:.2f}")

        # Renormaliser le Sky d'entrée pour Generic ETC
        # Si Generic intègre sur slit_size et Astropy sur pixel²,
        # alors Sky_generic_input = Sky_astropy × (pixel² / slit_size)

        sky_flux_renormalized = sky_flux * (pixel_area_arcsec2 / obs.slit_size_arcsec_after_slit)

        print(f"\nSky renormalization:")
        print(f"  Original Sky: {sky_flux:.2e} erg/cm²/s/arcsec²/Å")
        print(f"  Renormalized Sky: {sky_flux_renormalized:.2e} (÷ {area_ratio_sky:.2f})")

        # Modifier l'instrument
        idx_sky = list(instruments['Charact.']).index('Sky')
        original_sky = instruments[instrument_name][idx_sky]
        instruments[instrument_name][idx_sky] = sky_flux_renormalized

        obs2 = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=t_exp, SNR_res="per pix", IFS=False, test=True
        )

        # Restore
        instruments[instrument_name][idx_sky] = original_sky

        print(f"\nGeneric ETC (avec Sky renormalisé):")
        print(f"  sky: {obs2.sky[obs2.i]:.6e} e⁻")
        print(f"  sky × npix: {obs2.sky[obs2.i] * obs2.pixels_total_source:.6e} e⁻")

        print(f"\nAstropy (inchangé):")
        print(f"  sky_total: {sky_eps * t_exp:.6e} e⁻")

        print(f"\nRatio (Generic renormalisé × npix / Astropy):")
        sky_ratio_renormalized = (obs2.sky[obs2.i] * obs2.pixels_total_source) / (sky_eps * t_exp)
        print(f"  Sky: {sky_ratio_renormalized:.4f}")

        if abs(sky_ratio_renormalized - 1.0) < 0.15:
            print(f"  ✓ MATCH ! La renormalisation fonctionne !")
        else:
            print(f"  ✗ Encore un problème...")

        # TEST 3: Faire la même chose pour le signal
        print("\n" + "="*80)
        print("TEST 3: Renormalisation du Signal aussi")
        print("="*80)

        signal_flux_renormalized = signal_flux * (pixel_area_arcsec2 / obs.source_size_arcsec_after_slit)

        print(f"\nSignal renormalization:")
        print(f"  Original Signal: {signal_flux:.2e} erg/cm²/s/arcsec²/Å")
        print(f"  Renormalized Signal: {signal_flux_renormalized:.2e} (÷ {area_ratio_signal:.2f})")

        idx_signal = list(instruments['Charact.']).index('Signal')
        original_signal = instruments[instrument_name][idx_signal]
        instruments[instrument_name][idx_signal] = signal_flux_renormalized
        instruments[instrument_name][idx_sky] = sky_flux_renormalized

        obs3 = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=t_exp, SNR_res="per pix", IFS=False, test=True
        )

        # Restore
        instruments[instrument_name][idx_signal] = original_signal
        instruments[instrument_name][idx_sky] = original_sky

        print(f"\nGeneric ETC (Signal ET Sky renormalisés):")
        print(f"  Signal_el × npix: {obs3.Signal_el[obs3.i] * obs3.pixels_total_source:.6e} e⁻")
        print(f"  sky × npix: {obs3.sky[obs3.i] * obs3.pixels_total_source:.6e} e⁻")

        print(f"\nAstropy:")
        print(f"  signal_total: {source_eps * t_exp:.6e} e⁻")
        print(f"  sky_total: {sky_eps * t_exp:.6e} e⁻")

        signal_ratio_final = (obs3.Signal_el[obs3.i] * obs3.pixels_total_source) / (source_eps * t_exp)
        sky_ratio_final = (obs3.sky[obs3.i] * obs3.pixels_total_source) / (sky_eps * t_exp)

        print(f"\nRatios finaux:")
        print(f"  Signal: {signal_ratio_final:.4f}")
        print(f"  Sky: {sky_ratio_final:.4f}")

        # CONCLUSION
        print("\n" + "="*80)
        print("CONCLUSION")
        print("="*80)

        if abs(signal_ratio_final - 1.0) < 0.15 and abs(sky_ratio_final - 1.0) < 0.15:
            print("\n✓✓✓ HYPOTHÈSE CONFIRMÉE ✓✓✓")
            print("\nLe problème n'est PAS dans le calcul du SNR !")
            print("Le problème est dans l'INTERPRÉTATION des paramètres d'entrée :")
            print(f"  - Generic ETC intègre Signal sur {obs.source_size_arcsec_after_slit:.1f} arcsec²")
            print(f"  - Generic ETC intègre Sky sur {obs.slit_size_arcsec_after_slit:.1f} arcsec²")
            print(f"  - Astropy intègre sur {pixel_area_arcsec2:.2f} arcsec² (pixel²)")
            print(f"\nPour comparer, il faut:")
            print(f"  Sky_generic = Sky_astropy × {pixel_area_arcsec2 / obs.slit_size_arcsec_after_slit:.4f}")
            print(f"  Signal_generic = Signal_astropy × {pixel_area_arcsec2 / obs.source_size_arcsec_after_slit:.4f}")
        else:
            print("\n✗ Il reste encore un problème...")
            print(f"  Signal ratio: {signal_ratio_final:.4f} (attendu: 1.0)")
            print(f"  Sky ratio: {sky_ratio_final:.4f} (attendu: 1.0)")

    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

def simple_snr_comparison():
    """
    Test simple : comparer directement les SNR en injectant les valeurs de Generic ETC dans Astropy
    """
    import numpy as np
    from Observation import Observation, load_instruments
    from astropy.stats import signal_to_noise_oir_ccd

    instruments, database = load_instruments()
    instrument_name = "GALEX FUV"

    params = {}
    for i, charact in enumerate(instruments['Charact.']):
        if charact and not isinstance(charact, np.ma.core.MaskedConstant):
            value = instruments[instrument_name][i]
            if not isinstance(value, np.ma.core.MaskedConstant):
                params[charact] = value

    print("="*80)
    print("TEST SIMPLE : Comparaison SNR")
    print("="*80)

    t_exp = 1000.0

    # Generic ETC
    obs = Observation(
        instruments=instruments, instrument=instrument_name,
        exposure_time=t_exp, SNR_res="per pix", IFS=False, test=True
    )

    print(f"\nGeneric ETC:")
    print(f"  Signal_el: {obs.Signal_el[obs.i]:.6e} e⁻")
    print(f"  sky: {obs.sky[obs.i]:.6e} e⁻")
    print(f"  Dark_current_f: {obs.Dark_current_f[obs.i]:.6e} e⁻")
    print(f"  number_pixels_used: {obs.number_pixels_used}")
    print(f"  pixels_total_source: {obs.pixels_total_source}")
    print(f"  SNR: {obs.SNR[obs.i]:.6e}")

    # Convertir en taux pour Astropy (e⁻ → e⁻/s)
    source_eps_generic = obs.Signal_el[obs.i] / t_exp
    sky_eps_generic = obs.sky[obs.i] / t_exp
    dark_eps_generic = obs.Dark_current_f[obs.i] / t_exp

    print(f"\nTaux (Generic ETC):")
    print(f"  source_eps: {source_eps_generic:.6e} e⁻/s")
    print(f"  sky_eps: {sky_eps_generic:.6e} e⁻/s")
    print(f"  dark_eps: {dark_eps_generic:.6e} e⁻/s")

    # Astropy avec npix=1 (comme Generic qui a number_pixels_used=1 en mode "per pix")
    snr_astropy = signal_to_noise_oir_ccd(
        t=t_exp,
        source_eps=source_eps_generic,
        sky_eps=sky_eps_generic,
        dark_eps=dark_eps_generic,
        rd=params['RN'],
        npix=int(obs.number_pixels_used),
        gain=1.0
    )

    print(f"\nAstropy (avec taux Generic):")
    print(f"  SNR: {snr_astropy:.6e}")

    print(f"\n{'='*80}")
    print(f"COMPARAISON SNR")
    print(f"{'='*80}")
    print(f"  Generic: {obs.SNR[obs.i]:.6e}")
    print(f"  Astropy: {snr_astropy:.6e}")
    print(f"  Ratio: {obs.SNR[obs.i]/snr_astropy:.4f}")

    if abs(obs.SNR[obs.i]/snr_astropy - 1) < 0.15:
        print(f"\n✓ SNR MATCH !")
    else:
        print(f"\n✗ SNR différent : {obs.SNR[obs.i]/snr_astropy:.4f}")

    # Test variations
    print(f"\n{'='*80}")
    print(f"TEST VARIATIONS")
    print(f"{'='*80}")

    results = []

    # Variation exposure time
    for t in [100, 500, 1000, 2000, 5000]:
        obs = Observation(instruments=instruments, instrument=instrument_name,
                         exposure_time=t, SNR_res="per pix", IFS=False, test=True)

        source_eps = obs.Signal_el[obs.i] / t
        sky_eps = obs.sky[obs.i] / t
        dark_eps = obs.Dark_current_f[obs.i] / t

        snr_astropy = signal_to_noise_oir_ccd(
            t=t, source_eps=source_eps, sky_eps=sky_eps, dark_eps=dark_eps,
            rd=params['RN'], npix=int(obs.number_pixels_used), gain=1.0
        )

        results.append({
            'param': 'exp_time',
            'value': t,
            'snr_generic': obs.SNR[obs.i],
            'snr_astropy': snr_astropy,
            'ratio': obs.SNR[obs.i]/snr_astropy
        })

    # Variation sky
    sky_factors = [0.1, 0.5, 1.0, 2.0, 5.0]
    for factor in sky_factors:
        idx = list(instruments['Charact.']).index('Sky')
        original = instruments[instrument_name][idx]
        instruments[instrument_name][idx] = params['Sky'] * factor

        obs = Observation(instruments=instruments, instrument=instrument_name,
                         exposure_time=1000, SNR_res="per pix", IFS=False, test=True)

        source_eps = obs.Signal_el[obs.i] / 1000
        sky_eps = obs.sky[obs.i] / 1000
        dark_eps = obs.Dark_current_f[obs.i] / 1000

        snr_astropy = signal_to_noise_oir_ccd(
            t=1000, source_eps=source_eps, sky_eps=sky_eps, dark_eps=dark_eps,
            rd=params['RN'], npix=int(obs.number_pixels_used), gain=1.0
        )

        instruments[instrument_name][idx] = original

        results.append({
            'param': 'sky',
            'value': factor,
            'snr_generic': obs.SNR[obs.i],
            'snr_astropy': snr_astropy,
            'ratio': obs.SNR[obs.i]/snr_astropy
        })

    # Afficher résultats
    import pandas as pd
    df = pd.DataFrame(results)

    print("\nVariation exposure time:")
    print(df[df['param']=='exp_time'][['value', 'snr_generic', 'snr_astropy', 'ratio']])

    print("\nVariation sky:")
    print(df[df['param']=='sky'][['value', 'snr_generic', 'snr_astropy', 'ratio']])

    print(f"\n{'='*80}")
    print(f"CONCLUSION")
    print(f"{'='*80}")

    ratio_std = df['ratio'].std()
    ratio_mean = df['ratio'].mean()

    print(f"Ratio moyen: {ratio_mean:.4f}")
    print(f"Écart-type: {ratio_std:.4f}")

    if ratio_std < 0.01:
        print(f"✓ Ratio constant → les variations sont cohérentes")
    else:
        print(f"✗ Ratio varie → problème dans le calcul")

    if abs(ratio_mean - 1.0) < 0.15:
        print(f"✓ Ratio proche de 1 → calculs cohérents")
    else:
        print(f"✗ Ratio = {ratio_mean:.4f} → différence systématique")

if __name__ == "__main__":
    simple_snr_comparison()
