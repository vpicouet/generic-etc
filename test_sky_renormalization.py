#!/usr/bin/env python3
"""
Test de l'hypothèse : le problème vient de la différence d'aire pour le sky
Si on renormalise le Sky d'entrée, tout devrait matcher !
"""

import sys
sys.path.insert(0, 'notebooks')

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

if __name__ == "__main__":
    test_sky_renormalization()
