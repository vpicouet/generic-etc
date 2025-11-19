"""
Script de debug pour identifier le problème de cross-validation
À exécuter dans un environnement avec numpy/astropy installés
"""

import sys
sys.path.insert(0, 'notebooks')

def analyze_cross_validation():
    """Analyse le problème de cross-validation"""

    try:
        import numpy as np
        from Observation import Observation, load_instruments
        from astropy.stats import signal_to_noise_oir_ccd

        # Load instruments
        instruments, database = load_instruments()
        instrument_name = "GALEX FUV"
        t_exp = 1000.0

        # Generic ETC
        obs = Observation(
            instruments=instruments,
            instrument=instrument_name,
            exposure_time=t_exp,
            SNR_res="per pix",
            IFS=False,
            test=True
        )

        print("="*80)
        print("DEBUG: Valeurs internes Generic ETC")
        print("="*80)

        # Extraire params
        params = {}
        for i, charact in enumerate(instruments['Charact.']):
            if charact and not isinstance(charact, np.ma.core.MaskedConstant):
                value = instruments[instrument_name][i]
                if not isinstance(value, np.ma.core.MaskedConstant):
                    params[charact] = value

        print(f"\n1. PIXELS ET DIMENSIONS:")
        print(f"   pixels_total_source = {obs.pixels_total_source}")
        print(f"   elem_size = {obs.elem_size}")
        print(f"   number_pixels_used = {obs.number_pixels_used}")
        print(f"   source_size = {obs.source_size}")

        print(f"\n2. FACTEURS DE CONVERSION (CU/LU → e⁻):")
        print(f"   factor_CU2el = {obs.factor_CU2el:.6e}")
        print(f"   factor_CU2el_sky = {obs.factor_CU2el_sky:.6e}")
        print(f"   Ratio (sky/signal) = {obs.factor_CU2el_sky/obs.factor_CU2el:.2f}")

        # Calculer ce que serait le facteur SANS division par pixels_total_source
        effective_area = params['Collecting_area'] * 1e4  # m² → cm²
        throughput_total = params['Throughput'] * params['QE'] * params['Atmosphere']
        effective_area_total = effective_area * throughput_total
        arcsec2str = (np.pi/180/3600)**2
        pixel_area = params['pixel_scale']**2
        dispersion = params['dispersion']

        factor_per_pix_expected = effective_area_total * arcsec2str * pixel_area * dispersion

        print(f"\n3. FACTEUR ATTENDU (sans division par pixels_total_source):")
        print(f"   factor_per_pix_expected = {factor_per_pix_expected:.6e}")
        print(f"   factor_CU2el * pixels_total_source = {obs.factor_CU2el * obs.pixels_total_source:.6e}")
        print(f"   Ratio = {(obs.factor_CU2el * obs.pixels_total_source) / factor_per_pix_expected:.4f}")

        print(f"\n4. SIGNAL (e⁻/pix/exp):")
        print(f"   Signal_el (Generic) = {obs.Signal_el[obs.i]:.6e}")
        print(f"   Signal_el * pixels_total_source = {obs.Signal_el[obs.i] * obs.pixels_total_source:.6e}")

        print(f"\n5. SKY (e⁻/pix/exp):")
        print(f"   sky (Generic) = {obs.sky[obs.i]:.6e}")
        print(f"   Sky_noise (Generic) = {obs.Sky_noise[obs.i]:.6e}")
        print(f"   sky * pixels_total_source = {obs.sky[obs.i] * obs.pixels_total_source:.6e}")

        print(f"\n6. DARK (e⁻/pix/exp):")
        print(f"   Dark_current_f = {obs.Dark_current_f[obs.i]:.6e}")

        print(f"\n7. TAILLES ET AIRES:")
        print(f"   source_size_arcsec_after_slit = {obs.source_size_arcsec_after_slit} arcsec²")
        print(f"   slit_size_arcsec_after_slit = {obs.slit_size_arcsec_after_slit} arcsec²")
        print(f"   pixel_area = {pixel_area} arcsec²")
        print(f"   Ratio (slit/pixel) = {obs.slit_size_arcsec_after_slit / pixel_area:.2f}")

        # Calcul Astropy
        wavelength_nm = params['wavelength']
        E_photon = 1.986e-8 / wavelength_nm

        signal_flux = params['Signal']
        sky_flux = params['Sky']

        signal_photons = signal_flux / E_photon
        sky_photons = sky_flux / E_photon

        collecting_area_cm2 = params['Collecting_area'] * 1e4

        # OPTION 1: Par pixel (ce qu'Astropy attend)
        source_eps_per_pix = (signal_photons * collecting_area_cm2 *
                              pixel_area * dispersion * throughput_total)
        sky_eps_per_pix = (sky_photons * collecting_area_cm2 *
                          pixel_area * dispersion * throughput_total)
        dark_eps = params['Dark_current'] / 3600.0

        signal_total_per_pix = source_eps_per_pix * t_exp
        sky_total_per_pix = sky_eps_per_pix * t_exp
        dark_total = dark_eps * t_exp

        print(f"\n{'='*80}")
        print("ASTROPY (par pixel individuel)")
        print("="*80)
        print(f"   source_eps = {source_eps_per_pix:.6e} e⁻/s/pix")
        print(f"   sky_eps = {sky_eps_per_pix:.6e} e⁻/s/pix")
        print(f"   dark_eps = {dark_eps:.6e} e⁻/s/pix")
        print(f"\n   Signal total = {signal_total_per_pix:.6e} e⁻")
        print(f"   Sky total = {sky_total_per_pix:.6e} e⁻")
        print(f"   Dark total = {dark_total:.6e} e⁻")

        snr_astropy = signal_to_noise_oir_ccd(
            t=t_exp,
            source_eps=source_eps_per_pix,
            sky_eps=sky_eps_per_pix,
            dark_eps=dark_eps,
            rd=params['RN'],
            npix=1,
            gain=1.0
        )
        print(f"   SNR = {snr_astropy:.6e}")

        print(f"\n{'='*80}")
        print("COMPARAISON")
        print("="*80)
        print(f"\nSignal:")
        print(f"   Generic = {obs.Signal_el[obs.i]:.6e}")
        print(f"   Astropy = {signal_total_per_pix:.6e}")
        print(f"   Ratio Gen/Ast = {obs.Signal_el[obs.i]/signal_total_per_pix:.4f}")
        print(f"   Generic * pixels_total_source = {obs.Signal_el[obs.i] * obs.pixels_total_source:.6e}")
        print(f"   Ratio (Gen*pix)/Ast = {(obs.Signal_el[obs.i] * obs.pixels_total_source)/signal_total_per_pix:.4f}")

        print(f"\nSky:")
        print(f"   Generic = {obs.sky[obs.i]:.6e}")
        print(f"   Astropy = {sky_total_per_pix:.6e}")
        print(f"   Ratio Gen/Ast = {obs.sky[obs.i]/sky_total_per_pix:.4f}")
        print(f"   Generic * pixels_total_source = {obs.sky[obs.i] * obs.pixels_total_source:.6e}")
        print(f"   Ratio (Gen*pix)/Ast = {(obs.sky[obs.i] * obs.pixels_total_source)/sky_total_per_pix:.4f}")

        print(f"\nDark:")
        print(f"   Generic = {obs.Dark_current_f[obs.i]:.6e}")
        print(f"   Astropy = {dark_total:.6e}")
        print(f"   Ratio Gen/Ast = {obs.Dark_current_f[obs.i]/dark_total:.4f}")

        print(f"\nSNR:")
        print(f"   Generic = {obs.SNR[obs.i]:.6e}")
        print(f"   Astropy = {snr_astropy:.6e}")
        print(f"   Ratio Gen/Ast = {obs.SNR[obs.i]/snr_astropy:.4f}")

        print(f"\n{'='*80}")
        print("DIAGNOSTIC")
        print("="*80)

        signal_match_with_correction = abs((obs.Signal_el[obs.i] * obs.pixels_total_source) / signal_total_per_pix - 1) < 0.1
        sky_match_with_correction = abs((obs.sky[obs.i] * obs.pixels_total_source) / sky_total_per_pix - 1) < 0.1

        if signal_match_with_correction:
            print("✓ Signal: OK après multiplication par pixels_total_source")
            print("  → Generic ETC calcule le signal MOYEN par pixel")
            print("  → Correction: multiplier par pixels_total_source pour comparaison")
        else:
            print("✗ Signal: Problème persiste même après correction")
            print(f"  → Facteur restant: {(obs.Signal_el[obs.i] * obs.pixels_total_source) / signal_total_per_pix:.2f}")

        if sky_match_with_correction:
            print("\n✓ Sky: OK après multiplication par pixels_total_source")
            print("  → Generic ETC calcule le sky MOYEN par pixel")
            print("  → Correction: multiplier par pixels_total_source pour comparaison")
        else:
            print("\n✗ Sky: Problème persiste même après correction")
            print(f"  → Facteur restant: {(obs.sky[obs.i] * obs.pixels_total_source) / sky_total_per_pix:.2f}")
            print(f"  → Vérifier factor_CU2el_sky vs factor_CU2el")
            print(f"  → Ratio factors: {obs.factor_CU2el_sky / obs.factor_CU2el:.2f}")

        return obs, params

    except ImportError as e:
        print(f"Erreur: modules manquants - {e}")
        print("\nInstallez les dépendances avec:")
        print("  pip install numpy astropy pandas")
        return None, None

if __name__ == "__main__":
    analyze_cross_validation()
