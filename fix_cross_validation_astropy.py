"""
Corrections temporaires pour cross-validation avec Astropy
À utiliser dans le notebook de cross-validation
"""

def get_astropy_comparison_params(obs, t_exp, params):
    """
    Convertit les paramètres Generic ETC en paramètres Astropy compatibles

    Args:
        obs: Objet Observation de Generic ETC
        t_exp: Temps d'exposition (s)
        params: Dictionnaire des paramètres de l'instrument

    Returns:
        dict avec source_eps, sky_eps, dark_eps, signal_total, sky_total, dark_total
    """
    import numpy as np

    # Paramètres de base
    wavelength_nm = params['wavelength']
    E_photon = 1.986e-8 / wavelength_nm  # erg per photon

    signal_flux = params['Signal']  # erg/cm²/s/arcsec²/Å
    sky_flux = params['Sky']        # erg/cm²/s/arcsec²/Å

    # Conversion en photons
    signal_photons_rate = signal_flux / E_photon  # photons/cm²/s/arcsec²/Å
    sky_photons_rate = sky_flux / E_photon

    # Paramètres télescope/instrument
    collecting_area_cm2 = params['Collecting_area'] * 1e4  # m² → cm²
    pixel_area_arcsec2 = params['pixel_scale']**2
    wavelength_range = params['dispersion']  # Å/pix
    total_throughput = params['Throughput'] * params['QE'] * params['Atmosphere']

    # OPTION A: Taux par pixel individuel (ce qu'Astropy attend normalement)
    source_eps_per_pix = (signal_photons_rate *
                          collecting_area_cm2 *
                          pixel_area_arcsec2 *
                          wavelength_range *
                          total_throughput)

    sky_eps_per_pix = (sky_photons_rate *
                       collecting_area_cm2 *
                       pixel_area_arcsec2 *
                       wavelength_range *
                       total_throughput)

    dark_eps = params['Dark_current'] / 3600.0  # e⁻/pix/hour → e⁻/s/pix

    # OPTION B: Taux moyennés sur pixels_total_source (pour matcher Generic ETC)
    source_eps_averaged = source_eps_per_pix / obs.pixels_total_source
    sky_eps_averaged = sky_eps_per_pix / obs.pixels_total_source

    return {
        # Option A: Par pixel individuel
        'source_eps_per_pix': source_eps_per_pix,
        'sky_eps_per_pix': sky_eps_per_pix,
        'dark_eps': dark_eps,
        'signal_total_per_pix': source_eps_per_pix * t_exp,
        'sky_total_per_pix': sky_eps_per_pix * t_exp,
        'dark_total': dark_eps * t_exp,
        'npix_per_pix': 1,

        # Option B: Moyenné (pour matcher Generic ETC)
        'source_eps_averaged': source_eps_averaged,
        'sky_eps_averaged': sky_eps_averaged,
        'signal_total_averaged': source_eps_averaged * t_exp,
        'sky_total_averaged': sky_eps_averaged * t_exp,
        'npix_averaged': obs.pixels_total_source,

        # Métadonnées
        'pixels_total_source': obs.pixels_total_source,
    }


def compare_generic_vs_astropy(obs, t_exp, params, use_averaged=False):
    """
    Compare Generic ETC avec Astropy SNR

    Args:
        obs: Objet Observation
        t_exp: Temps d'exposition
        params: Paramètres instrument
        use_averaged: Si True, utilise les valeurs moyennées sur pixels_total_source

    Returns:
        dict avec résultats de comparaison
    """
    from astropy.stats import signal_to_noise_oir_ccd

    astropy_params = get_astropy_comparison_params(obs, t_exp, params)

    if use_averaged:
        # Comparer avec Generic ETC "moyen"
        source_eps = astropy_params['source_eps_averaged']
        sky_eps = astropy_params['sky_eps_averaged']
        npix = astropy_params['npix_averaged']
        signal_total = astropy_params['signal_total_averaged']
        sky_total = astropy_params['sky_total_averaged']
    else:
        # Comparer par pixel individuel (nécessite de corriger Generic ETC)
        source_eps = astropy_params['source_eps_per_pix']
        sky_eps = astropy_params['sky_eps_per_pix']
        npix = astropy_params['npix_per_pix']
        signal_total = astropy_params['signal_total_per_pix']
        sky_total = astropy_params['sky_total_per_pix']

    dark_eps = astropy_params['dark_eps']
    dark_total = astropy_params['dark_total']

    # Calcul SNR Astropy
    snr_astropy = signal_to_noise_oir_ccd(
        t=t_exp,
        source_eps=source_eps,
        sky_eps=sky_eps,
        dark_eps=dark_eps,
        rd=params['RN'],
        npix=npix,
        gain=1.0
    )

    # Valeurs Generic ETC
    if use_averaged:
        # Utiliser directement les valeurs Generic ETC (déjà moyennées)
        signal_generic = obs.Signal_el[obs.i]
        sky_generic = obs.sky[obs.i]
    else:
        # Corriger Generic ETC pour obtenir valeurs par pixel
        signal_generic = obs.Signal_el[obs.i] * obs.pixels_total_source
        sky_generic = obs.sky[obs.i] * obs.pixels_total_source

    dark_generic = obs.Dark_current_f[obs.i]
    snr_generic = obs.SNR[obs.i]

    return {
        'mode': 'averaged' if use_averaged else 'per_pixel',
        'pixels_total_source': obs.pixels_total_source,

        # Generic ETC
        'signal_generic': signal_generic,
        'sky_generic': sky_generic,
        'dark_generic': dark_generic,
        'snr_generic': snr_generic,

        # Astropy
        'signal_astropy': signal_total,
        'sky_astropy': sky_total,
        'dark_astropy': dark_total,
        'snr_astropy': snr_astropy,

        # Ratios
        'ratio_signal': signal_generic / signal_total if signal_total > 0 else np.nan,
        'ratio_sky': sky_generic / sky_total if sky_total > 0 else np.nan,
        'ratio_dark': dark_generic / dark_total if dark_total > 0 else np.nan,
        'ratio_snr': snr_generic / snr_astropy if snr_astropy > 0 else np.nan,
    }


def print_comparison(results):
    """Affiche les résultats de comparaison de manière claire"""

    print(f"\n{'='*80}")
    print(f"COMPARAISON - Mode: {results['mode']}")
    print(f"pixels_total_source = {results['pixels_total_source']}")
    print("="*80)

    print(f"\n{'Paramètre':<20} {'Generic ETC':>15} {'Astropy':>15} {'Ratio':>10}")
    print("-"*80)
    print(f"{'Signal (e⁻)':<20} {results['signal_generic']:>15.6e} {results['signal_astropy']:>15.6e} {results['ratio_signal']:>10.4f}")
    print(f"{'Sky (e⁻)':<20} {results['sky_generic']:>15.6e} {results['sky_astropy']:>15.6e} {results['ratio_sky']:>10.4f}")
    print(f"{'Dark (e⁻)':<20} {results['dark_generic']:>15.6e} {results['dark_astropy']:>15.6e} {results['ratio_dark']:>10.4f}")
    print(f"{'SNR':<20} {results['snr_generic']:>15.6e} {results['snr_astropy']:>15.6e} {results['ratio_snr']:>10.4f}")

    # Diagnostic
    print(f"\n{'='*80}")
    print("DIAGNOSTIC")
    print("="*80)

    tolerance = 0.15  # 15% tolerance

    checks = [
        ('Signal', results['ratio_signal']),
        ('Sky', results['ratio_sky']),
        ('Dark', results['ratio_dark']),
        ('SNR', results['ratio_snr']),
    ]

    for name, ratio in checks:
        if abs(ratio - 1.0) < tolerance:
            print(f"✓ {name}: OK (ratio = {ratio:.4f})")
        else:
            print(f"✗ {name}: Discrepancy (ratio = {ratio:.4f})")


# Exemple d'utilisation dans un notebook :
"""
from fix_cross_validation_astropy import compare_generic_vs_astropy, print_comparison
from Observation import Observation, load_instruments

instruments, database = load_instruments()
obs = Observation(instruments=instruments, instrument="GALEX FUV",
                  exposure_time=1000, SNR_res="per pix", IFS=False, test=True)

params = {charact: instruments["GALEX FUV"][i]
          for i, charact in enumerate(instruments['Charact.'])
          if not isinstance(charact, np.ma.core.MaskedConstant)}

# Test 1: Comparaison "moyennée" (devrait matcher)
results_avg = compare_generic_vs_astropy(obs, 1000, params, use_averaged=True)
print_comparison(results_avg)

# Test 2: Comparaison "per pixel" (avec correction Generic ETC)
results_pix = compare_generic_vs_astropy(obs, 1000, params, use_averaged=False)
print_comparison(results_pix)
"""
