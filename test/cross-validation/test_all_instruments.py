#!/usr/bin/env python3
"""
Test cross-validation Generic ETC vs Astropy pour différents instruments
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
    # Already in notebooks or subdirectory
    sys.path.insert(0, 'notebooks')

import numpy as np
from Observation import Observation, load_instruments
from astropy.stats import signal_to_noise_oir_ccd

def test_instrument(instruments, instrument_name, verbose=False):
    """
    Test cross-validation for a single instrument
    Returns: dict with test results
    """
    # Extract parameters
    params = {}
    for i, charact in enumerate(instruments['Charact.']):
        if charact and not isinstance(charact, np.ma.core.MaskedConstant):
            value = instruments[instrument_name][i]
            if not isinstance(value, np.ma.core.MaskedConstant):
                params[charact] = value

    # Test with default exposure time
    t_exp = 1000.0

    try:
        # Get readout time for this instrument
        readout_time = params.get('Readout', 0.0)

        # Generic ETC - set acquisition_time so that N_images = 1
        # N_images = acquisition_time * 3600 / (exposure_time + readout_time)
        # For N_images = 1: acquisition_time = (exposure_time + readout_time) / 3600
        acq_time = (t_exp + readout_time) / 3600.0

        obs = Observation(
            instruments=instruments, instrument=instrument_name,
            exposure_time=t_exp, SNR_res="per pix", IFS=False, test=True,
            acquisition_time=acq_time
        )

        # Convert to rates
        source_eps = obs.Signal_el[obs.i] / t_exp
        sky_eps = obs.sky[obs.i] / t_exp
        dark_eps = obs.Dark_current_f[obs.i] / t_exp

        # Astropy
        snr_astropy = signal_to_noise_oir_ccd(
            t=t_exp,
            source_eps=source_eps,
            sky_eps=sky_eps,
            dark_eps=dark_eps,
            rd=params['RN'],
            npix=int(obs.number_pixels_used),
            gain=1.0
        )

        ratio = obs.SNR[obs.i] / snr_astropy

        if verbose:
            print(f"\n{instrument_name}:")
            print(f"  SNR Generic: {obs.SNR[obs.i]:.6e}")
            print(f"  SNR Astropy: {snr_astropy:.6e}")
            print(f"  Ratio: {ratio:.6f}")
            print(f"  N_images_true: {obs.N_images_true}")

        return {
            'instrument': instrument_name,
            'snr_generic': obs.SNR[obs.i],
            'snr_astropy': snr_astropy,
            'ratio': ratio,
            'n_images': obs.N_images_true,
            'success': True,
            'error': None
        }

    except Exception as e:
        if verbose:
            print(f"\n{instrument_name}: ERROR - {str(e)}")
        return {
            'instrument': instrument_name,
            'success': False,
            'error': str(e)
        }

def main():
    """
    Test cross-validation on multiple instruments
    """
    print("="*80)
    print("TEST CROSS-VALIDATION SUR PLUSIEURS INSTRUMENTS")
    print("="*80)

    instruments, database = load_instruments()

    # Get list of instruments (column names except metadata columns)
    metadata_cols = {'.', 'Charact.', 'Unit'}
    instrument_names = [key for key in instruments.keys() if key not in metadata_cols]

    print(f"\n{len(instrument_names)} instruments trouvés\n")

    # Test each instrument
    results = []
    for inst_name in instrument_names:
        result = test_instrument(instruments, inst_name, verbose=True)
        results.append(result)

    # Summary
    print("\n" + "="*80)
    print("RÉSUMÉ")
    print("="*80)

    import pandas as pd
    df = pd.DataFrame([r for r in results if r['success']])

    if len(df) > 0:
        print(f"\nInstruments testés avec succès: {len(df)}/{len(results)}")
        print(f"\nRatio moyen: {df['ratio'].mean():.6f}")
        print(f"Écart-type: {df['ratio'].std():.6f}")
        print(f"Min: {df['ratio'].min():.6f}")
        print(f"Max: {df['ratio'].max():.6f}")

        print("\nDétails par instrument:")
        print(df[['instrument', 'snr_generic', 'snr_astropy', 'ratio']].to_string(index=False))

        # Check for issues
        bad_ratios = df[(df['ratio'] < 0.95) | (df['ratio'] > 1.05)]
        if len(bad_ratios) > 0:
            print("\n⚠ Instruments avec ratio suspect (< 0.95 ou > 1.05):")
            print(bad_ratios[['instrument', 'ratio']].to_string(index=False))
        else:
            print("\n✓ Tous les instruments ont un ratio entre 0.95 et 1.05")

    # Failed tests
    failed = [r for r in results if not r['success']]
    if len(failed) > 0:
        print(f"\n❌ {len(failed)} instruments ont échoué:")
        for r in failed:
            print(f"  - {r['instrument']}: {r['error']}")

if __name__ == "__main__":
    main()
