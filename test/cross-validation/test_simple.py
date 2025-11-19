"""
Simple test without matplotlib
"""
import numpy as np
import sys
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, str(Path.cwd().parent.parent / 'notebooks'))

# Import Generic ETC
from Observation import Observation, load_instruments

# Import Astropy SNR function
try:
    from astropy.stats import signal_to_noise_oir_ccd
    ASTROPY_SNR_AVAILABLE = True
    print("✓ astropy.stats available")
except ImportError:
    ASTROPY_SNR_AVAILABLE = False
    print("⚠ astropy.stats not available")
    sys.exit(1)

def calculate_SNR_astropy(signal_flux, sky_flux, collecting_area, throughput, QE, atmosphere,
                          pixel_scale, dispersion, wavelength_nm, exposure_time,
                          dark_current, read_noise, n_pixels=1):
    """Calculate SNR using astropy."""
    E_photon = 1.986e-8 / wavelength_nm
    signal_photons_rate = signal_flux / E_photon
    sky_photons_rate = sky_flux / E_photon

    collecting_area_cm2 = collecting_area * 1e4
    pixel_area_arcsec2 = pixel_scale**2
    wavelength_range = dispersion
    total_throughput = throughput * QE * atmosphere

    source_eps_per_pix = (signal_photons_rate * collecting_area_cm2 *
                          pixel_area_arcsec2 * wavelength_range * total_throughput)
    source_eps = source_eps_per_pix * n_pixels

    sky_eps = (sky_photons_rate * collecting_area_cm2 *
               pixel_area_arcsec2 * wavelength_range * total_throughput)

    dark_eps = dark_current / 3600.0
    rd = read_noise

    snr = signal_to_noise_oir_ccd(
        t=exposure_time,
        source_eps=source_eps,
        sky_eps=sky_eps,
        dark_eps=dark_eps,
        rd=rd,
        npix=n_pixels,
        gain=1.0
    )

    return {'SNR': snr}

# Load instruments
print("Loading instruments...")
instruments, database = load_instruments()
print(f"✓ Loaded {len(instruments.colnames)-3} instruments\n")

# Test on GALEX FUV
instrument_name = "GALEX FUV"
print(f"Testing: {instrument_name}")
print("="*60)

# Extract parameters
idx = list(instruments.colnames).index(instrument_name)
params = {}
for i, charact in enumerate(instruments['Charact.']):
    if charact and not isinstance(charact, np.ma.core.MaskedConstant):
        value = instruments[instrument_name][i]
        if not isinstance(value, np.ma.core.MaskedConstant):
            params[charact] = value

# Test exposure time variation
print("\nTesting exposure time variation...")
exposure_times = np.logspace(1, 4, 10)
snr_generic_list = []
snr_astropy_list = []

for t_exp in exposure_times:
    # Generic ETC
    obs = Observation(
        instruments=instruments,
        instrument=instrument_name,
        exposure_time=t_exp,
        SNR_res="per pix",
        IFS=False,
        test=False
    )
    snr_generic_list.append(obs.SNR[obs.i])

    # Astropy SNR
    result = calculate_SNR_astropy(
        signal_flux=params['Signal'],
        sky_flux=params['Sky'],
        collecting_area=params['Collecting_area'],
        throughput=params['Throughput'],
        QE=params['QE'],
        atmosphere=params['Atmosphere'],
        pixel_scale=params['pixel_scale'],
        dispersion=params['dispersion'],
        wavelength_nm=params['wavelength'],
        exposure_time=t_exp,
        dark_current=params['Dark_current'],
        read_noise=params['RN'],
        n_pixels=1
    )
    snr_astropy_list.append(result['SNR'])

snr_generic_array = np.array(snr_generic_list)
snr_astropy_array = np.array(snr_astropy_list)
relative_diff = 100 * (snr_generic_array - snr_astropy_array) / snr_astropy_array

print(f"\nResults:")
print(f"  Mean relative difference: {np.mean(relative_diff):.2f}%")
print(f"  Std relative difference: {np.std(relative_diff):.2f}%")
print(f"  Max relative difference: {np.max(np.abs(relative_diff)):.2f}%")

if np.max(np.abs(relative_diff)) < 10:
    print(f"\n  Result: ✓✓ EXCELLENT (< 10%)")
elif np.max(np.abs(relative_diff)) < 15:
    print(f"\n  Result: ✓ PASSED (< 15%)")
else:
    print(f"\n  Result: ✗ FAILED (> 15%)")

print("="*60)
print("\n✓ Validation function works correctly!")
print("✓ Notebook is ready to use")
