"""
Simple SNR Calculator with Interactive Widgets

Demonstrates the relationship between:
- Charges (electrons accumulated)
- Noise (sigma of Poisson distributions)
- SNR evolution with exposure time, number of exposures, and pixels
"""

import numpy as np
import matplotlib.pyplot as plt
from ipywidgets import FloatSlider, IntSlider, interact, VBox, HBox, Layout
from IPython.display import display

def calculate_snr_vs_exptime(CIC=0.001, RN=0.01, dark_rate=0.1, signal_rate=0.02,
                             max_exptime=100, acquisition_time=1.0, n_pixels=1,
                             current_exptime=50, readout_time=0):
    """
    Calculate SNR vs exposure time for given detector parameters.

    Parameters
    ----------
    CIC : float
        Clock Induced Charge [e-/pix/frame] - constant per frame
    RN : float
        Read Noise [e-/pix/frame] - constant per frame
    dark_rate : float
        Dark current rate [e-/pix/hour]
    signal_rate : float
        Signal rate [e-/pix/hour]
    max_exptime : float
        Maximum exposure time for plot [seconds]
    acquisition_time : float
        Total acquisition time [hours]
    n_pixels : int
        Number of pixels to sum over
    current_exptime : float
        Current exposure time for vertical bar [seconds]
    readout_time : float
        Readout time between exposures [seconds]

    Returns
    -------
    exptime : array
        Exposure times [seconds]
    charges_per_exp : dict
        Charges per exposure per pixel
    noises_per_exp : dict
        Noises per exposure per pixel
    snr_per_exp : array
        SNR per exposure per pixel
    n_exposures : array
        Number of exposures for each exptime
    snr_total : array
        Total SNR (with n_exposures and n_pixels)
    """
    # Create exposure time array
    exptime = np.linspace(0.1, max_exptime, 200)  # seconds

    # Convert rates from e-/pix/hour to e-/pix/second
    dark_rate_per_sec = dark_rate / 3600
    signal_rate_per_sec = signal_rate / 3600

    # Calculate number of exposures for each exposure time
    acquisition_time_sec = acquisition_time * 3600  # convert hours to seconds
    n_exposures = acquisition_time_sec / (exptime + readout_time)

    # Calculate CHARGES per exposure per pixel (nombre d'electrons)
    signal_charges = signal_rate_per_sec * exptime  # e-/pix/exp
    dark_charges = dark_rate_per_sec * exptime      # e-/pix/exp
    CIC_charges = CIC * np.ones_like(exptime)       # e-/pix/exp (constant per frame)
    RN_charges = RN * np.ones_like(exptime)         # e-/pix/exp (constant per frame)

    # Calculate NOISES per exposure per pixel (sigma, bruit)
    signal_noise = np.sqrt(signal_charges)          # sigma [e-]
    dark_noise = np.sqrt(dark_charges)              # sigma [e-]
    CIC_noise = np.sqrt(CIC_charges)                # sigma [e-]
    RN_noise = RN_charges                           # sigma [e-] (deja un sigma, pas Poisson)

    # Total noise per exposure per pixel (quadrature sum)
    total_noise_per_exp = np.sqrt(signal_noise**2 + dark_noise**2 + CIC_noise**2 + RN_noise**2)

    # SNR per exposure per pixel
    snr_per_exp = signal_charges / total_noise_per_exp

    # ============================================================================
    # SNR TOTAL = SNR_per_exp x sqrt(n_exposures) x sqrt(n_pixels)
    # ============================================================================
    # Note: This is the simplified formula. More accurate would be:
    # Signal_total = signal_charges x n_exposures x n_pixels
    # Noise_total = sqrt(signal_noise^2 + dark_noise^2 + CIC_noise^2 + RN_noise^2) x sqrt(n_exposures) x sqrt(n_pixels)
    # But for display purposes, we use the scaling factor

    factor = np.sqrt(n_exposures) * np.sqrt(n_pixels)
    signal_total = signal_charges * n_exposures * n_pixels
    noise_total = total_noise_per_exp * np.sqrt(n_exposures) * np.sqrt(n_pixels)
    snr_total = signal_total / noise_total

    # Package results
    charges_per_exp = {
        'signal': signal_charges,
        'dark': dark_charges,
        'CIC': CIC_charges,
        'RN': RN_charges
    }

    noises_per_exp = {
        'signal': signal_noise,
        'dark': dark_noise,
        'CIC': CIC_noise,
        'RN': RN_noise,
        'total': total_noise_per_exp
    }

    return exptime, charges_per_exp, noises_per_exp, snr_per_exp, n_exposures, snr_total, signal_total, noise_total


def plot_snr_analysis(CIC=0.001, RN=0.01, dark_rate=0.1, signal_rate=0.02,
                     max_exptime=100, acquisition_time=1.0, n_pixels=1,
                     current_exptime=50, readout_time=0):
    """
    Interactive plotting function for SNR analysis.
    """
    # Calculate
    results = calculate_snr_vs_exptime(
        CIC=CIC, RN=RN, dark_rate=dark_rate, signal_rate=signal_rate,
        max_exptime=max_exptime, acquisition_time=acquisition_time,
        n_pixels=n_pixels, current_exptime=current_exptime, readout_time=readout_time
    )
    exptime, charges, noises, snr_per_exp, n_exposures, snr_total, signal_total, noise_total = results

    # Create figure
    fig, axes = plt.subplots(3, 1, figsize=(12, 11), sharex=True)

    # ============================================================================
    # SUBPLOT 1: CHARGES per exposure per pixel (nombre d'electrons)
    # ============================================================================
    ax1 = axes[0]
    ax1.plot(exptime, charges['signal'], 'b-', linewidth=2, label='Signal')
    ax1.plot(exptime, charges['dark'], 'r-', linewidth=2, label='Dark')
    ax1.plot(exptime, charges['CIC'], 'g--', linewidth=2, label='CIC')
    ax1.plot(exptime, charges['RN'], 'm--', linewidth=2, label='RN')

    # Vertical line at current exposure time
    ax1.axvline(x=current_exptime, color='black', linestyle='-', linewidth=2, alpha=0.5)

    ax1.set_ylabel('Charges [e-/pix/exp]', fontsize=12, fontweight='bold')
    ax1.set_title('Charges par exposition par pixel', fontsize=14, fontweight='bold')
    ax1.legend(loc='best', fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    ax1.set_ylim(1e-5, max(charges['signal'].max(), charges['dark'].max()) * 2)

    # ============================================================================
    # SUBPLOT 2: NOISES per exposure per pixel (sigma, bruit)
    # ============================================================================
    ax2 = axes[1]
    ax2.plot(exptime, noises['signal'], 'b-', linewidth=2, label='Signal noise (sqrt(N))')
    ax2.plot(exptime, noises['dark'], 'r-', linewidth=2, label='Dark noise (sqrt(N))')
    ax2.plot(exptime, noises['CIC'], 'g--', linewidth=2, label='CIC noise (sqrt(N))')
    ax2.plot(exptime, noises['RN'], 'm--', linewidth=2, label='RN (constant)')
    ax2.plot(exptime, noises['total'], 'k-', linewidth=3, label='Total noise', alpha=0.7)

    # Vertical line at current exposure time
    ax2.axvline(x=current_exptime, color='black', linestyle='-', linewidth=2, alpha=0.5)

    ax2.set_ylabel('Noise sigma [e-/pix/exp]', fontsize=12, fontweight='bold')
    ax2.set_title('Bruit par exposition par pixel (sigma)', fontsize=14, fontweight='bold')
    ax2.legend(loc='best', fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.set_yscale('log')
    ax2.set_ylim(1e-4, noises['total'].max() * 2)

    # ============================================================================
    # SUBPLOT 3: SNR TOTAL (with n_exposures and n_pixels)
    # ============================================================================
    ax3 = axes[2]

    # Plot SNR per exposure (1 pixel, 1 exposure)
    ax3.plot(exptime, snr_per_exp, 'gray', linewidth=2, linestyle='--',
             label='SNR (1 pix, 1 exp)', alpha=0.5)

    # Plot total SNR
    ax3.plot(exptime, snr_total, 'k-', linewidth=3,
             label='SNR total ({} pix, {:.1f}h)'.format(n_pixels, acquisition_time))

    # Vertical line at current exposure time
    ax3.axvline(x=current_exptime, color='black', linestyle='-', linewidth=2, alpha=0.5,
                label='Current: {:.1f}s'.format(current_exptime))

    # Horizontal line at SNR=5
    ax3.axhline(y=5, color='orange', linestyle='--', linewidth=2, alpha=0.7, label='SNR = 5')

    # Mark current SNR value
    idx_current = np.argmin(np.abs(exptime - current_exptime))
    current_snr = snr_total[idx_current]
    current_n_exp = n_exposures[idx_current]
    ax3.plot(current_exptime, current_snr, 'ro', markersize=12, zorder=10)
    ax3.text(current_exptime, current_snr * 1.1,
             'SNR={:.1f}\n{:.0f} exp'.format(current_snr, current_n_exp),
             fontsize=10, ha='center', fontweight='bold',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax3.set_xlabel('Exposure Time [seconds]', fontsize=12, fontweight='bold')
    ax3.set_ylabel('SNR', fontsize=12, fontweight='bold')
    ax3.set_title('Signal-to-Noise Ratio (Total: {:.1f}h, {} pixels)'.format(acquisition_time, n_pixels),
                  fontsize=14, fontweight='bold')
    ax3.legend(loc='best', fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim(0, max_exptime)
    ax3.set_ylim(0, max(snr_total.max() * 1.1, 10))

    plt.tight_layout()
    plt.show()

    # Print summary at current exposure time
    idx = idx_current
    print("=" * 80)
    print("SUMMARY at exposure time = {:.1f}s".format(current_exptime))
    print("=" * 80)
    print("Acquisition time:  {:.1f} hours = {:.0f} seconds".format(acquisition_time, acquisition_time*3600))
    print("Number of pixels:  {}".format(n_pixels))
    print("Number of exposures: {:.0f}".format(n_exposures[idx]))
    print("Readout time:      {:.1f} seconds".format(readout_time))
    print("-" * 80)
    print("PER EXPOSURE PER PIXEL:")
    print("  Signal:        {:.4f} e-  (noise: {:.4f} e-)".format(charges['signal'][idx], noises['signal'][idx]))
    print("  Dark:          {:.4f} e-  (noise: {:.4f} e-)".format(charges['dark'][idx], noises['dark'][idx]))
    print("  CIC:           {:.4f} e-  (noise: {:.4f} e-)".format(charges['CIC'][idx], noises['CIC'][idx]))
    print("  RN:            {:.4f} e-  (noise: {:.4f} e-)".format(charges['RN'][idx], noises['RN'][idx]))
    print("  Total noise:   {:.4f} e-".format(noises['total'][idx]))
    print("  SNR (1 pix, 1 exp): {:.2f}".format(snr_per_exp[idx]))
    print("-" * 80)
    print("TOTAL (all exposures, all pixels):")
    print("  Signal total:  {:.2f} e-".format(signal_total[idx]))
    print("  Noise total:   {:.2f} e-".format(noise_total[idx]))
    print("  SNR total:     {:.2f}".format(snr_total[idx]))
    print("  Scaling factor: sqrt({:.0f} exp) x sqrt({} pix) = {:.1f}".format(
        n_exposures[idx], n_pixels, np.sqrt(n_exposures[idx] * n_pixels)))
    print("=" * 80)


# ============================================================================
# Create interactive widget
# ============================================================================
if __name__ == "__main__":
    # Define widget sliders
    style = {'description_width': '180px'}
    layout = Layout(width='550px')

    CIC_slider = FloatSlider(
        value=0.001, min=0, max=0.01, step=0.0001,
        description='CIC [e-/pix/frame]',
        style=style, layout=layout,
        readout_format='.4f',
        continuous_update=False
    )

    RN_slider = FloatSlider(
        value=0.01, min=0, max=1, step=0.01,
        description='RN [e-/pix/frame]',
        style=style, layout=layout,
        readout_format='.3f',
        continuous_update=False
    )

    dark_slider = FloatSlider(
        value=0.1, min=0, max=10, step=0.1,
        description='Dark [e-/pix/hour]',
        style=style, layout=layout,
        readout_format='.2f',
        continuous_update=False
    )

    signal_slider = FloatSlider(
        value=0.02, min=0, max=1, step=0.01,
        description='Signal [e-/pix/hour]',
        style=style, layout=layout,
        readout_format='.3f',
        continuous_update=False
    )

    max_exptime_slider = FloatSlider(
        value=100, min=10, max=1000, step=10,
        description='Max Exptime [s] (plot)',
        style=style, layout=layout,
        readout_format='.0f',
        continuous_update=False
    )

    current_exptime_slider = FloatSlider(
        value=50, min=0.1, max=1000, step=1,
        description='Current Exptime [s]',
        style=style, layout=layout,
        readout_format='.1f',
        continuous_update=True  # Real-time update for vertical bar
    )

    acquisition_time_slider = FloatSlider(
        value=1.0, min=0.1, max=100, step=0.1,
        description='Total Acq. Time [hours]',
        style=style, layout=layout,
        readout_format='.1f',
        continuous_update=False
    )

    n_pixels_slider = IntSlider(
        value=1, min=1, max=100, step=1,
        description='Number of pixels',
        style=style, layout=layout,
        continuous_update=False
    )

    readout_time_slider = FloatSlider(
        value=0, min=0, max=60, step=1,
        description='Readout time [s]',
        style=style, layout=layout,
        readout_format='.0f',
        continuous_update=False
    )

    # Create interactive plot
    interact(plot_snr_analysis,
             CIC=CIC_slider,
             RN=RN_slider,
             dark_rate=dark_slider,
             signal_rate=signal_slider,
             max_exptime=max_exptime_slider,
             current_exptime=current_exptime_slider,
             acquisition_time=acquisition_time_slider,
             n_pixels=n_pixels_slider,
             readout_time=readout_time_slider)
