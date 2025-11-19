# Solution finale pour la cross-validation

## Problème identifié

Le Generic ETC calcule :
- `Signal_el` = signal **moyen** sur `pixels_total_source` pixels
- `sky` = sky **moyen** sur `pixels_total_source` pixels, MAIS avec une aire spatiale différente

Astropy attend :
- Signal et sky **par pixel individuel**

## Différences clés

### 1. Division par `pixels_total_source`

```python
# Generic ETC (ligne 2481, 2504)
factor_CU2el = ... / pixels_total_source
```

Cela moyenne sur plusieurs pixels. Pour obtenir la valeur par pixel :
```python
Signal_el_per_pixel = Signal_el × pixels_total_source
```

### 2. Aires spatiales différentes

```python
# Signal (ligne 2481)
factor_CU2el_tot ∝ source_size_arcsec_after_slit

# Sky (ligne 2483)
factor_CU2el_sky_tot ∝ slit_size_arcsec_after_slit
```

Si `slit_size > source_size`, alors `factor_sky > factor_signal`.

### 3. Couvertures spectrales différentes

```python
# Signal
spectral_coverage_signal = min(Line_width, Bandwidth)

# Sky
spectral_coverage_sky = max(min(Line_width, Bandwidth), dispersion)
```

Le sky intègre au minimum sur 1 pixel spectral (`dispersion`).

## Solution : Ajouter des variables "per pixel" pour cross-validation

Insérer ce code dans `Observation.py` après la ligne 2545 :

```python
#####################################
# FOR CROSS-VALIDATION WITH ASTROPY
# Compute per-pixel values (not averaged over pixels_total_source)
#####################################
if self.SNR_res == "per pix":
    # Factor per individual pixel (not divided by pixels_total_source)
    # This matches what Astropy expects
    self.factor_CU2el_astropy = (self.effective_area * self.arcsec2str *
                                  self.pixel_scale**2 * self.dispersion)

    # Use same factor for sky (same pixel area)
    self.factor_CU2el_sky_astropy = self.factor_CU2el_astropy

    # Signal electrons per pixel
    self.Signal_el_astropy = (self.Signal_LU * self.factor_CU2el_astropy *
                              self.exposure_time * self.flux_fraction_slit_applied)

    # Sky electrons per pixel
    self.sky_astropy = (self.Sky_CU * self.factor_CU2el_sky_astropy *
                        self.exposure_time)

    # For reference: show the correction factor
    self.astropy_correction_factor = self.pixels_total_source
```

## Utilisation pour cross-validation

Après cette modification, comparer :

```python
from astropy.stats import signal_to_noise_oir_ccd

# Generic ETC (with new Astropy-compatible variables)
obs = Observation(instruments=instruments, instrument="GALEX FUV",
                  exposure_time=1000, SNR_res="per pix", IFS=False, test=True)

# Astropy
wavelength_nm = params['wavelength']
E_photon = 1.986e-8 / wavelength_nm

signal_photons = params['Signal'] / E_photon
sky_photons = params['Sky'] / E_photon

collecting_area_cm2 = params['Collecting_area'] * 1e4
pixel_area_arcsec2 = params['pixel_scale']**2
wavelength_range = params['dispersion']
total_throughput = params['Throughput'] * params['QE'] * params['Atmosphere']

source_eps = (signal_photons * collecting_area_cm2 * pixel_area_arcsec2 *
              wavelength_range * total_throughput)
sky_eps = (sky_photons * collecting_area_cm2 * pixel_area_arcsec2 *
           wavelength_range * total_throughput)
dark_eps = params['Dark_current'] / 3600.0

snr_astropy = signal_to_noise_oir_ccd(
    t=1000, source_eps=source_eps, sky_eps=sky_eps,
    dark_eps=dark_eps, rd=params['RN'], npix=1, gain=1.0
)

# Compare
print(f"Signal (Generic): {obs.Signal_el_astropy:.6e} e⁻")
print(f"Signal (Astropy): {source_eps * 1000:.6e} e⁻")
print(f"Ratio: {obs.Signal_el_astropy / (source_eps * 1000):.4f}")  # Should be ~1.0

print(f"Sky (Generic): {obs.sky_astropy:.6e} e⁻")
print(f"Sky (Astropy): {sky_eps * 1000:.6e} e⁻")
print(f"Ratio: {obs.sky_astropy / (sky_eps * 1000):.4f}")  # Should be ~1.0
```

## Alternative : Correction temporaire sans modifier Observation.py

Si vous ne voulez pas modifier Observation.py, utilisez dans votre notebook :

```python
# Correction factor
correction = obs.pixels_total_source

# Corrected values
signal_corrected = obs.Signal_el * correction
sky_corrected = obs.sky * correction

# MAIS ATTENTION: sky_corrected peut encore être faux si
# slit_size_arcsec_after_slit != source_size_arcsec_after_slit

# Pour le sky, recalculer avec la bonne aire:
if hasattr(obs, 'slit_size_arcsec_after_slit') and hasattr(obs, 'source_size_arcsec_after_slit'):
    sky_area_correction = obs.source_size_arcsec_after_slit / obs.slit_size_arcsec_after_slit
    sky_corrected_final = sky_corrected * sky_area_correction
else:
    # Fallback: recalculate from scratch
    sky_corrected_final = (obs.Sky_CU * obs.effective_area * obs.arcsec2str *
                          obs.pixel_scale**2 * obs.dispersion * obs.exposure_time)
```

## Test rapide

Pour vérifier que la correction fonctionne :

```bash
python analyze_factor_formulas.py
```

Cela affichera les 3 formules et identifiera laquelle est la plus proche d'Astropy.

## Commit

Après vérification, commiter les fichiers d'analyse :

```bash
git add DETAILED_ANALYSIS.md FINAL_SOLUTION.md analyze_factor_formulas.py run_comprehensive_tests.py
git commit -m "Add comprehensive analysis and solution for cross-validation

- Identified 3 overlapping formulas in Observation.py (lines 2481-2527)
- Formula 1 (_tot) and Formula 3 (_average) divide by pixels_total_source
- Formula 2 (line 2488-2492) is overwritten and never used
- Signal uses source_size, Sky uses slit_size (different areas!)

Solution: Add *_astropy variables that compute per-pixel values
matching Astropy's expectations (pixel_scale² × dispersion)
"
git push
```
