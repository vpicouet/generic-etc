# Cross-validation Generic ETC vs Astropy

## Résumé du problème et solution

**Problème identifié** : Generic ETC calcule SNR pour des expositions **stackées** sur `acquisition_time`, tandis qu'Astropy calcule SNR pour une **exposition unique**.

**Solution** : Passer `acquisition_time=exposure_time/3600.0` lors de la création de l'objet `Observation` pour obtenir `N_images=1`.

## Test simple : comparaison directe des SNR

```bash
python test_sky_renormalization.py
```

### Principe

1. **Calculer avec Generic ETC** avec `acquisition_time=exposure_time/3600` → obtenir Signal_el, sky, Dark_current_f (en e⁻) pour une seule exposition
2. **Convertir en taux** → e⁻/s pour Astropy
3. **Injecter dans Astropy** → signal_to_noise_oir_ccd()
4. **Comparer les SNR** → Generic vs Astropy
5. **Tester les variations** → exposure_time, sky, signal, dark

### Résultats

**Ratio SNR_Generic/SNR_Astropy = 0.9993** → les calculs matchent à 0.07% près ✓

**Le ratio reste constant pour toutes les variations** → cohérence parfaite ✓

#### Variations testées

| Paramètre | Facteurs testés | Ratio | Variation |
|-----------|-----------------|-------|-----------|
| **Exposure time** | 100s à 5000s | 0.9993 | Constant |
| **Sky background** | 0.1x à 5.0x | 0.996 - 0.999 | Très faible |
| **Dark current** | 0.1x à 10x | 0.9993 | Constant |
| **Signal flux** | 0.5x à 10x | 0.999 - 0.999 | Très faible |

**Statistiques finales** :
- Ratio moyen : **0.9992**
- Écart-type : **0.0007** (excellent!)

### Code du test

```python
# Generic ETC - IMPORTANT: set acquisition_time = exposure_time/3600 for single exposure
t_exp = 1000.0  # seconds
obs = Observation(instruments, instrument="GALEX FUV",
                  exposure_time=t_exp,
                  acquisition_time=t_exp/3600.0,  # Convert to hours, ensures N_images=1
                  SNR_res="per pix", IFS=False, test=True)

# Convertir en taux (e⁻ → e⁻/s)
source_eps = obs.Signal_el[obs.i] / t_exp
sky_eps = obs.sky[obs.i] / t_exp
dark_eps = obs.Dark_current_f[obs.i] / t_exp

# Astropy
snr_astropy = signal_to_noise_oir_ccd(
    t=t_exp,
    source_eps=source_eps,
    sky_eps=sky_eps,
    dark_eps=dark_eps,
    rd=RN,
    npix=obs.number_pixels_used,  # 1 en mode "per pix"
    gain=1.0
)

# Comparer
ratio = obs.SNR[obs.i] / snr_astropy  # Should be ~0.999
```

## Explication technique

### Le facteur N_images

Generic ETC calcule :
```python
N_images = acquisition_time * 3600 / (exposure_time + readout_time)
N_images_true = N_images * (1 - cosmic_ray_loss)
factor = sqrt(number_pixels_used) × sqrt(N_resol_element_A) × sqrt(N_images_true)
SNR = Signal_el × factor / sqrt(noise²)
```

Si `acquisition_time >> exposure_time`, alors `N_images > 1` et le SNR est multiplié par `sqrt(N_images)` car on stack plusieurs expositions.

Pour comparer avec Astropy (qui calcule pour une seule exposition), il faut `N_images = 1`, donc `acquisition_time = exposure_time / 3600`.

## Remarques sur les différences observées initialement

Les facteurs ~15x pour signal et ~565x pour sky viennent de :

1. **Division par pixels_total_source** (~15x)
   - Generic ETC calcule signal/sky **moyens** sur N pixels
   - Astropy calcule **par pixel individuel**

2. **Aires différentes pour signal vs sky** (~35x pour sky)
   - Signal : intégré sur `source_size_arcsec_after_slit`
   - Sky : intégré sur `slit_size_arcsec_after_slit` (plus grand si fente > source)

3. **Couvertures spectrales différentes**
   - Signal : `min(Line_width, Bandwidth)`
   - Sky : `max(min(Line_width, Bandwidth), dispersion)` (au moins 1 pixel spectral)

Ces différences sont **physiques** pour un spectrographe :
- Le bruit de fond vient de toute la fente, pas juste du pixel de la source
- C'est correct pour Generic ETC (design pour spectrographes)
- Pour comparer avec Astropy, il faut injecter les valeurs calculées par Generic ETC

## Si le test échoue

1. **Ratio constant mais ≠ 1** → vérifier les unités, ENF, facteurs
2. **Ratio qui varie** → problème dans une des formules (signal, sky, ou SNR)
3. **Variations différentes** → vérifier les dépendances (sky², signal², etc.)
