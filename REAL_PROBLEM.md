# Cross-validation Generic ETC vs Astropy

## Test simple : comparaison directe des SNR

```bash
python test_sky_renormalization.py
```

### Principe

1. **Calculer avec Generic ETC** → obtenir Signal_el, sky, Dark_current_f (en e⁻)
2. **Convertir en taux** → e⁻/s pour Astropy
3. **Injecter dans Astropy** → signal_to_noise_oir_ccd()
4. **Comparer les SNR** → Generic vs Astropy
5. **Tester les variations** → exposure_time, sky, signal, dark

### Ce qu'on vérifie

**Si le ratio SNR_Generic/SNR_Astropy est constant** → les variations sont cohérentes

**Si le ratio ≈ 1.0** → les calculs matchent parfaitement

**Si le ratio ≠ 1.0 mais constant** → différence systématique (convention, unités, etc.) mais calcul correct

### Code du test

```python
# Generic ETC
obs = Observation(instruments, instrument="GALEX FUV", exposure_time=1000,
                  SNR_res="per pix", IFS=False, test=True)

# Convertir en taux (e⁻ → e⁻/s)
source_eps = obs.Signal_el[obs.i] / obs.exposure_time
sky_eps = obs.sky[obs.i] / obs.exposure_time
dark_eps = obs.Dark_current_f[obs.i] / obs.exposure_time

# Astropy
snr_astropy = signal_to_noise_oir_ccd(
    t=obs.exposure_time,
    source_eps=source_eps,
    sky_eps=sky_eps,
    dark_eps=dark_eps,
    rd=RN,
    npix=obs.number_pixels_used,  # 1 en mode "per pix"
    gain=1.0
)

# Comparer
ratio = obs.SNR[obs.i] / snr_astropy
```

## Résultats attendus

Si tout est correct, on devrait voir :
- Ratio constant pour toutes les variations
- Évolutions identiques (si sky ×2 → SNR baisse de même façon)

Si problème :
- Ratio qui varie → erreur dans le calcul
- Évolutions différentes → problème physique (mauvaise formule)

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
