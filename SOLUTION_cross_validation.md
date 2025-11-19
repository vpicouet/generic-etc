# Solution au problème de cross-validation

## Problème identifié

Dans `Observation.py` lignes **2481-2483** et **2523-2524** :

```python
self.factor_CU2el_tot = ... / self.pixels_total_source
self.factor_CU2el_sky_tot = ... / self.pixels_total_source
```

Ces facteurs **divisent par `pixels_total_source`**, ce qui fait que :
- `Signal_el` = signal **moyen** sur N pixels
- `sky` = sky **moyen** sur N pixels

Mais **Astropy attend des valeurs par pixel individuel**.

## Résultats observés expliqués

Pour GALEX FUV avec `pixels_total_source ≈ 15` :

| Variable | Observé | Explication |
|----------|---------|-------------|
| Signal Generic/Astropy | 0.07 (15x moins) | Divisé par `pixels_total_source` ≈ 15 ✓ |
| Sky Generic/Astropy | 565x plus | **Problème additionnel** avec `factor_CU2el_sky` |
| Dark | 1.00 ✓ | OK (pas affecté par factor_CU2el) |

## Diagnostic du problème Sky

Le sky utilise `factor_CU2el_sky` qui multiplie par `slit_size_arcsec_after_slit` au lieu de `pixel_scale²`.

Si la fente est grande, `slit_size >> pixel_area`, donc sky est surestimé.

## Solutions

### Solution 1 : Pour la cross-validation (TEMPORAIRE)

Modifier les valeurs comparées côté **Astropy** :

```python
# Au lieu de npix=1
snr_astropy = signal_to_noise_oir_ccd(
    t=t_exp,
    source_eps=source_eps / obs.pixels_total_source,  # Diviser par pixels_total_source
    sky_eps=sky_eps / obs.pixels_total_source,        # Diviser par pixels_total_source
    dark_eps=dark_eps,
    rd=RN,
    npix=obs.pixels_total_source,  # Utiliser N pixels
    gain=1.0
)
```

**OU** multiplier les valeurs Generic ETC :

```python
signal_generic_corrected = obs.Signal_el * obs.pixels_total_source
sky_generic_corrected = obs.sky * obs.pixels_total_source
```

### Solution 2 : Créer des versions "per pix" des facteurs

Pour la cross-validation propre, ajouter dans `Observation.py` :

```python
# Après la ligne 2527, ajouter :
if SNR_res == "per pix":
    # Pour la cross-validation, on veut les valeurs par pixel individuel
    self.factor_CU2el_per_pix = (self.effective_area * self.arcsec2str *
                                  pixel_area_arcsec2 * wavelength_range_angstrom)
    self.factor_CU2el_sky_per_pix = self.factor_CU2el_per_pix

    # Variables pour cross-validation
    self.Signal_el_per_pix = self.Signal_LU * self.factor_CU2el_per_pix * self.exposure_time
    self.sky_per_pix = self.Sky_CU * self.factor_CU2el_sky_per_pix * self.exposure_time
```

Puis comparer `obs.Signal_el_per_pix` et `obs.sky_per_pix` avec Astropy.

## Test rapide

Exécuter le script de debug :

```bash
python debug_cross_validation.py
```

Cela affichera toutes les valeurs intermédiaires et confirmera le diagnostic.

## Recommandation

**Pour l'instant** (sans modifier l'ETC) :
1. Multiplier `obs.Signal_el` par `obs.pixels_total_source` avant de comparer
2. Multiplier `obs.sky` par `obs.pixels_total_source` avant de comparer
3. Vérifier que les ratios sont maintenant ~1.0

**Pour le sky**, il y a probablement un problème additionnel avec `slit_size_arcsec_after_slit`.
À investiguer en regardant les valeurs de :
- `obs.slit_size_arcsec_after_slit`
- `obs.source_size_arcsec_after_slit`
- `params['pixel_scale']²`
