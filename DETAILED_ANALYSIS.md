# Analyse détaillée du problème de cross-validation

## Structure du code (Observation.py lignes 2481-2527)

Le code a **3 formules différentes** qui s'écrasent mutuellement :

### Formula 1: `factor_CU2el_tot` (lignes 2481-2485)

```python
# Line 2481
self.factor_CU2el_tot = (effective_area * arcsec2str *
                         min(Line_width, Bandwidth) *
                         source_size_arcsec_after_slit /
                         pixels_total_source)

# Line 2483
self.factor_CU2el_sky_tot = (effective_area * arcsec2str *
                             max(min(Line_width, Bandwidth), dispersion) *
                             slit_size_arcsec_after_slit /
                             pixels_total_source)

# Line 2484-2485
self.factor_CU2el = self.factor_CU2el_tot
self.factor_CU2el_sky = self.factor_CU2el_sky_tot
```

**Caractéristiques:**
- Divise par `pixels_total_source`
- Utilise `source_size_arcsec_after_slit` pour signal
- Utilise `slit_size_arcsec_after_slit` pour sky (DIFFÉRENT !)
- Donne des valeurs **moyennées sur plusieurs pixels**

### Formula 2: Lignes 2488-2492 (ÉCRASE Formula 1)

```python
# Line 2488 - OVERWRITES!
self.factor_CU2el = (effective_area * arcsec2str *
                     min(Slitwidth, Size_source) *
                     min(dispersion, min(Line_width, Bandwidth)/2.35) *
                     pixel_scale)

# Line 2491-2492 - OVERWRITES!
sky_spectral_coverage = min(dispersion, min(Line_width, Bandwidth)/2.35)
self.factor_CU2el_sky = (effective_area * arcsec2str *
                         Slitwidth *
                         sky_spectral_coverage *
                         pixel_scale)
```

**Caractéristiques:**
- Ne divise PAS par `pixels_total_source`
- Multiplie par `pixel_scale` (arcsec/pix)
- Utilise `Slitwidth * pixel_scale` au lieu de `pixel_scale²`
- **Cette formule est écrasée par Formula 3 ci-dessous !**

### Formula 3: `factor_CU2el_average` (lignes 2504-2527)

```python
# Line 2497-2504
source_spatial_arcsec = min(Slitwidth, Size_source)
source_spectral_angstrom = min(Line_width, Bandwidth)
self.factor_CU2el_average = (effective_area * arcsec2str *
                             source_spatial_arcsec *
                             source_spectral_angstrom /
                             pixels_total_source)

# Line 2506-2514
sky_spatial_arcsec = Slitwidth  # Could also be slit_size_arcsec_after_slit
sky_spectral_angstrom = min(Line_width, Bandwidth)
self.factor_CU2el_sky_average = (effective_area * arcsec2str *
                                 sky_spatial_arcsec *
                                 sky_spectral_angstrom /
                                 pixels_total_source)

# Line 2522-2527 - FINAL OVERWRITE!
if self.test:
    self.factor_CU2el = self.factor_CU2el_tot        # Use Formula 1
    self.factor_CU2el_sky = self.factor_CU2el_sky_tot
else:
    self.factor_CU2el = self.factor_CU2el_average    # Use Formula 3
    self.factor_CU2el_sky = self.factor_CU2el_sky_average
```

**Caractéristiques:**
- Divise par `pixels_total_source` (comme Formula 1)
- Formule plus simple que Formula 1
- **C'est cette formule qui est utilisée au final** (sauf Formula 2 qui est calculée entre-temps puis écrasée)

## Quelle formule est réellement utilisée ?

D'après le flux du code :

1. **Ligne 2484-2485** : Assigne Formula 1 → `factor_CU2el = factor_CU2el_tot`
2. **Ligne 2488-2492** : **ÉCRASE** avec Formula 2
3. **Ligne 2504-2514** : Calcule Formula 3 (sans assigner)
4. **Ligne 2522-2527** : **ÉCRASE À NOUVEAU** :
   - Si `test=True` → utilise Formula 1
   - Si `test=False` → utilise Formula 3

**Donc Formula 2 (lignes 2488-2492) n'est JAMAIS utilisée car elle est écrasée !**

## Comparaison avec Astropy

### Ce qu'Astropy attend (par pixel individuel) :

```python
factor_per_pixel = (collecting_area_cm2 * throughput_total *
                   (π/180/3600)² *
                   pixel_scale² *
                   dispersion)
```

**Unités** :
- `collecting_area_cm2` : cm²
- `throughput_total` : sans dimension
- `(π/180/3600)²` : conversion arcsec² → str²
- `pixel_scale²` : arcsec²/pix²
- `dispersion` : Å/pix

**Résultat** : e⁻ par (LU × seconde × pixel)

### Différences avec Generic ETC :

| Formule | Spatial | Spectral | Division | Match Astropy? |
|---------|---------|----------|----------|----------------|
| Astropy | `pixel_scale²` | `dispersion` | - | ✓ (référence) |
| Formula 1 (_tot) | `source_size_arcsec²` | `Line_width` | `÷ pixels_total_source` | ✗ × N |
| Formula 2 (ligne 2488) | `Slitwidth × pixel_scale` | `dispersion` | - | ✗ (mauvaise unité) |
| Formula 3 (_average) | `Slitwidth × Size_source` | `Line_width` | `÷ pixels_total_source` | ✗ × N |

**Aucune formule ne matche Astropy !**

## Problèmes identifiés

### Problème 1 : Division par `pixels_total_source`

Formulas 1 et 3 divisent par `pixels_total_source`, ce qui donne un signal **moyen** sur N pixels au lieu du signal **par pixel**.

**Impact** : Facteur ~15x pour GALEX FUV

### Problème 2 : Aire spatiale différente pour signal vs sky

- **Signal** : utilise `source_size_arcsec_after_slit` (limité par source)
- **Sky** : utilise `slit_size_arcsec_after_slit` (taille de la fente)

Si `slit >> source`, alors `factor_CU2el_sky >> factor_CU2el`

**Exemple GALEX FUV** :
- Source : 4 arcsec × 4 arcsec = 16 arcsec²
- Slit : 4 arcsec × 20 arcsec = 80 arcsec²
- Ratio : 80/16 = 5x

Combiné avec la division par `pixels_total_source` qui varie, cela crée des facteurs différents.

### Problème 3 : Couverture spectrale différente

- **Signal** : `min(Line_width, Bandwidth)` (largeur de la raie)
- **Sky** : `max(min(Line_width, Bandwidth), dispersion)` (au moins 1 pixel spectral)

Pour une raie fine (Line_width < dispersion), le sky intègre sur plus de Å que le signal.

### Problème 4 : Formula 2 a des unités incorrectes

Formula 2 (ligne 2488) multiplie :
```python
Slitwidth [arcsec] × pixel_scale [arcsec/pix] = [arcsec²/pix]
```

Mais Astropy attend :
```python
pixel_scale² [arcsec²/pix²] × 1 [pix] = [arcsec²/pix]
```

Les unités sont les mêmes mais la valeur est différente si `Slitwidth ≠ pixel_scale`.

## Solution pour la cross-validation

Pour comparer Generic ETC avec Astropy, il faut :

### Option 1 : Corriger Generic ETC pour obtenir valeurs par pixel

```python
# Multiplier par pixels_total_source
signal_per_pixel = obs.Signal_el * obs.pixels_total_source
sky_per_pixel = obs.sky * obs.pixels_total_source

# Puis comparer avec Astropy (npix=1)
```

**Problème** : Cela ne fonctionne que si le facteur `pixels_total_source` est constant, ce qui n'est pas le cas si la source ou les paramètres changent.

### Option 2 : Créer des facteurs "per pixel" dans Generic ETC

Ajouter dans Observation.py après ligne 2527 :

```python
# For cross-validation with Astropy
if SNR_res == "per pix":
    self.factor_CU2el_per_pix = (self.effective_area * self.arcsec2str *
                                  self.pixel_scale**2 * self.dispersion)
    self.factor_CU2el_sky_per_pix = self.factor_CU2el_per_pix

    self.Signal_el_per_pix = (self.Signal_LU * self.factor_CU2el_per_pix *
                              self.exposure_time * self.flux_fraction_slit_applied)
    self.sky_per_pix = (self.Sky_CU * self.factor_CU2el_sky_per_pix *
                        self.exposure_time)
```

Puis comparer `obs.Signal_el_per_pix` et `obs.sky_per_pix` avec Astropy.

### Option 3 : Modifier les paramètres Astropy

Au lieu de `npix=1`, utiliser :

```python
snr_astropy = signal_to_noise_oir_ccd(
    t=t_exp,
    source_eps=source_eps / obs.pixels_total_source,
    sky_eps=sky_eps / obs.pixels_total_source,
    dark_eps=dark_eps,
    rd=RN,
    npix=obs.pixels_total_source,
    gain=1.0
)
```

**Problème** : Le ratio sky/signal dans Generic ETC est faussé par les aires différentes (`source_size` vs `slit_size`).

## Recommandation

**Pour fixer la cross-validation immédiatement** :

1. Utiliser **Option 2** : Créer `Signal_el_per_pix` et `sky_per_pix`
2. Corriger le calcul du sky pour utiliser la même aire que le signal

**Pour fixer Generic ETC à long terme** :

1. Nettoyer les formulas 1/2/3 qui se chevauchent
2. Décider quelle formule est la bonne et supprimer les autres
3. Documenter clairement si les valeurs sont "par pixel" ou "moyennées sur N pixels"
4. Utiliser la même aire spatiale pour signal et sky (ou documenter pourquoi c'est différent)
