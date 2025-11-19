# Analyse du problème de cross-validation ETC

## Problème observé

D'après les résultats du notebook de debug pour GALEX FUV (t_exp=1000s):

| Paramètre | Generic ETC | Astropy | Ratio (Gen/Ast) |
|-----------|-------------|---------|-----------------|
| Signal (e⁻) | 0.000061 | 0.000912 | ~0.07 (15x moins) |
| Sky (e⁻) | 0.140896 | 0.000249 | ~565 (565x plus) |
| Dark (e⁻) | 0.013611 | 0.013611 | 1.00 ✓ |
| SNR | 0.002570 | 0.007504 | ~0.34 |

## Cause identifiée

### Dans Observation.py (lignes 2481-2483)

```python
self.factor_CU2el_tot = 1*self.effective_area * self.arcsec2str * \
    np.minimum(self.Line_width,self.Bandwidth) * \
    self.source_size_arcsec_after_slit / self.pixels_total_source

self.factor_CU2el_sky_tot = 1*self.effective_area * self.arcsec2str * \
    np.maximum(np.minimum(self.Line_width,self.Bandwidth),self.dispersion) * \
    self.slit_size_arcsec_after_slit / self.pixels_total_source
```

**Les deux facteurs divisent par `pixels_total_source`** !

Cela signifie que :
- `Signal_el` = signal **moyen** par pixel (moyenné sur `pixels_total_source` pixels)
- `sky` = sky **moyen** par pixel (moyenné sur `pixels_total_source` pixels)

### Comparaison avec Astropy

Astropy attend des taux d'électrons **par pixel individuel** :
- `source_eps` : e⁻/s/pix (pour un pixel spécifique)
- `sky_eps` : e⁻/s/pix (pour un pixel spécifique)

## Analyse détaillée

### Pourquoi Signal est trop PETIT dans Generic ETC ?

`Signal_el` utilise `factor_CU2el` qui est divisé par `pixels_total_source`.

Pour GALEX FUV, si `pixels_total_source` ≈ 15, alors:
- Signal Generic = Signal correct / 15
- Ratio observé ≈ 0.07 ✓ (correspond!)

### Pourquoi Sky est trop GRAND dans Generic ETC ?

Attendez... sky utilise `factor_CU2el_sky` qui est **aussi** divisé par `pixels_total_source`.
Donc sky devrait aussi être trop petit, pas trop grand !

Il y a donc un **deuxième problème** avec le calcul du sky.

Regardons la ligne 2539 :
```python
self.sky = self.Sky_CU * self.factor_CU2el_sky * self.exposure_time
```

Et ligne 2540 :
```python
self.Sky_noise = np.sqrt(self.sky * self.ENF)
```

**HYPOTHÈSE** : Le notebook de debug affiche peut-être `obs.Sky_noise` au lieu de `obs.sky` ?
- Si c'est le cas, on comparerait un **bruit** (√electrons) avec des **électrons** !

Vérification :
- Si sky electrons ≈ 0.000249, alors Sky_noise = √0.000249 ≈ 0.0158
- Mais on observe 0.140896, ce n'est donc pas ça.

**AUTRE HYPOTHÈSE** : Problème dans `factor_CU2el_sky`

Le calcul à la ligne 2483 utilise :
- `slit_size_arcsec_after_slit` pour le sky
- `source_size_arcsec_after_slit` pour le signal

Si `slit_size >> source_size`, alors factor_CU2el_sky >> factor_CU2el.

## Solutions possibles

### Option 1 : Modifier les paramètres Astropy (temporaire)

Pour la comparaison, on peut multiplier `source_eps` par `pixels_total_source` dans Astropy :

```python
# Au lieu de comparer par pixel
source_eps_astropy = ... # e⁻/s/pix

# Comparer sur pixels_total_source pixels
source_eps_for_comparison = source_eps_astropy / obs.pixels_total_source
```

### Option 2 : Corriger les valeurs Generic ETC comparées

```python
# Multiplier les valeurs Generic ETC pour obtenir le signal sur 1 pixel
Signal_el_per_single_pixel = obs.Signal_el * obs.pixels_total_source
sky_per_single_pixel = obs.sky * obs.pixels_total_source
```

### Option 3 : Utiliser les bons attributs pour la comparaison

Il faut vérifier qu'on compare bien :
- `obs.sky` (électrons) et pas `obs.Sky_noise` (bruit)
- Les valeurs **avant** ou **après** multiplication par `factor^2`

## Investigation requise

1. **Imprimer les valeurs de debug** :
   - `obs.pixels_total_source`
   - `obs.source_size_arcsec_after_slit`
   - `obs.slit_size_arcsec_after_slit`
   - `obs.factor_CU2el` vs `obs.factor_CU2el_sky`

2. **Vérifier ce que le notebook de debug affiche** :
   - Est-ce `obs.sky` ou `obs.Sky_noise` ?
   - Est-ce `obs.Signal_el` ou autre chose ?

3. **Tester avec test=True vs test=False** :
   - `test=True` utilise `factor_CU2el_tot` (divise par pixels_total_source)
   - `test=False` utilise `factor_CU2el_average` (divise aussi par pixels_total_source)

## Recommandation

**Pour la cross-validation**, créer des versions "per pixel" des facteurs :

```python
# Version "per pixel" sans division par pixels_total_source
factor_CU2el_per_pix = effective_area * arcsec2str * pixel_scale^2 * dispersion

# Pour le sky
factor_CU2el_sky_per_pix = effective_area * arcsec2str * pixel_scale^2 * dispersion
```

Puis utiliser ces facteurs pour calculer les taux e⁻/s/pix à comparer avec Astropy.
