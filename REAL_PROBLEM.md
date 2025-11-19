# Le VRAI problème : définition des paramètres d'entrée

## TL;DR

**Le calcul du SNR est probablement correct dans les deux cas.**

Le problème est que **Generic ETC et Astropy n'interprètent PAS les paramètres d'entrée de la même façon** !

## Définition de "Sky" dans les deux ETCs

### Astropy
```python
# Sky = surface brightness en erg/cm²/s/arcsec²/Å
# Intégration sur 1 pixel :
sky_electrons = Sky × pixel_scale² × dispersion × ...
#                     └─ arcsec²/pix
```

**Aire d'intégration** : `pixel_scale²` (ex: 1.5² = 2.25 arcsec²)

### Generic ETC
```python
# Sky = surface brightness en erg/cm²/s/arcsec²/Å (même unité!)
# Intégration sur la fente :
sky_electrons = Sky × slit_size_arcsec_after_slit × dispersion × ...
#                     └─ arcsec² (taille de la fente)
```

**Aire d'intégration** : `slit_size_arcsec_after_slit` (ex: 4 × 20 = 80 arcsec²)

## Le ratio qu'on observe

Pour GALEX FUV :
- `pixel_scale²` ≈ 2.25 arcsec²
- `slit_size_arcsec_after_slit` ≈ 80 arcsec²
- **Ratio** : 80 / 2.25 ≈ **35x**

Combiné avec `pixels_total_source` ≈ 15, ça explique le facteur ~565x observé !

```
565 ≈ 15 (pixels_total_source) × 35 (ratio aires)
```

## C'est quoi le "bon" comportement ?

**Les deux sont corrects !** C'est juste une question de convention :

### Convention Astropy (par pixel)
- **Input** : `Sky` = surface brightness
- **Sortie** : SNR **par pixel**
- **Usage** : Empiler plusieurs pixels pour améliorer le SNR

### Convention Generic ETC (par source)
- **Input** : `Sky` = surface brightness
- **Sortie** : SNR **pour toute la source** (sur N pixels)
- **Usage** : Voir directement le SNR de détection

## Pour faire une vraie cross-validation

Il faut comparer la **même chose** :

### Option 1 : Comparer "per pixel"

Dans Generic ETC, renormaliser les inputs :

```python
# Aires d'intégration
pixel_area = pixel_scale²
source_area = source_size_arcsec_after_slit
slit_area = slit_size_arcsec_after_slit

# Renormaliser
Sky_for_generic = Sky_astropy × (pixel_area / slit_area)
Signal_for_generic = Signal_astropy × (pixel_area / source_area)

# Créer observation
obs = Observation(..., Sky=Sky_for_generic, Signal=Signal_for_generic)

# Multiplier par pixels_total_source pour avoir "per pixel"
signal_per_pixel = obs.Signal_el × obs.pixels_total_source
sky_per_pixel = obs.sky × obs.pixels_total_source

# Comparer avec Astropy (npix=1)
```

### Option 2 : Comparer "per source"

Dans Astropy, intégrer sur N pixels :

```python
# Calculer combien de pixels couvre la source dans Generic ETC
npix = obs.pixels_total_source

# Ajuster les eps pour Astropy
source_eps_total = source_eps / npix
sky_eps_total = sky_eps × (slit_area / pixel_area) / npix

# SNR Astropy sur N pixels
snr_astropy = signal_to_noise_oir_ccd(
    ...,
    source_eps=source_eps_total,
    sky_eps=sky_eps_total,
    npix=npix
)
```

## Test simple

Si l'hypothèse est correcte, alors en renormalisant les inputs, tout devrait matcher :

```bash
python test_sky_renormalization.py
```

### Ce que fait le test

**PARTIE 1: Calcul manuel**
- Reproduit le calcul de Observation.py ligne par ligne
- Calcule `source_size_arcsec_after_slit` et `slit_size_arcsec_after_slit`
- Calcule `factor_CU2el_tot` et `factor_CU2el_sky_tot`
- Compare brut vs Astropy
- Compare avec correction `× pixels_total_source`
- **Teste la renormalisation** : `Signal_renorm = Signal / (area_ratio × spectral_ratio) × npix`

**PARTIE 2: Test avec Observation.py**
- Crée des observations réelles
- Renormalise les inputs
- Vérifie que les outputs matchent

### Résultats attendus

Si tout match après renormalisation → **Le calcul du SNR est correct !**

Les ratios doivent être ~1.0 (± 15%) :
- `Signal_renorm_final / signal_astropy ≈ 1.0`
- `sky_renorm_final / sky_astropy ≈ 1.0`

## Quelle convention utiliser ?

Ça dépend de ce que tu veux tester :

### Pour un spectrographe à fente
- **Generic ETC** a du sens : le sky est intégré sur toute la fente
- C'est réaliste : le bruit de fond vient de toute la fente, pas juste du pixel où est la source

### Pour comparer avec Astropy
- Faut renormaliser pour comparer des pommes avec des pommes
- Ou utiliser les mêmes conventions d'aire

## Conclusion

**Ce n'est PAS un bug**, c'est une différence de **définition** !

Generic ETC répond à la question : "Quel SNR j'obtiens pour détecter ma source dans un spectre ?"
Astropy répond à : "Quel SNR par pixel ?"

Pour cross-valider, il faut choisir quelle question on pose et renormaliser en conséquence.
