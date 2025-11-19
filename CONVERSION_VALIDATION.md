# Validation des conversions Generic ETC

Ce document vérifie que toutes les conversions de Generic ETC sont cohérentes avec les formules physiques attendues.

## Méthodologie

Pour chaque paramètre, nous vérifions :
1. La formule théorique attendue
2. La formule utilisée par Generic ETC
3. La comparaison numérique sur GALEX FUV

## 1. Signal (flux source → électrons)

### Formule théorique

Le nombre d'électrons collectés depuis une source:

```
Signal_el = Flux × Collecting_area × Throughput × QE × Atmosphere
            × pixel_scale × dispersion × exposure_time
```

Où :
- `Flux` : flux source [erg/cm²/s/arcsec²/Å]
- `Collecting_area` : surface collectrice [m²] → [cm²]
- `Throughput` : transmission optique [sans unité]
- `QE` : efficacité quantique [sans unité]
- `Atmosphere` : transmission atmosphérique [sans unité]
- `pixel_scale` : échelle spatiale [arcsec/pix]
- `dispersion` : dispersion spectrale [Å/pix]
- `exposure_time` : temps d'exposition [s]

Conversion photons → électrons :
```
E_photon = h × c / λ  [erg]
N_photons = Flux × Area × Δλ × Δθ × Δt / E_photon
N_electrons = N_photons × QE × Throughput × Atmosphere
```

### Formule Generic ETC

Generic ETC calcule (voir Observation.py:2545):
```python
self.Signal_el = self.Signal_LU * self.factor_CU2el * self.exposure_time * self.flux_fraction_slit_applied
```

Avec `factor_CU2el` calculé selon le mode `test=True` (ligne 2522):
```python
self.factor_CU2el = self.factor_CU2el_tot
```

Où `factor_CU2el_tot` (ligne 2481):
```python
self.factor_CU2el_tot = (effective_area * arcsec2str *
                         np.minimum(Line_width, Bandwidth) *
                         source_size_arcsec_after_slit / pixels_total_source)
```

### Test à vérifier

- [ ] Vérifier que `effective_area` inclut bien `Collecting_area × Throughput × QE × Atmosphere`
- [ ] Vérifier la conversion flux → photons → électrons
- [ ] Vérifier l'intégration spatiale (arcsec²)
- [ ] Vérifier l'intégration spectrale (Å)
- [ ] Comparer avec calcul manuel

## 2. Sky background (flux ciel → électrons)

### Formule théorique

```
sky = Sky_flux × Collecting_area × Throughput × QE × Atmosphere
      × slit_spatial × slit_spectral × exposure_time
```

Où :
- `Sky_flux` : flux de fond de ciel [erg/cm²/s/arcsec²/Å]
- `slit_spatial` : dimension spatiale de la fente [arcsec]
- `slit_spectral` : largeur spectrale [Å]

**Important** : Le ciel est intégré sur toute la fente, pas juste la source!

### Formule Generic ETC

Generic ETC calcule (ligne 2539):
```python
self.sky = self.Sky_CU * self.factor_CU2el_sky * self.exposure_time
```

Avec `factor_CU2el_sky` (ligne 2534):
```python
self.factor_CU2el_sky = (effective_area * arcsec2str *
                         np.maximum(np.minimum(Line_width, Bandwidth), dispersion) *
                         slit_size_arcsec_after_slit / pixels_total_source)
```

**Différence clé** :
- Signal utilise `source_size_arcsec_after_slit` (taille de la source)
- Sky utilise `slit_size_arcsec_after_slit` (taille de la fente)

### Test à vérifier

- [ ] Vérifier que `slit_size_arcsec_after_slit` est bien la surface de la fente
- [ ] Vérifier la largeur spectrale (max avec dispersion pour au moins 1 pixel)
- [ ] Comparer avec calcul manuel
- [ ] Vérifier que sky > signal si la fente est plus grande que la source

## 3. Dark current (e⁻/pix/hour → électrons totaux)

### Formule théorique

```
Dark_current_f = Dark_current × exposure_time × n_pixels
```

Où :
- `Dark_current` : courant d'obscurité [e⁻/pix/hour]
- `exposure_time` : temps d'exposition [s]
- `n_pixels` : nombre de pixels (1 en mode "per pix")

Conversion hour → s :
```
Dark_current_f = Dark_current × exposure_time / 3600 × n_pixels
```

### Formule Generic ETC

Generic ETC calcule (ligne 2542):
```python
self.Dark_current_f = self.Dark_current * self.exposure_time / 3600 * self.number_pixels_used
```

### Test à vérifier

- [x] Formule correcte
- [ ] Vérifier `number_pixels_used = 1` en mode "per pix"
- [ ] Comparer avec calcul manuel

## 4. Read noise (e⁻ → électrons totaux)

### Formule théorique

Le bruit de lecture est déjà en électrons. En mode "per pix" avec 1 pixel:
```
RN_final = RN
```

Mais Generic ETC peut appliquer un facteur EMCCD.

### Formule Generic ETC

Generic ETC calcule (ligne 2547-2560):
```python
if self.EMCCD:
    # EMCCD mode avec EM_gain
    self.RN_final = self.RN_EMCCD / self.EM_gain
else:
    self.RN_final = self.RN
```

### Test à vérifier

- [ ] Vérifier que `RN_final = RN` pour instruments normaux (EM_gain=1)
- [ ] Vérifier le facteur EMCCD si applicable
- [ ] Comparer avec calcul manuel

## 5. Photon energy et conversion erg → photons

### Formule théorique

```
E_photon = h × c / λ  [erg]
```

Où :
- h = 6.62607015e-27 erg·s (constante de Planck)
- c = 2.99792458e10 cm/s (vitesse de la lumière)
- λ : longueur d'onde [cm]

Pour convertir flux [erg/cm²/s/arcsec²/Å] en photons/s/cm²/arcsec²/Å :
```
Flux_photons = Flux_erg / E_photon
```

### Formule Generic ETC

Generic ETC utilise (ligne 2406-2412):
```python
h = 6.62606896e-27  # erg s
c_light = 299792.458  # km/s
c = c_light * 1e5  # cm/s
wavelength_cm = self.wavelength / 1e7  # nm → cm
photon_energy = h * c / wavelength_cm  # erg
arcsec2str = ((np.pi / 180.0) / 3600.0) ** 2  # arcsec² → sr
```

Et la conversion (ligne 2447-2453):
```python
effective_area = (self.Collecting_area * 1e4 * self.Throughput *
                  self.QE * self.Atmosphere / photon_energy)
```

**Important** : `effective_area` inclut déjà la conversion erg → photons !

### Test à vérifier

- [ ] Vérifier les constantes physiques
- [ ] Vérifier la conversion wavelength nm → cm
- [ ] Vérifier que effective_area = A × T × QE × Atm / E_photon
- [ ] Comparer avec calcul manuel

## 6. Aires et intégration spatiale

### Source size vs Slit size

Generic ETC calcule (ligne 2472-2473):
```python
source_size_arcsec_after_slit = (np.minimum(Size_source * fwhm_sigma_ratio, Slitlength) *
                                 np.minimum(Size_source, Slitwidth))
slit_size_arcsec_after_slit = Slitwidth * Slitlength
```

**Comportement** :
- Pour le signal : on intègre sur `min(source, slit)` → ne capte que la source
- Pour le sky : on intègre sur toute la fente → capte tout le bruit de fond

### Test à vérifier

- [ ] Vérifier `source_size_arcsec_after_slit < slit_size_arcsec_after_slit` si source < slit
- [ ] Vérifier le facteur `fwhm_sigma_ratio = 2.35`
- [ ] Comparer les aires calculées avec les valeurs attendues

## 7. Division par pixels_total_source

Generic ETC divise les facteurs de conversion par `pixels_total_source` (ligne 2395-2399):
```python
sigma_x_pix = Size_source / pixel_scale / fwhm_sigma_ratio
sigma_y_pix = Line_width / dispersion / fwhm_sigma_ratio
pixels_total_source = 2 * np.pi * sigma_x_pix * sigma_y_pix
```

**Interprétation** : Generic ETC calcule le signal/sky **moyenné** sur les pixels de la source.

### Test à vérifier

- [ ] Vérifier le calcul de `pixels_total_source`
- [ ] Vérifier que cette division est cohérente avec le mode "per pix"
- [ ] Comparer avec calcul manuel

## Tests à effectuer

### Test 1: Calcul manuel complet pour GALEX FUV

Paramètres:
- wavelength = 152.8 nm
- Collecting_area = 0.196 m²
- pixel_scale = 1.5 arcsec/pix
- dispersion = 2.5 Å/pix
- Throughput = 0.08
- QE = 0.12
- Atmosphere = 1.0
- Signal = 5.6e-19 erg/cm²/s/arcsec²/Å
- Sky = 1.53e-19 erg/cm²/s/arcsec²/Å
- Dark_current = 0.049 e⁻/pix/hour
- Size_source = 4.0 arcsec
- Slitwidth = 4.0 arcsec
- Slitlength = 20.0 arcsec
- Line_width = 2.0 Å
- Bandwidth = 250.0 Å
- exposure_time = 1000.0 s

Calculer manuellement:
1. E_photon
2. effective_area
3. source_size_arcsec_after_slit
4. slit_size_arcsec_after_slit
5. pixels_total_source
6. factor_CU2el_tot
7. factor_CU2el_sky
8. Signal_el
9. sky
10. Dark_current_f

Comparer avec Generic ETC.

### Test 2: Variations paramétriques

Vérifier que les conversions varient correctement avec:
- Collecting_area (linéaire)
- Throughput (linéaire)
- QE (linéaire)
- wavelength (1/λ pour E_photon)
- pixel_scale (linéaire spatiale)
- dispersion (linéaire spectrale)

## Résultats

*À compléter après les tests*

## Conclusion

*À compléter après validation*
