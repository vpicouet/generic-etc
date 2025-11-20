# Validation Generic ETC vs Keck NIRES ETC

## Configuration de test

**Instrument: NIRES**
- Type: Long slit spectrograph
- Wavelength: 1500 nm (1.5 μm)
- Spectral resolution: R = 2700
- Slit: 0.55" × 18"
- Pixel scale: 0.123"/pix

**Observation**
- Exposure time: 3600 s
- Seeing: 0.75"
- Signal: 5.6e-19 erg/cm²/s/arcsec²/Å → AB mag 23

## Résultats

| Paramètre | Generic ETC | NIRES ETC | Ratio |
|-----------|-------------|-----------|-------|
| Signal | 138 e⁻ | 1947 e⁻ | 0.071 |
| Noise | 127 e⁻ | 36663 e⁻ | 0.0035 |
| SNR | 1.09 | 0.075 | 14.5 |

## Différences identifiées

### 1. Throughput

- Generic ETC: 0.31 (peak value)
- NIRES ETC file at 1.5 μm: 0.1117
- Ratio: 2.8×

### 2. Intégration spatiale

**Generic ETC:**
```python
source_size_arcsec = Size_source² = 0.49² = 0.24 arcsec²
```

**NIRES ETC:**
```python
spatial_cov = (slit_width * seeing) / (π * seeing²)
            = (0.55 * 0.75) / (π * 0.5625)
            = 0.233
```
Cette valeur est une fraction, pas une aire !

### 3. Erreur possible dans le calcul de bruit NIRES

Le code NIRES semble avoir une erreur :
```python
noise_int = (sig_src*exp_time + num_pix_slt*dark_current*exp_time +
             (num_pix_slt-num_pix_src)*sky_bkg_interp*exp_time +
             num_pix_slt*read_noise**2*exp_time)**(1/2)
```

Le terme `read_noise**2 * exp_time` est incorrect - le read noise ne devrait pas être multiplié par le temps d'exposition.

### 4. Nombre de pixels

- Generic ETC: pixels_total_source = 9.52 pix
- NIRES ETC: num_pix_slt = 18 * 0.55 / 0.123² = 655 pix

NIRES utilise toute la fente (655 pixels) pour le bruit, pas juste la source.

## Conclusion

La comparaison n'est pas directe car :

1. **Unités différentes** : Generic ETC utilise flux en erg/cm²/s/arcsec²/Å, NIRES utilise magnitudes
2. **Modèle d'intégration différent** : Generic intègre sur la source, NIRES sur toute la fente
3. **Throughput différent** : Peak (0.31) vs valeur à λ (0.11)
4. **Possible bug** dans le calcul de bruit NIRES (RN² × t_exp)

## Recommandation

Pour une comparaison valide :
- Utiliser la même valeur de throughput (0.11 à 1.5 μm)
- Vérifier le calcul de bruit NIRES avec les développeurs
- Comparer avec les observations réelles plutôt qu'un autre ETC

## Statut

⚠ **Non validé** - Différences de modélisation trop importantes pour comparaison directe
