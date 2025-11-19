# Validation Images

Ce dossier contient les images de validation pour tous les instruments testés.

## Format des images

Chaque image (`{Instrument_Name}.png`) contient 5 tests de validation :

### Row 1 - Comparaison SNR
- **Col 1** : SNR vs Temps d'exposition
- **Col 2** : SNR vs Flux signal
- **Col 3** : SNR vs Fond de ciel
- **Col 4** : SNR vs Bruit de lecture

### Row 2 - Différences relatives (%)
- **Col 1** : Différence vs Temps d'exposition
- **Col 2** : Différence vs Flux signal
- **Col 3** : Différence vs Fond de ciel
- **Col 4** : Différence vs Bruit de lecture

### Row 3 - Dark current
- **Col 1-2** : SNR vs Courant d'obscurité
- **Col 3-4** : Différence vs Courant d'obscurité

## Interprétation

- Courbes **Generic ETC** (cercles) vs **Astropy SNR** (carrés pointillés)
- Lignes orange à ±15% = limite acceptable
- Titre de chaque subplot de différence montre le max de l'écart

## Critères

- ✓✓ EXCELLENT : < 10%
- ✓ PASSED : < 15%
- ✗ FAILED : ≥ 15%

## Génération

Images générées par la fonction `validate_instrument()` dans le notebook `ETC_cross_validation.ipynb`.
