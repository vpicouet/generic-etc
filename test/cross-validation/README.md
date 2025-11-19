# Generic ETC Cross-Validation

Ce dossier contient le notebook de validation croisée de l'**ETC générique** (Exposure Time Calculator).

## Fichiers

- **`ETC_cross_validation.ipynb`** : Notebook principal de validation
- **`README.md`** : Ce fichier (documentation de la validation)

---

## Méthodologie de Validation

Nous comparons **deux approches** pour le calcul du SNR (Signal-to-Noise Ratio) :

### 1. **Generic ETC** (notre implémentation)
- Classe `Observation` dans `notebooks/Observation.py`
- Modélise toutes les sources de bruit : photon noise, dark current, read noise, CIC, cosmic rays
- Support EMCCD avec gain et excess noise factor
- Pertes d'ouverture (slit/fiber losses) via PSF
- Résolution spectrale prise en compte

### 2. **Astropy SNR** 📚 (référence validée)
- Fonction `astropy.stats.signal_to_noise_oir_ccd()`
- **Standard de la communauté** bien testé et documenté
- Calcule le SNR pour détecteurs CCD optical/infrared
- **Documentation** : https://docs.astropy.org/en/stable/api/astropy.stats.signal_to_noise_oir_ccd.html

---

## Comparaison des Paramètres

### Paramètres Astropy SNR

| Paramètre Astropy | Description | Unités |
|-------------------|-------------|--------|
| `t` | Temps d'exposition | secondes |
| `source_eps` | Taux d'électrons de la source/sec (dans l'ouverture) | e⁻/s |
| `sky_eps` | Taux d'électrons du ciel/sec/pixel | e⁻/s/pix |
| `dark_eps` | Dark current | e⁻/s/pix |
| `rd` | Read noise | e⁻ |
| `npix` | Nombre de pixels dans l'ouverture | pixels |
| `gain` | Gain CCD (défaut=1.0) | e⁻/DN |

### Formule Astropy

$$
\text{SNR} = \frac{t \cdot \text{source\_eps}}{\sqrt{t \cdot (\text{source\_eps} + \text{npix} \cdot (\text{sky\_eps} + \text{dark\_eps})) + \text{npix} \cdot \text{rd}^2}}
$$

---

## Table de Comparaison : Generic ETC vs Astropy SNR

| Fonctionnalité | Generic ETC | Astropy SNR |
|----------------|-------------|-------------|
| **Signal source** | ✅ `Signal` (erg/cm²/s/arcsec²/Å) | ✅ Converti en `source_eps` |
| **Sky background** | ✅ `Sky` (erg/cm²/s/arcsec²/Å) | ✅ Converti en `sky_eps` |
| **Dark current** | ✅ `Dark_current` (e⁻/pix/hour) | ✅ `dark_eps` (e⁻/s/pix) |
| **Read noise** | ✅ `RN` (e⁻/pix) | ✅ `rd` (e⁻) |
| **Temps d'exposition** | ✅ `exposure_time` | ✅ `t` |
| **Nombre de pixels** | ✅ Calculé automatiquement | ✅ `npix` |
| **Throughput** | ✅ `Throughput × QE × Atmosphere` | ⚠️ Appliqué manuellement |
| **Aire de collection** | ✅ `Collecting_area` | ⚠️ Appliqué manuellement |
| **Pixel scale** | ✅ `pixel_scale` | ⚠️ Appliqué manuellement |
| **Dispersion spectrale** | ✅ `dispersion` | ⚠️ Appliqué manuellement |
| **Résolution spectrale** | ✅ `Spectral_resolution` | ⚠️ Via `npix` |
| | | |
| **EMCCD - EM gain** | ✅ `EM_gain` | ❌ Non supporté |
| **EMCCD - CIC** | ✅ `CIC_charge` | ❌ Non supporté |
| **EMCCD - Smearing** | ✅ `smearing` | ❌ Non supporté |
| **EMCCD - Excess noise** | ✅ ENF = √2 | ❌ Non supporté |
| **PSF losses** | ✅ `flux_fraction` via erf() | ❌ Non supporté |
| **Cosmic ray losses** | ✅ `cosmic_ray_loss_per_sec` | ❌ Non supporté |
| **Extra background** | ✅ `extra_background` | ❌ Non supporté |

**Légende** :
- ✅ : Supporté nativement
- ⚠️ : Nécessite conversion/calcul manuel
- ❌ : Non supporté

---

## Résumé des Différences

### 🎯 **Generic ETC** - Forces

1. **Support EMCCD complet**
   - EM gain, CIC, smearing
   - Excess noise factor (√2)

2. **Pertes d'ouverture**
   - PSF (Gaussien via fonction erf)
   - Slit/fiber losses

3. **Cosmic rays**
   - Perte de pixels prise en compte

4. **Modes de calcul**
   - SNR per pixel
   - SNR per resolution element
   - SNR per source

5. **Interface unifiée**
   - Imaging + spectroscopie

### 📚 **Astropy SNR** - Forces

1. **Standard validé**
   - Fonction de référence de la communauté
   - Bien testée et documentée

2. **Simplicité**
   - Paramètres clairs
   - Formule transparente

3. **Fiabilité**
   - Partie d'Astropy (projet mature)
   - Large base d'utilisateurs

### ⚠️ **Limitations d'Astropy SNR**

- ❌ Pas de support EMCCD
- ❌ Pas de pertes PSF/slit
- ❌ Pas de cosmic rays
- ⚠️ Nécessite conversion manuelle flux → électrons

---

## Stratégie de Validation

### Phase 1 : Détecteur Simple ✅ (implémenté)

**Instrument** : GALEX FUV, UVEX FUV (MCP/CMOS, pas d'EMCCD)

**Tests** :
- SNR vs temps d'exposition
- SNR per pixel (`npix=1`)

**Critère de succès** : Accord ±10-15%

### Phase 2 : Résolution Spectrale (à implémenter)

**Tests** :
- SNR per resolution element (`npix > 1`)
- Vérifier calcul de `elem_size`

**Critère de succès** : Accord ±10-15%

### Phase 3 : EMCCD (limitation)

**Problème** : Astropy SNR ne supporte pas EMCCD

**Solutions possibles** :
1. Désactiver EM_gain dans Generic ETC pour comparer
2. Comparer uniquement le signal (pas le SNR complet)
3. Trouver un autre ETC de référence pour EMCCD

### Phase 4 : Analyse des Écarts

**Si écarts > 15%**, identifier la source :
- Conversion flux → électrons ?
- Pertes d'ouverture (`flux_fraction`) ?
- Cosmic ray losses ?
- Calcul du nombre de pixels ?

---

## Installation

```bash
pip install astropy numpy matplotlib pandas
```

## Usage

```bash
cd /Users/Vincent/Github/generic-etc/test/cross-validation
jupyter notebook ETC_cross_validation.ipynb
```

Exécuter les cellules séquentiellement.

---

## Résultats Attendus

### Accord attendu
- **±10-15%** pour détecteurs simples (CCD, MCP, CMOS sans EM gain)
- **±15-25%** pour cas complexes (PSF non-Gaussien, transmission variable)

### Écarts acceptables
Les écarts jusqu'à 25% sont normaux pour :
- Sources très faibles (proche limite de détection)
- Pertes PSF/slit (Generic ETC les modélise, Astropy non)
- Cosmic ray losses (Generic ETC les modélise, Astropy non)

---

**Auteur** : Vincent
**Date** : Novembre 2024
**Repository** : `/Users/Vincent/Github/generic-etc/`
