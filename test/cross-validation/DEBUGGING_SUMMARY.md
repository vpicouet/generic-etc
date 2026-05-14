# KCWI ETC Validation - Debugging Summary

## Context
Validating Generic ETC's new "per Source per 2λpix" mode against KCWI's official ETC. Goal: Perfect SNR agreement when using identical inputs.

## Current Status: ~12% Systematic Difference

**Latest results** (1 config: Small-BL, matched inputs):
- SNR ratio (Generic/KCWI): **0.88 (constant across all parameters)**
- Generic ETC gives **~12% lower SNR** than KCWI
- Difference is **extremely constant** (-10.5% to -14.7%) → suggests systematic factor, not calculation error

## What We Fixed

### 1. Spectral Integration ✓
**Problem**: Generic ETC was integrating over `Line_width` instead of fixed spectral bins.
**Solution**: Modified `Observation.py` to use adaptive spectral bins: `2^(floor(log2(Slitwidth/pixel_scale)))` pixels
- Small slicer: 2 pixels
- Medium slicer: 4 pixels
- Large slicer: 8 pixels

### 2. Spatial Integration ✓
**Problem**: Generic ETC used FWHM factor (2.35) for spatial integration.
**Solution**: Changed to use `Size_source²` directly (= seeing² for point sources) in `Observation.py:2514-2519`

### 3. Throughput Calibration ✓
**Problem**: Generic ETC used fixed throughput (0.34×0.85×0.9), KCWI uses grating-specific values.
**Solution**: Pass KCWI's empirical efficiency directly via `Throughput` parameter.

### 4. Sky Background Input ⚠️
**CRITICAL FINDING**: KCWI ETC **ignores the input sky parameter** and uses its own sky model from `mk_sky.dat`!
- Test showed: varying `sky_erg` input 0.5× to 2× → KCWI SNR unchanged
- Must extract KCWI's sky value and pass it to Generic ETC for fair comparison
**Solution**: Modified comparison script to use KCWI's sky model:
```python
c_s_array = sky_cts(w_kcwi, grating, 1.0)
sky_per_s = c_s_array[idx]
sky_erg = (sky_per_s / (eff * A_geo)) * E_photon
```

### 5. Dark Current ⚠️
**FINDING**: KCWI ETC has **NO dark current** in its calculation (line 233 of `etc_kcwi.py`: `snr = c_o/sqrt(c_s+c_o+c_r)` - no dark term).
**Solution**: Set `Dark_current=0.0` in Generic ETC for validation.

## Current Comparison Setup

**File**: `/Users/Vincent/Github/generic-etc/test/cross-validation/plot_kcwi_snr_ratios_fast.py`

**Matched parameters**:
- ✓ Throughput: Use KCWI's empirical efficiency per grating
- ✓ Sky: Use KCWI's sky model (extracted from `sky_cts()`)
- ✓ Dark: 0.0 (KCWI doesn't have dark)
- ✓ Read noise: 2.7 e⁻
- ✓ Pixel counting: Matches exactly (29.155 pixels for Small-BL-seeing=0.75)

**Test results** (`kcwi_snr_comparison.csv`):
```
config  param   value      pct_diff    snr_gen  snr_keck  signal_total  sky_total
S-BL    t_exp   1000.0     -14.67%     0.654    0.767     14.64         265.16
S-BL    t_exp   3600.0     -11.83%     1.504    1.706     52.72         954.59
S-BL    t_exp   10000.0    -10.78%     2.663    2.985     146.45        2651.64
S-BL    flux    2.8e-19    -11.86%     0.760    0.863     26.36         954.59
S-BL    flux    5.6e-19    -11.83%     1.504    1.706     52.72         954.59
S-BL    rn      1.0        -10.53%     1.632    1.824     52.72         954.59
S-BL    rn      2.7        -11.83%     1.504    1.706     52.72         954.59
S-BL    rn      5.0        -14.07%     1.261    1.468     52.72         954.59
```

## The Mystery: Where Does the 12% Come From?

### What We Verified ✓
1. **Pixel counting matches exactly**: 29.155 pixels (both ETCs)
2. **Object signal matches**: Generic=52.72 e⁻ vs KCWI=65.44 e⁻ → ratio 0.806
3. **Read noise² matches exactly**: 212.54 e² (both ETCs)
4. **Sky electrons** from debug script:
   - KCWI: 1192.966 e⁻ (from sky_cts)
   - Generic: 954.59 e⁻ (from CSV) → ratio 0.800

### The Smoking Gun: Sky Total
From `debug_counts_comparison.py` output:
```
KCWI sky:     1192.966 e⁻
Generic sky:  4,295,652 e⁻  ← 3600× higher!!!
```

**Ratio**: 4295652 / 1192.966 = **3600.8** = exactly `t_exp`!

**BUT** in the CSV from the fast script:
```
Generic sky_total: 954.59 e⁻  (much better!)
```

This suggests the sky calculation was partially fixed but there's still something wrong.

### Detailed Counts Comparison (t_exp=3600s)

**KCWI counts**:
- Signal: 65.441 e⁻
- Sky: 1192.966 e⁻
- RN²: 212.536 e²
- SNR = 65.441 / sqrt(65.441 + 1192.966 + 212.536) = 1.706

**Generic counts** (from CSV):
- Signal: 52.722 e⁻ (0.806× KCWI)
- Sky: 954.589 e⁻ (0.800× KCWI)
- RN²: 212.536 e² (perfect match)
- Expected SNR = 52.722 / sqrt(52.722 + 954.589 + 212.536) = 1.509
- Actual SNR from CSV: 1.504 ✓ matches calculation

**Signal ratio**: 52.722 / 65.441 = **0.806**
**Sky ratio**: 954.589 / 1192.966 = **0.800**
**Both are ~0.80** → suggests a common ~20% reduction factor

## Hypotheses to Investigate

### Hypothesis 1: Factor_CU2el Calculation Error
The signal/sky conversion factor may be missing a multiplicative term.

**Check**:
```python
# In Observation.py, "per Source per 2λpix" mode around line 2539-2553
# Look at factor_CU2el_tot calculation
# Should be: A_T × Throughput × arcsec2str × Ω_spatial × Δλ_spectral
```

**Specific issue to check**: Is `arcsec2str` correctly applied? Is there a factor of `t_exp` missing or added?

### Hypothesis 2: Sky Calculation Still Wrong
Even though we extract KCWI's sky and convert to erg, the conversion may be incorrect.

**The conversion chain**:
1. KCWI gives: `sky_cts = eff × A_geo × t_exp × sky_photons × spatial_bin × spectral_bin`
2. For 1 second: `sky_per_s = sky_cts(t=1.0)` (electrons/s/Å/arcsec²)
3. Reverse to photons: `sky_photons_per_s = sky_per_s / (eff × A_geo)`
4. Convert to erg: `sky_erg = sky_photons_per_s × E_photon`

**Check**: Is this conversion correct? Compare `sky_erg` value used in Generic ETC vs expected.

### Hypothesis 3: Spatial/Spectral Bin Area Mismatch
KCWI uses:
- Spatial bin: `seeing × seeing` = 0.5625 arcsec²
- Spectral bin: `2 × 0.625` = 1.25 Å
- **Total**: 0.5625 × 1.25 = 0.703 arcsec²·Å

Generic ETC should use equivalent via `factor_CU2el`.

**Check in Observation.py around line 2540**:
```python
spectral_bin_angstrom = spectral_bin_pixels * self.dispersion * 10.0
self.factor_CU2el_tot = (self.effective_area * self.arcsec2str *
                         self.source_size_arcsec_after_slit * spectral_bin_angstrom /
                         self.pixels_total_source)
```

Does `source_size_arcsec_after_slit × spectral_bin_angstrom` equal 0.703 arcsec²·Å?

## Next Steps for Debugging

1. **Add debug prints** to `Observation.py` in "per Source per 2λpix" mode:
   ```python
   print(f"DEBUG factor_CU2el: {self.factor_CU2el[self.i]:.6e}")
   print(f"DEBUG spatial_bin: {self.source_size_arcsec_after_slit:.6f} arcsec²")
   print(f"DEBUG spectral_bin: {spectral_bin_angstrom:.6f} Å")
   print(f"DEBUG total_area: {self.source_size_arcsec_after_slit * spectral_bin_angstrom:.6f} arcsec²·Å")
   print(f"DEBUG expected_area: {0.5625 * 1.25:.6f} arcsec²·Å")
   ```

2. **Compare signal conversion**:
   - KCWI: `signal = eff × A_geo × t_exp × flux_photons × 0.703`
   - Generic should give equivalent via `Signal_el × pixels_total`

3. **Check if `Sky` parameter units are correct**:
   - Expected: erg/cm²/s/Å/arcsec²
   - Verify the conversion from KCWI's sky model matches this

4. **Look for factor of 0.8 or 1.25 somewhere**:
   - 0.8 = 4/5
   - 1.25 = 5/4
   - Might be related to spectral bin calculation or pixel counting

## Files Modified

1. **`/Users/Vincent/Github/generic-etc/notebooks/Observation.py`**
   - Lines 2514-2519: Spatial integration (use Size_source² directly)
   - Lines 2539-2553: Adaptive spectral bin calculation
   - Lines 2445-2468: Adaptive pixel counting for SNR

2. **`/Users/Vincent/Github/exposureTimeCalculator/etc_kcwi.py`**
   - Lines 298-369: Added `calculate_snr()` function for easy access to KCWI's SNR

3. **`/Users/Vincent/Github/generic-etc/test/cross-validation/plot_kcwi_snr_ratios_fast.py`**
   - Modified to use matched inputs (KCWI's sky, dark=0)
   - Quick test mode: 1 config, 3 points per parameter
   - Saves CSV with detailed breakdown

## Key Insight

The **constant ~12% difference** across all parameter variations means:
- ✓ The SNR **calculation formula** is correct (shape of curves matches)
- ✗ There's a **normalization/scaling factor** wrong by ~0.88
- The factor affects both signal AND sky equally (~0.80 ratio for both)
- This points to `factor_CU2el` or the area calculation

**Most likely culprit**: The calculation of spatial_bin × spectral_bin area in Generic ETC doesn't exactly match KCWI's 0.703 arcsec²·Å.
