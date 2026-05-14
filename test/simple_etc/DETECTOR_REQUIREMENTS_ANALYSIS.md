# Detector Requirements Analysis for Ultra-Faint Sources

## Source Parameters

| Parameter | Value | Unit |
|-----------|-------|------|
| Surface brightness | 5.6e-19 | erg/cm²/s/arcsec²/Å |
| Flux (Lyman units) | 2400 | LU/Å |
| Line width | 1 | Å |
| Aperture | 10×10 = 100 | arcsec² |
| Acquisition time | 90 | hours |

## Expected Signal

- Total signal: 1-2 electrons distributed over 50 pixels
- Signal per pixel: ~0.03 e-/pix over 90h
- Signal rate: ~0.0003 e-/pix/hour

## Detector Requirements for SNR = 3

### Noise Budget

To achieve SNR~2.5 with this extremely faint signal, the following maximum values are required:

| Component | Maximum Value | Unit |
|-----------|---------------|------|
| Dark + Sky | 0.01 | e-/pix/hour |
| CIC | 0.001 | e-/pix/frame |
| RN | 50/2000=0.025 | e-/pix/frame |

Which would basically give 30% of noise for each source + 10% for schot noise

Reducing further the RN could allow to reduce the gain while keeping the necessary G~50xRN which can improve CIC and CTE


read noise / 
0.5 (10 secondes)
0.01 (2 minutes)

## Optimistic Assumptions Not Included

This analysis does not account for several degrading factors:

| Factor | Estimated Impact |
|--------|-----------------|
| Cosmic ray loss | -30% duty cycle |
| Excess Noise Factor (ENF) | ×sqrt(2) on noise (EMCCD) |
| Photon counting threshold loss | -5 to -30% photons |
| Smearing | reduced pc loss |
| Stacking inefficiencies | Signal loss / noise increase |
| Higher CR rate in flight | Additional duty cycle loss |
| Contamination (zodiacal, airglow) | Increased sky background |


## Recommendations

### 1. Required Detector Specifications

```
- ultra Dark current: << 0.01 e-/pix/hour 
- CIC: << 0.001 e-/pix/frame (high-quality EMCCD)
- Read noise: < 0.01 e-/pix (with high EM gain)
- Background +straylight lower than dark
```
Laboratory validation required before flight?
- Cosmic ray rejection efficiency?



**Key Finding**: With nominal parameters (Dark+Sky = 0.01 e-/pix/h), SNR = 3 is not achievable in 90h. Required improvements:
- Reduce background by factor of 10 (Dark+Sky ≤ 0.001 e-/pix/h), OR
- Increase acquisition time by factor of 9 (810h), OR
- Improve light collection (larger telescope)


qCCD:
- Would need better than <0.2, what does it mean ? 0.1 is much better!!
- Clearly would mean going for as long exposure time as possible almost
- If we assume similar dark (actually we don'; t need lower than 0.1 if RN is 0.5, 0.01 if RN is 0.1) as what they send (with still good CTE) their second option is much better for a source of noise we have. we could take 1h images, if it is physically possible
   -  of course currently is not feasible due to sky rotation
-  Could be interesting to know the cosmic ray impact, number of pixels loss per impact, but it's really second order
-  QE would be similar if we coat them ? We can also have a rejection filter?
- My fear is really the readiness level of such a detector for a space mission.
Frame transfer CCD, 2min 2min 