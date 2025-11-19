# ETC Development TODO List

## High Priority - Core Functionality

### IFS (Integral Field Spectrograph) Mode
- [ ] Verify SNR calculation is correct for IFS vs slits (1st plot: real SNR integrating several fibers/slicers; 2nd plot: SNR for single fiber)
- [ ] Ensure slit width and length have different impacts on instruments (one dimension affects resolution and mixes with sky)
- [ ] IFS mode should not cut flux by slit based on PSF at slit level (flux goes to next slit) - except when source << slit
- [ ] Verify flux integration: sky and flux integrated based on slit size (same impact as binning)
- [ ] Check: smaller fiber = less flux per fiber, but flux per arcsec² on image stays the same
- [ ] IFS mode should not impact SNR when slit width >> sigma source (but flux loss occurs when cutting)
- [ ] Verify SNR per spaxel makes sense for IFS
- [ ] Fix issue with fiber spectro where spaxel should be same in both directions
- [ ] Slicer does not work for MUSE

### SNR Calculations
- [ ] Understand why SNR is low for MUSE narrow (even when calculated per source)
- [x] Stacking: compute noise by dividing by sqrt(N) instead of computing N images
- [ ] Verify number of pixels calculation: total flux divided by number of pixels vs flux per pixel
- [ ] Ensure total flux and flux per pixel calculations give same result
- [ ] Should measure resolution for diffuse source (accounting for slit width)
- [ ] Understand factor 1.5 difference between fit and sigma source

### Throughput and Wavelength Dependencies
- [x] Account for QE even when throughput curve is given (normalize curve first, then multiply by QE × Atm × Throughput)
- [x] Verify counts don't change significantly when adding/removing QE
- [x] Account for Throughput FWHM
- [ ] Change file system to check for throughput curve when changing instrument
- [ ] Throughput not taken into account in the image (needs fixing)

## Medium Priority - Features & Improvements

### Spectra and Data Input
- [ ] Add ability to upload custom spectrum
- [ ] Understand why emission line evolution is weird
- [ ] Fix: flux changes when switching between cube and non-cube source
- [ ] Ensure spectra overlap correctly (continuum subtraction when summing pixels)
- [ ] Take only emission from CGM
- [ ] Subtract continuum option
- [ ] Check spectra in annulus
- [ ] Add adaptive smoothing? (may be impossible)

### Imager Mode
- [x] Remove unnecessary parameters for imagers (spectra, throughput(λ), sky lines, atm(λ), delta(λ))
- [x] Remove equivalent width, lambda stack, observed lambda, mask PSF for imagers
- [x] Always remove spectrograph design parameters (R, slit dims, dispersion, IFS) for imagers
- [x] Hide X-axis parameters when in imager mode
- [x] Replace "Spectral image" with "image" tab name for imagers
- [ ] Add ability to use any filter for imagers
- [ ] Add option to increase filter size to see observation improvements
- [ ] Add all filters from LePhare
- [ ] Find what image to show
- [ ] Add magnitude display
- [ ] Add simulated image from COSMOS or simulation (cosmic web?)

### Display and UI
- [x] Show units: e-/pix, photons, or ADU for images and datacubes
- [x] Fix: preserve tab name when switching from imager to instrument
- [x] Use ylog for all plots
- [x] Fix log scale issues in histogram and filters
- [x] Fix filter visibility issue (e.g., GALEX)
- [x] Last subplot shows surface brightness limit per pixel and per resolution element
- [ ] Add ability to have two different pixel sizes (requires changing slider type, slit display, etc.)
- [ ] Add dispersion in title when it disappears

## Low Priority - Specific Cases

### Instrument-Specific Issues
- [ ] Fix: passage to LUMOS does not work
- [ ] Fix: starting with LUMOS fails: `FB = ExposureTimeCalulator(instrument="LUMOS")`
- [ ] Verify sky lines and atmosphere can be applied to all instruments (currently buggy)
- [ ] Add HSC instrument

### Edge Cases & Testing
- [ ] Avoid divergence with pixel scale
- [ ] When dispersion increases beyond line width, handle correctly
- [ ] When dispersion decreases significantly, verify summing behavior
- [ ] Verify Throughput_FWHM size appears correctly on plot

## Completed ✓

- [x] Fix error with slicer displaying circle (x=y issue)
- [x] Show e-/pix, photons, or ADU units for datacube
- [x] Preserve tab name when switching instruments
- [x] Remove unnecessary parameters for imagers
- [x] Hide parameters in X-axis for imagers
- [x] Resolve all issues and verify plot updates
- [x] Use throughput for band integration
- [x] Apply ylog to all plots
- [x] Fix log scale in histogram and filters
- [x] Fix filter visibility
- [x] Add surface brightness limit plot
- [x] Stacking implementation using noise division

## Notes & Ideas

### Parameters to Consider
- exposure_time (affects cosmic ray, CIC)
- Spectral_resolution
- dispersion
- equivalent width
- pixel_scale
- Slitwidth
- EM_gain
- Bandwidth

### Research Questions
- Fabry-Perot interferometer implementation?
- Different source profiles?
- Change stack size dynamically?

---

**Last Updated**: 2025-01-05
**Priority Legend**: High = Core functionality issues | Medium = Features & UX | Low = Edge cases & nice-to-haves
