
# A generic spectrograph model for exposure time calculations, trade studies, and instrument/observations optimization

> Vincent Picouet
The California Institute of Technology, 1200 E. California Blvd., Pasadena, CA, 91125, USA
> 

<aside>

**Abstract**:

This paper introduces a basic tool at the frontier between an exposure time calculator (ETC) and an instrument model, aiming for broader applicability, with a high level of genericity. Its goal goes slightly beyond the nominal use of traditional ETCs which predict the signal-to-noise ratio (SNR) of a source given some instrument parameter. As others ETC, it estimates the expected SNR and predicts observations to improve the reduction pipeline or adapt detection strategy, but more importantly, it enables analysis of the evolution of the SNR with all the instrument/observation parameters, allowing to examine the instrument efficiency, explore the SNR evolution under different scenarios and run different trade studies. The ETC is linked to an online database to allow any scientist to add their own spectrograph instrument. The current version already encompasses more than 20 instruments, some with several channels or configurations. This tool can be valuable for instrument comparisons and trade-off analysis to either optimize the instrument design or the observation strategy. Despite being a personal initiative with modest resources, it serves as an illustrative example of development simplicity with collaborative database utilization. Observations predictions have been cross-checked with ETC-42 based on several spectrograph designs. This article briefly outlines its development philosophy and role in facilitating trade-off analyses for future instrument developments.

</aside>

---

> Keywords: Spectrographs, Astronomical Instruments, Astronomical Observations, Exposure Time Calculator, Simulation, Spectroscopy
> 

---

<aside>

Links: [GitHub](https://github.com/vpicouet/spectro-imager-etc), [Database](https://docs.google.com/spreadsheets/d/1Ox0uxEm2TfgzYA6ivkTpU4xrmN5vO5kmnUPdCSt73uU/edit?usp=sharing), [Online-ETC](https://spectro-imager-etc-fbcad4e5849d.herokuapp.com)

</aside>

# Introduction

Exposure Time Calculators (ETCs) are critical tools for optimizing observational strategies in astronomy. Traditionally, these calculators have been tailored to specific instruments, however, there has been some development of universal ETCs capable of accommodating a variety of imagers or spectrographs \citep{ETC42}. 

These ETCs are mostly released after the commissioning phase and addressed to observers. In the development phase, instrument builders and scientist usually rely on more or less complex instrument models to optimize their designs or some of their subsystems. Even though this development brings a lot of interesting considerations, it requires significant resources. Creating an initial model during the exploration phase, evolving it into a more complex tool during the instrument development, and finally building an ETC post-commissioning can be resource-intensive, especially for experimental spectrographs or low-budget projects like suborbital missions.

To address these challenges, we present a simple generic tool that can allow to address several of these specific needs bridging some of the gaps between an ETC and an instrument model. This straightforward ETC enables instrument builders to both explore how spectrograph designs influence observations and analyze how some more specific changes in instrument parameters impact the SNR , but also allows observers to optimize their observations strategies when the instrument is finalized.

The ETC offers insights into the evolution of signal and noise contributions, as well as SNR, with respect to various instrument parameters. This not only provides a quick overview of what drives SNR but also highlights potential trade-offs and optimization opportunities. The ETC is designed to accommodate any spectro-imager, and scientists can easily contribute new instruments or configurations via an online database for direct integration with the tool.

[In](https://arxiv.org/pdf/2401.08152) this paper, we describe the ETC's goals, input and output in Section 2 and its deployment philosophy in Section 3. In Section 4 we describe the ETC’s design and architecture. Section 5 describes the FIREBall case, performed to analyze different SNR tradeoffs and optimize the instrument design and observations strategy accordingly. Finally, we discuss our conclusions and potential future developments in Section 6.

# Spectrograph model

## **Important features and goals**

The ETC distinguishes itself from a specific ETC or complete instrument simulator, although some parts of its functionality are similar. Its goal is to enable tradeoff analysis in instrument design or observation strategy and compare different instruments capability. The tool's development inevitably represents a compromise between realism, cost, and user experience. 

The key features of the spectro-imager ETC are summarize below:

**Database:**

- Simplicity to contribute any instrument in the [database](https://docs.google.com/spreadsheets/d/1Ox0uxEm2TfgzYA6ivkTpU4xrmN5vO5kmnUPdCSt73uU/edit?usp=sharing) ([Figure](https://www.notion.so/A-generic-spectrograph-model-for-exposure-time-calculations-trade-studies-and-instrument-observati-10493411461380889a97d23f7e382d14?pvs=21)).
- Allows to add several instrument configurations (see [Section](https://www.notion.so/A-generic-spectrograph-model-for-exposure-time-calculations-trade-studies-and-instrument-observati-10493411461380889a97d23f7e382d14?pvs=21)).
- Allows to compare instruments’ capability easily ([Figure](https://www.notion.so/A-generic-spectrograph-model-for-exposure-time-calculations-trade-studies-and-instrument-observati-10493411461380889a97d23f7e382d14?pvs=21)) for different science cases.

**Exposure time calculator:**

- Works for almost any type of spectrograph or spectro-imager: Slit spectrograph, slitless spectrograph, prism spectrograph, fiber integral field spectrograph, slicer integral field spectrograph.
- Facilitates the visualization of the evolution of (noise, contributions, SNR, limiting surface brightness) with any instrument parameter, enabling trade-off analysis or instrument/observation optimization.
- Rapid computation: ETCs that take too long to execute can prevent users from properly assessing their observation parameter space and understand what parameters drive the final SNR. This ETC specifically allows this. In the image simulation mode, it will only consider part of the FOV to perform the calculation ($100\times 500$ pixels, spatially and spectrally respectively). This allows to explore more easily the instrument or observation parameter space to try to optimize it.
- Users can easily upload their own λ-dependent throughput spectra (Section 2), atmospheric transmission, or sky/background.
- For simplicity, all unused parameters are automatically hidden in the GUI.
- Possibility to export simulated images to get used to them, develop or improve their instrument’s reduction pipeline.

Its key attributes lie in its versatility, intuitive user interface, and diverse plotting options. These different features can then can make this tool appropriate for relatively small projects that lack the resources to develop their own ETC (like suborbital projects) or any larger looking to crosscheck their results.The ETC is useful for both observers and instrument scientists, and it also has educational value for teaching instrument design or observational spectroscopy astronomy.

## Parameters input

### λ-independent parameters

All instrument and observation parameters are shown in Figure 1 and are described in appendix. Their description and unit can also be found in the instrument spreadsheet as well as when the user put the cursor on each widget name. 

The different instrument parameters used can all be used in the $x$-axis of the SNR plot ([Figure](https://www.notion.so/A-generic-spectrograph-model-for-exposure-time-calculations-trade-studies-and-instrument-observati-10493411461380889a97d23f7e382d14?pvs=21)) to analyze the impact on the SNR. No other parameters than the one accessible in the spreadsheet (and the λ-dependant one discussed in the following Section) enter the SNR calculation and image/cube simulator. This choice prevents from accounting complex behaviour such as field distortion, special PSF shape or its evolution in the field.

![All the parameters are accessible for the exposure calculator. For simplicity, all the parameters that should not be used disappear automatically](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/ba2bc345-4044-4488-9d51-572ae9b92384/Parameters.jpg)

All the parameters are accessible for the exposure calculator. For simplicity, all the parameters that should not be used disappear automatically

### Atmosphere and other λ-dependant parameters

The ETC supports both flat and λ-dependent models for atmospheric absorption, sky emission, and end-to-end instrument throughput (including quantum efficiency). These models are automatically convolved with the instrument's spectral resolution. All input files must adhere to the comma-separated values (CSV) format, with wavelengths provided in nanometers and efficiency values between 0 and 1.

**Atmospheric absorption**

A standard atmospheric transmission model (λ-dependent) is applied for ground-based instruments using the python package pwv_kpno. Users also have the option to introduce custom transmission curves to better match their specific observation bandpass. For instance, FIREBall’s specific UV transmission curve at 40 km altitude is loaded from a file located at `Instruments/Instrument_name/Atmosphere_transmission.csv`. This feature enhances the ETC's flexibility across different observational platforms, from ground to suborbital.

![Atmospheric transmission assessed for ground and suborbital instrument.](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/ebece987-35c5-4946-89bb-ac2ded7482e1/image.png)

Atmospheric transmission assessed for ground and suborbital instrument.

**Sky emission**

A flat or generic sky emission model based on [UVES](https://www.eso.org/observing/dfo/quality/UVES/uvessky/sky_3460_1.html) observations (see following Figure) is natively available. Similarly as the atmospheric absorption, users can easily push their sky emission model in the instrument folder: `Instruments/Instrument_name/Sky_emission_lines.csv`

![Sky emission lines model used in the ETC when Sky(λ) is selected. ](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/29461ef0-2c60-40b8-917a-794cabb6598f/image.png)

Sky emission lines model used in the ETC when Sky(λ) is selected. 

The sky emission line button appears for any ground instrument or any other instrument that has a sky emission line. When it is checked, it will use the user-given sky model or, if there is none, the generic UVES observations. 

This flexibility extends to any background emission, not just sky. For instance, even orbital spectrographs can integrate non-flat background emission spectra, allowing the ETC to treat this data as it would a sky spectrum, allowing simulating in a range of environments.

**Throughput and quantum efficiency**

The ETC allows users to input throughput curves either as user-provided data or modeled as Gaussian profiles (based on the `Throughput_FWHM` parameter in the database). Users can easily upload an instrument-specific λ-dependent throughput curve, by pushing a `.csv` file to the GitHub repository : `Instruments_throughput/Instrument_name/Throughput.csv` 

If done, the `Throughput_FWHM` will be hidden and the curve will be maximum-normalized to the value of  $Throughput \times QE$, so that changing any of these parameters modifies the curve.

Users are encouraged to upload their instrument throughput, but if no table is added, the code will default to using the Throughput_FWHM value in the spreadsheet.

## Output plots

Currently, there are two main visualization options: the SNR evolution visualization and an image or cube (integral field units: slicers, fiber IFU, etc) simulation, both capable of accommodating any instruments stored in the [database](https://docs.google.com/spreadsheets/d/1Ox0uxEm2TfgzYA6ivkTpU4xrmN5vO5kmnUPdCSt73uU/edit?pli=1#gid=2066284077).

These visualizations are dynamic, adjusting to any parameter changes made via the intuitive interface widgets.

### SNR evolution

The SNR plot is crafted to offer a straightforward depiction of noise budgets concerning various variables that can be fine-tuned or mitigated to enhance instrument sensitivity. In the top panel, the noise from distinct sources (Signal, Dark, Sky, Read-noise, clock induced charge, extra -background) is presented in electrons per pixel. The middle panel provides the average electron-per-pixel value for each component (pre-stacking). The last plot outlines the relative fractions of all noise sources per N frames over the total acquisition time and the resulting SNR. All plots can be adjusted to display per-pixel or per-resolution element data, offering users a clear view of how SNR evolves across different configurations.

When the SNR exhibit a local optimal it is then straightforward to analyze the evolution of the different contribution and noises to see what drives this SNR optimum.

[Evolution of noise (top plot), individual contributions (second), SNR (third) and surface brightness limit (bottom) with respect to exposure time. Any instrument/source parameter can be used to analyze the evolution of these metrics and potentially optimize the instrument or observation strategy.](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/83dcfc55-499f-4b4a-9e52-1f9556713955/Enregistrement_decran_le_2024-09-16_a_19.31.15.mp4)

Evolution of noise (top plot), individual contributions (second), SNR (third) and surface brightness limit (bottom) with respect to exposure time. Any instrument/source parameter can be used to analyze the evolution of these metrics and potentially optimize the instrument or observation strategy.

### Image simulator

The image simulator utilizes the various parameters and a specific source (continuum, emission line or galaxy/stellar spectra) to simulate individual images and a final stacked image. The visualization include:

- **single & stacked image:** Presented in the upper left and right corners of this [Figure](https://www.notion.so/A-generic-spectrograph-model-for-exposure-time-calculations-trade-studies-and-instrument-observati-10493411461380889a97d23f7e382d14?pvs=21), these images are 100 × 500 pixels (resulting in distinct physical FOVs for different instruments). The spectral direction is horizontal, while the spatial one is vertical. The slit size is incorporated (in arc-second and kpc if the redshift is informed in the Imaging ), along with contributions from different noise sources.
- **Wavelength dependant informations**: Evolution of atmospheric transmission, throughput and sky/background emission with wavelength (see previous Section).
- **Histogram:** Lower left plot displays the histogram for both individual and stacked images.
- **Profiles:** Lower right plot offers profiles in both spatial and spectral directions for the single (large & transparent) and stacked images.

![Observation of a single slit/fiber spectra (single and stacked image on top respectively). The bottom left plots show respectively the λ-dependant characteristics (possible sky emission lines, atmospheric transmission and instrument throughput)](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/d7b5e3d8-8dd4-491c-bb06-d0b8b614468e/Capture_decran_le_2024-09-23_a_17.30.32.jpg)

Observation of a single slit/fiber spectra (single and stacked image on top respectively). The bottom left plots show respectively the λ-dependant characteristics (possible sky emission lines, atmospheric transmission and instrument throughput)

### Cube simulator

For integral field spectrographs (dimensions=3 in the spreadsheet, slicer or fiber IFU), the code uses the generated spectra (presented in the previous subsection) to generate the spatio-spectral data cube. For sake of execution speed, the code only converts the 2D spectra into a 3D cube by convolving it by the spatial extension of the source (added quadratically to the instrument resolution). Source signal and background (sky, dark, addition background) are added independantly and readnoise pixels values are shuffled to prevent noise correlation.

![Cube and spectra visualization of IFS instruments. This output tab is only available when IFS is selected. ](attachment:ea8bf951-2089-4783-ac48-ace35ffefdc6:Capture_decran_le_2025-02-12_a_16.25.06.jpg)

Cube and spectra visualization of IFS instruments. This output tab is only available when IFS is selected. 

# Deployment

## Instrument data base

To facilitate the addition of new instruments, a database has been developed on a Google spreadsheet. This database organizes all instrument and observation parameters, with descriptions of each parameter accessible by clicking on the respective cell. Units are listed in the third column, and users can add remarks or comments to notify others of changes or important considerations. This method ensures essential versioning, tracks changes, and allows for easy user interaction via comments.

Different instrument or source configurations are also possible to add by just creating a dropdown button and use the options of this button to modify the other parameter cells accordingly. This allows users to efficiently modify instrument subsystems or observation strategies, such as:

- **Spectro design mode:** to analyze for instance the impact of using the low, medium or high resolution gratings (with the possibility to also account for the narrow of wide field of view).
- **Detector type:** Allows to simply examine the SNR evolution with different detector setup (EMCCD, CMOS, CCD, qCCD, MCP) which can automatically change the detector parameters.
- **Observation strategy:** The user can also define different science case by changing the source parameter based on its name.

For example, users can implement the following formula:

```jsx
Slitlength=if(G11="Slicer",120,4)
Spectral_resolution=int(REGEXEXTRACT (R11, "\d+"))*if(R8=0.18,1,if(R8=0.09,2,4))
```

![ Instrument database on Google sheet. Any user can easily enter his instrument configuration in an “in-progress” tab. Once done, it will be added to this main page and will be directly accessible from the ETC.](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/3b8bf5e1-315f-4a09-a149-7e3ac0ccbf50/Capture_decran_le_2024-09-12_a_16.51.08.jpg)

 Instrument database on Google sheet. Any user can easily enter his instrument configuration in an “in-progress” tab. Once done, it will be added to this main page and will be directly accessible from the ETC.

Once completed, configurations can be accessed directly by the ETC. Additionally, the Google spreadsheet includes a tab that enables users to compare instrument performance. Selecting a performance metric or instrument for comparison generates a radar plot of the most generic performance metrics.

![Instrument performance comparison plot. **Left:** Control panel for the user to change parameters, instrument or scale options. **Middle:** Scatter plot to compare all instrument in some parameter space. **Right:** Radar plot that encompasses the most generic performance metric for up to 5 instruments.](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/94590f13-1d42-4d9d-aa19-0e7a87debaa6/Capture_decran_le_2024-09-23_a_17.39.24.jpg)

Instrument performance comparison plot. **Left:** Control panel for the user to change parameters, instrument or scale options. **Middle:** Scatter plot to compare all instrument in some parameter space. **Right:** Radar plot that encompasses the most generic performance metric for up to 5 instruments.

Groups or scientist who work on designs, or instruments characteristics they do not wish to share, can use the local instrument database instead of the online one (simple widget on the ETC). The local `.csv` database, in the GitHub repository can be modified locally, allowing to perform private work or to work without an internet connection.

## Notebook and python code

One significant interest of the ETC is its simplicity of use. It is as simple as a Jupyter notebook. We have designed the algorithm of the ETC and have written the associated code in Python with a few basic external package dependencies (the use of only well-maintained packages such as `Astropy` or `NumPy`  allows to prevent rapid obsolescence). The ETC is built using Python 3.12 and features a graphical user interface (GUI) that allows for cross-platform use. Accordingly, the ETC is able to run on Windows, Linux, and Mac OS systems both locally (usually much faster) and virtually thro–––ugh any web browser. The use of ipywidget allows an extremely easy interaction through sliders and menus, organized into tabs.

The source code is made available on GitHub for public use, the code is easy to maintain with or without collaboration. The GitHub [repository](https://github.com/vpicouet/spectro-imager-etc) is not only used for version tracking/management bus also for automatic code deployment.

The testing was conducted using a Mac laptop with an M2 processor. As a result, the average processing time is less than 5 second to initiate the ETC and about one second to change any parameter or view. 

## Local or server use through Heroku+Voila

The ETC is easy to use locally on VS-code. Using it with *Zen* mode, makes it almost seem as a desktop application. Alternatively, it can directly be used on the web without having to install anything. The application is automatically deployed on Heroku server accessible [here](https://spectro-imager-etc-fbcad4e5849d.herokuapp.com), each time the code is modified and pushed through GitHub. Heroku’s integration with Voila allows a seamless transition from code to an interactive dashboard, simplifying the user experience. Additionally, self-hosting options are available for those who prefer to manage their own servers, ensuring that the ETC can be used conveniently by any team, regardless of their technical infrastructure.

![ Diagram of how the database, notebook, server and user interact.](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/ae247bd0-5ea1-4ae1-b472-836c3b746f61/image.png)

 Diagram of how the database, notebook, server and user interact.

# SNR computation

## Contribution

### Source’s and background/sky signal

The sky and signal contributions are first converted from $ergs/cm^2/s/asec^2/Å$ to continuum unit (photons/cm $^2$/s/sr/Å) :

$F_{CU} =   \frac{F_{ergs/cm^2/s/asec^2/Å}}{\frac{h c}{ \lambda} \times \frac{\pi}{ 180 \times 3600}^2 }$

The spatial shape is modelled by a gaussian distribution (converging towards a punctual source or a flat surface brightness source for respectively small and large FWHM)

The average spectral flux will be based on the user defined flux. The spectral shape though, will depend on the user's choice:

1. Emission line modelled by gaussian distribution (converging towards a sharp emission line or a flat continuum for respectively small and equivalent width). We decided deliberately to use flux per Angstrom with a gaussian profile so that the user can generate realistic simulations of both continuum and emission line sources.
2. Rest frame or observed frame spectra
3. Rest frame blackbody spectra

Multiple factors influence the throughput determination, which include light transmission through the telescope structure, the spectrograph, and the atmosphere. The total throughput is given by:

$Throughput(\lambda)=Atm_{\%}(\lambda)\times T_{\%}(\lambda)  \times QE_{\%}(\lambda)$

Then, both contributions are converted similarly into electrons per pixels:

$F_{e-/pix/exp} = S_{CU} \times Sw_{str}  \times d_{Å/pix}  \times t_{exp_{s}}  \times  Area_{T} \times Throughput$  

$Sky_{e-/pix/exp} = S_{CU} \times min(Sw_{str},σx_{str})  \times d_{Å/pix} \times texp_{s} \times  Area_{T} \times Throughput$

With $S_{CU}$ the source flux in continuum unit,$Sw_{str}$ the slit width in steradian, d the dispersion.

We also account for the flux cut by the spectrograph input (slit/fiber/slice) based on the PSF size there versus the size of the slit/fiber. For fiber systems, an additional π/4 correction is applied to account for flux loss at the fiber edges.

### **Other regular noise contributions**

Other contributions (dark current, read-noise, CIC) are straighforward to account for.

- Dark and extra-background sources (eg. straylight, IR thermal emission, light leak) are used with the same unit ($e^-$/pix/hour). Therefore $D_{e^-/pix/exp} = D_{e^-/pix/hour} \times \frac{ texp_{s} }{3600}$
- Read noise (RN) is directly given in $e^-$/pix/exp. It is simulated as a normal distribution centered on 0 and is applied after all sources of noise have been added.

## Noise model / Conversion to noise

This part describes the method used by the ETC for calculating the signal to noise ratio.

### Regular detector model

Each contribution is then converted to noise by taking the square root of the contribution and accounting for the effective number of frames and the element resolution size.

The number of effective images is affected by both the exposure time and the readout time, as well as the image-loss fraction due to cosmic rays:

$N_{images} = \frac{Ttot_{s}}{texp_{s} + tro_{s}} \times (1-CR_{\%})$  

The evolution of the SNR 

$SNR(t_{exp}) \frac{F \times t_{exp}}{\sqrt{ (\sigma^2_F + \sigma^2_S + \sigma^2_D+\sigma^2_B) \times t_{exp} + \sigma_{RN}^2  }}   =\frac{F \times t_{exp}}{\sqrt{ (F + S+D+B) \times t_{exp} + \sigma_{RN}^2  }}$ 

### EMCCD noise model

In the case of electron-amplified CCDs, some considerations must be taken into account:

- Clock induced charges (CIC) are charges induced in electron amplified CCD (already given in e-/pix/exp)
- All electrons before readout (flux, dark, CIC) are going to be amplified by the amplification gain.
- The read noise is not amplified and its impact must then be divided by the amplification gain: $RN_{e-/pix/exp} = \frac{ RN_{e-/pix/exp} }{G_{e-/e-}}$
- An excess noise factor of $\sqrt{2}$ must be used to account for the stochastic amplification if no thresholding algorithm (thresholding or Bayesian) is applied

We use the model detailed in \citet{Harpsoe} which leads to this SNR without thresholding:

$SNR(t_{exp}) = \frac{G \times F \times t_{exp}}{\sqrt{(\sigma^2_F + \sigma^2_S + \sigma^2_D+\sigma^2_B) \times t_{exp}  + \sigma_{CIC}^2 + \sigma_{RN}^2}} = \frac{G \times F \times t_{exp}}{\sqrt{ENF^2 \times G^2 \times [(F + S + B + D) \times t_{exp}  + CIC] + \sigma_{RN}^2}}$

The excess noise factor can be partly suppressed when some photon-counting algorithm like thresholding or Bayesian Inference is applied.

## Design

The SNR calculation is deterministic so that it always produces the same output, given identical input. The image and cube simulation draw pixels based on stochastic processes (poisson, gamma, normal distributions) and therefore are not deterministic.

The SNR calculation is designed to generate an array $SNR(P)=\frac{S(P)}{N(P)}$ value based on any instrumental observation parameter (P) which easily allows to visualize the evolution of S(P), N(P) and SNR(P) with any given parameter. The next [Figure](https://www.notion.so/A-generic-spectrograph-model-for-exposure-time-calculations-trade-studies-and-instrument-observati-10493411461380889a97d23f7e382d14?pvs=21) displays the detailed flowchart for this mode.

At the initialization phase, the ETC imports and reads the Google sheet data base.

![Diagram of the ETC’s design. ](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/a9c95789-c214-42e5-9ec9-bc2c59ed1155/image.png)

Diagram of the ETC’s design. 

The ETC also includes an image/cube simulation code that calculates the two or three dimensional pixel-by-pixel signal and noise properties of any of the selected instrument. A detailed explanation of this simulation is provided in the annex.

# Adapting the ETC for more specific trades analysis : FIREBall-2 case

| **Trade** | **Low limiter** | **High limiter** | **Optimum** | **Dev** |
| --- | --- | --- | --- | --- |
| **Exposure time** | Readout time, CIC and read noise 
(& thresholding efficiency) | Cosmic rays | ~60 seconds depending on the CR rate |  |
| **Slit size** | Loss of flux | Sky background | ~70 mu, depending on source size |  |
| **EMCCD temperature** | Charge transfer efficiency | Dark current | -90, -100 C | ✓ |
| **Photon counting threshold** | Read noise | Reduces flux | ~5σ, or more in the presence of CTI  | ✓ |

## Native/simple tradeoff analysis

### Exposure time

This tradeoff results from an important read noise contribution at short exposure times, important surface loss due to the ~10/sec CR impact rate for long exposures, and the photon counting thresholding efficiency that decreases significantly when the number of electrons per pixel per exposure gets too high/low.

### Slit size

This tradeoff result from the substantial flux loss when the slit width or length get significantly smaller than the PSF size at the slit level. Conversely, the source’s flux get negligeable compared to the sky or background when the slit becomes too large (in addition of reducing rapidly the spectral resolution)

## Specific tradeoff analysis for FIREBall

To meet FIREBall-2's specific needs, we added several hidden features to the ETC to analyze more nuanced SNR trade-offs. These additional tools can serve as examples for other projects that require exploration of parameter spaces beyond the default instrument database.

Indeed, most of the parameters in the database are fixed for the observed or even for the instrument scientist (noise contributions of the detector, slit/fiber size, etc). But in some case, these parameters are driven by some other parameters that could be controlled (eg. the detector temperature that will drive the dark current and the charge transfer efficiency, the detector controller frequency that could impact the detector readout noise or clock induced charges).  If such correlations are known, the user could implement them (Dark(T), CTE(T)) to directly observed the evolution of the SNR with temperature. This is what was performed for FIREball-2 for the temperature and photon counting threshold optimization.

### Detector temperature

While dark current increases exponentially with temperature, charge transfer efficiency decreases as the temperature drops (below ~190 K). Parameterizing both the dark current and CTE as functions of temperature allows for determining the optimal temperature trade-off between minimizing dark current and maximizing CTE.

### Photon counting thresholding

Thresholding in EMCCDs helps reduce the excess noise factor but depends heavily on certain conditions. A low threshold will result in an increased false event rate, while a high threshold leads to significant flux loss, negating the benefits of thresholding. The optimal SNR(threshold) is influenced by the amplification gain (which is affected by smearing) and the ratio of read noise to input flux. This ratio must remain significantly greater than 10 for efficient photon counting but is also impacted by other noise sources.

# Discussion and conclusion

## Work in progress

**Addition of simulated data cubes :**

An improvement of the work is the possibility to directly add data simulated data cubes to the 

## Limitations and possible upgrades

The goal of this ETC is not to compete with instrument-specific ETCs, which are often far more precise. Here, we outline several limitations:

**Limitations:**

- Limited parameters: The ETC does not account for field distortion, PSF evolution, adaptive optics, non-Gaussian PSF shapes, correlated detector noise, or read noise floors. Only parameters entered in the Google Sheet are used.
- No diffraction effects are included in the model.
- Limited validation: The performance predictions have only been validated on a couple instruments and therefore nothing can be guaranteed

**Possible upgrade;**

- Direct spectral uploads: Users will soon have the possibility to directly upload spectra in $ergs/cm^2/s/asec^2/Å$ in the GitHub repository (under notebooks/Spectra, λ in nanometers on the first column and with no column name).

## Conclusion

In this paper, we have described the development of a generic tool at the frontier between an exposure time calculator and an instrument model. Combining these two strenghs offer some possibilities that are currently not available on current ETCs (instruments specific or not). Its flexibility makes it easily adaptable to multiple instruments and configurations, allowing astronomers to quickly estimate signal-to-noise ratios (SNR), run trade-off studies, and make informed decisions on observation strategies. Its user-friendly interface, rapid computation times, and the ability to add custom instruments make it a practical tool to better optimize future instrument designs. It mostly serves as an example to potentially develop other kind of ETCs that would not only focus on observations but also on arbitrating trades analysis for instrument design.

Due to its simplicity, other projects could easily fork or reproduce this ETC to adapt it to their own need or specificity. Such an example was presented in the last section where different trades analysis were presented, some of them necessitating minor code changes to explore a different parameter space. 

# Acknowledgements

This work was supported by the Centre National d'Etudes Spatiales (CNES) and by the NASA Astrophysics Research and Analysis Program (APRA) under grant 80NSSC20K0395. We are grateful for NASA’s continued support in the development of tools that enable innovative astronomical research.

We would like to acknowledge the developers of `Matplotlib`, `ipywidgets`, and `Astropy`, whose open-source software libraries were instrumental in the development of the Exposure Time Calculator.

Special thanks go to Bruno Milliard and David Schiminovich for their valuable feedback and insights, which greatly improved the design and functionality of this work.

# References

1. K. M. Pontoppidan, T. E. Pickering, V. G. Laidler, *et al.*, “Pandeia: a multi-mission exposure
time calculator for JWST and WFIRST,” in *Observatory Operations: Strategies, Processes,
and Systems VI*, A. B. Peck, R. L. Seaman, and C. R. Benn, Eds., *SPIE* **9910**, 991016 (2016).
[doi:10.1117/12.2231768].
2. L. D. Nielsen, P. Ferruit, G. Giardino, *et al.*, “The JWST/NIRSpec exoplanet exposure time
calculator,” in *Space Telescopes and Instrumentation 2016: Optical, Infrared, and Millimeter
Wave*, H. A. MacEwen, G. G. Fazio, M. Lystrup, *et al.*, Eds., *SPIE* **9904**, 99043O (2016).
[doi:10.1117/12.2231624].
3. Kremic T. (2015) [Stratospheric Balloons for Planetary Science and the Balloon Observation Platform for Planetary Science](https://ntrs.nasa.gov/api/citations/20150010715/downloads/20150010715.pdf)
4. Gill A., et al. (2020) Optical night sky brightness measurements from the stratosphere, AJ

# Appendices

## Full code design

![Flow chart representation of the full ETC and image/cube simulator. The flowchart shows the functions and parameters for S/N calculation mode.](https://prod-files-secure.s3.us-west-2.amazonaws.com/134ff447-34fb-48c7-aa71-0175c4110b5d/6c77f39f-5bcf-4fdd-9bbd-5147ee768d89/image.png)

Flow chart representation of the full ETC and image/cube simulator. The flowchart shows the functions and parameters for S/N calculation mode.

## Full list of parameters

- **Source:**
    - **Flux** (F): $erg/cm^2/s/arcsec^2/Å$
    - **Sky** (S): Level of sky background illumination in $ergs/cm^2/s/arcsec^2/Å$
    - **Source extension** (σ$_x$): Spatial extension of the source in arcseconds
    - **Source's line width** (σ$_λ$): Spectral extension of the source/emission line in Å
- **Observing strategy:**
    - **Observed wavelength** (λ): Effective wavelength of the instrument in nanometer
    - **Exposure time** ($t_{exp}$): Exposure time of individual frames in seconds
    - **Total acquisition time** ($Ttot$): Total acquisition time in hours
    - **Atmospheric transmission** (Atm$_{\%}$): %/100
    - **Distance to source/line center** ($Δ_x, Δ_y$): Distance to the source being analyzed ['' or Å]
- **Instrument design:**
    - **Collecting area** (A): Physical collecting area of the instrument in m$^2$
    - **Plate scale** (P): Pixel plate scale in ''/pix
    - **Throughput** (T$_{\%}$): Instrument throughput at effective wavelength (not accounting for detector quantum efficiency and atmospheric transmission)
    - **Spatial resolution** (R$_{mask}$*,R$_{det}$*): Respectively at the mask (Encompassess seing / guiding / Telescope PSF) and at the detector
- **Spectrograph design:**
    - **Spectral resolution** (R$_λ$): Spectrograph spectral resolution λ[Å]/dλ[Å]
    - **Slit width** (Sw): Width of the slit/slice/fiber in ''
    - **Slit length** (Sl): Length of the slit/slice/fiber in ''
    - **Dispersion** (d): Dispersion at the detector Å/pix
    - **IFS**: Convert the slit or fiber design into a slicer or fiber-IFS. This allows a 3D visualization. The impact on the noise is that the flux (sky or source) that is usually cut by a slit or single fiber is then transfered to the neighboring fibers/slices suppressing this loss factor.
- **Detector parameters:**
    - **Quantum efficiency** (QE$_\%$): Detector quantum efficiency in %/100
    - **Dark current** (D): Detector dark current [e-/pix/hour]
    - **Read noise** (RN): Detector readout noise in electrons/pixel
    - **Readout time** (t$_{RO}$): Time in seconds that it takes for images to get read. Use 0 for MCPs or rolling shutter
    - **Pixel size** (P$_s$): Pixel size in microns
    - **Image loss due to cosmic ray** (CR$_\%$): Cosmic ray loss per second. eg. 0.01 would mean that 1 sec image looses 1% pixels due to cosmic rays
    - **Extra background** (B): Any additional source of straylight (could be any internal or external contribution/leak), in addition to the nominal sky background. (e-/pix/hour)
- **emCCD additional parameters:**
    - **EM gain** (G): EMCCD amplification gain in e-/e-
    - **CIC** (CIC): EMCCD spurious charges due to amplification in electrons [e-/pix]
    - **Smearing exponential length** (CTE): Smearing length of the EMCCD (exponential length in pixels). This length, representing the charge transfer efficiency is fixed by the temperature when the Temp checkbox is checked.
    - **Thresholding**: [Yes/No] and the possibility to choose the threshold value in the presence of smearing (which reduces photon counting efficiency)
    - **Temperature** (T): if you check it (based on a first rough evolution of smearing and dark current with temperature, therefore changing the temperature will change smearing and dark accordingly.)
- **Image simulator-related parameters:** Full well of the detector
    - **Conversion gain** (c$_g$): to convert e- into detector ADU (ADU/e-)
    - **Throughput FWHM** (FWHM$_{Th}$): taking into account all optics and QE, not atmosphere, to add the λ dependency.
    - **Atmospheric transmission** [Yes/No] to add a λ-depend transmission model (based on [pwv_kpno](https://mwvgroup.github.io/pwv_kpno/1.0.0/documentation/html/atmospheric_modeling.html)) for ground instruments (this only applies to the source, not to the sky emission)
    - **Atmospheric emission lines**: [Yes/No] replace a flat sky continuum by sky emission lines based on [UVES estimates](https://www.eso.org/observing/dfo/quality/UVES/pipeline/sky_spectrum.html)