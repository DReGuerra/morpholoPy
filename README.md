# wrinkleAnalysis
Fast Fourier Transform analysis for nanowrinkled surface feature analysis

André Guerra \
April, 2024 \
andre.guerra@mail.mcgill.ca  

---
<b>Description:</b> \
This repository contains a python implementation of image analysis algorithms to examine nanostructure on material surfaces. The scripts take in image files (.png), conduct the analyses and output summary figures.

---
## Core Contents
1. `tests/` $\rightarrow$ collection of test cases to validate and demonstrate usage.
2. `_wrinklelib.py` $\rightarrow$ contains functions to be used by `wrinkleAnalysis.py`
3. `wrinkleAnalysis.py` $\rightarrow$ main script executing the image analyses and surface feature examinations.

## References
1. [Stimpson et al., 2020](https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf)
2. [De France et al., 2019](https://pubs.acs.org/doi/full/10.1021/acsami.8b16232)

---

## Tests

### vertical_lines

This test reproduces the results from Stimpson et al., 2020 for the vertical lines in Figure 1A first panel.<br>

Source: [Stimpson et al., 2020](https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf)<br>
Figure 1<br>

![Stimpson et al 2020 Figure 1](./stimpsonetal2020_fig1.png)

#### Vertical lines test results:
![Vertical lines](tests/vertical_lines/figures/verticalLines_summary.png)

### chevron

This test reproduces the results from Stimpson et al., 2020 for the chevron lines in Figure 1A second panel.<br>

Source: [Stimpson et al., 2020](https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf)<br>
Figure 1<br>

##### Chevron test results:
![Chevron](tests/chevron/figures/chevron_summary.png)

### jigsaw and fragmented_jigsaw

This test reproduces the results from Stimpson et al., 2020 for the jigsaw and fragmented_jigsaw lines in Figure 1A third and forth panels.<br>

Source: [Stimpson et al., 2020](https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf)<br>
Figure 1<br>

The jigsaw and fragmented jigsaw images were made using Photopea at https://www.photopea.com/. The original vertical lines image was used to create these transformed images.<br>

#### Jigsaw test results:
![Jigsaw](tests/jigsaw/figures/jigsaw_summary.png)

#### Fragmented jigsaw test results:
![Fragmented jigsaw](tests/fragmented_jigsaw/figures/frag_jigsaw_summary.png)

### kdf_biaxial_20um

This test uses a biaxial wrinkled surface with a scale bar of 20 microns (um). The scale bar is 170 pixels in length.<br>

Source: [De France et al., 2019](https://pubs.acs.org/doi/full/10.1021/acsami.8b16232)<br>
ACS Appl. Mater. Interfaces 2019, 11, 6, 6325–6335<br>
Figure 6F

![KDF biaxial wrinkles](tests/kdf_biaxial_20um/figures/kdf_biaxial_20um_summary.png)

### kdf_uniaxial_20um

This test uses a uniaxial wrinkled surface with a scale bar of 20 microns (um). The scale bar is 170 pixels in length.<br>

Source: [De France et al., 2019](https://pubs.acs.org/doi/full/10.1021/acsami.8b16232)<br>
ACS Appl. Mater. Interfaces 2019, 11, 6, 6325–6335<br>
Figure 6J

![KDF uniaxial wrinkles](tests/kdf_uniaxial_20um/figures/kdf_uniaxial_20um_summary.png)

## Usage

### `wrinkleAnalysis.py`

The workflow of this script is as follows:<br>
1. Input parameters - to be manually changed to reflect the desired analysis, define the image file, and size parameters which depend on the scale of the image. The size parameters can be determined through trial and error by examining the curve fitting and wrinkle size estimate in the panel 6 of the summary figure produced.
2. Input image file (`imgFile`)
3. Band-pass filtering of the image using Gaussian differences
4. Canny edge detection - `sigma` value is modified to achieve the desired granularity in edge detection
5. Square the image (if not already squared)
6. Determine image scale and pixel relationship - if the image contains an SEM banner with a size scale, the function `measureScaleBar()` can be used to extract the pixel length of the scale bar. The representative physical length of the scale bar must be input manually (`barLength`).
7. Center-shifted 2D FFT (`fft2_shiftd`) of the image to obtain the 2D power spectral density (PSD) (`psd2D`) of the image.
8. Radially averaged 2D PSD (`radPSD`).
9. Frenquency conversion from pixel to spatial length.
10. High and low filters.
11. Curve fitting of the most prominent peak in the `radPSD` to extract the feature size.
12. Visualization of the results output to `figures/`.