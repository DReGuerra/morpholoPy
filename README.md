# wrinkleAnalysis
Fast Fourier Transform analysis for nanowrinkled surface feature analysis

André Guerra \
April, 2024 \
andre.guerra@mail.mcgill.ca  

---
Description: \
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

Source: https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf<br>
Figure 1<br>

![Stimpson et al 2020 Figure 1](./stimpsonetal2020_fig1.png)

#### Vertical lines test results:
![Vertical lines](tests/vertical_lines/figures/verticalLines_summary.png)

### chevron

This test reproduces the results from Stimpson et al., 2020 for the chevron lines in Figure 1A second panel.<br>

Source: https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf<br>
Figure 1<br>

##### Chevron test results:
![Chevron](tests/chevron/figures/chevron_summary.png)

### jigsaw and fragmented_jigsaw

This test reproduces the results from Stimpson et al., 2020 for the jigsaw and fragmented_jigsaw lines in Figure 1A third and forth panels.<br>

Source: https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf<br>
Figure 1<br>

The jigsaw and fragmented jigsaw images were made using Photopea at https://www.photopea.com/. The original vertical lines image was used to create these transformed images.<br>

#### Jigsaw test results:
![Jigsaw](tests/jigsaw/figures/jigsaw_summary.png)

#### Fragmented jigsaw test results:
![Fragmented jigsaw](tests/fragmented_jigsaw/figures/frag_jigsaw_summary.png)

### kdf_biaxial_20um

This test uses a biaxial wrinkled surface with a scale bar of 20 microns (um). The scale bar is 170 pixels in length.<br>

Source: https://pubs.acs.org/doi/full/10.1021/acsami.8b16232<br>
ACS Appl. Mater. Interfaces 2019, 11, 6, 6325–6335<br>
Figure 6F

![KDF biaxial wrinkles](tests/kdf_biaxial_20um/figures/kdf_biaxial_20um_summary.png)

### kdf_uniaxial_20um

This test uses a uniaxial wrinkled surface with a scale bar of 20 microns (um). The scale bar is 170 pixels in length.<br>

Source: https://pubs.acs.org/doi/full/10.1021/acsami.8b16232<br>
ACS Appl. Mater. Interfaces 2019, 11, 6, 6325–6335<br>
Figure 6J

![KDF uniaxial wrinkles](tests/kdf_uniaxial_20um/figures/kdf_uniaxial_20um_summary.png)