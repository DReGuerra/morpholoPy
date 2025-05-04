"""
    André Guerra
    andre.guerra@mail.mcgill.ca

    This script reproduces the results from Stimpson et al., 2020 fir the vertical lines in Fig. 1
    Ref: https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf
"""

import numpy as np
from scipy.stats import kurtosis
import matplotlib.pyplot as plt
from surfacetools.periodicfeatures import radially_averaged_PSD

# Generate an image of alternating black and white vertical stripes
# aray size NxN
N = 900             # pixels
# scale
L = N               # um
x_scale = L/N       # um/pixel
# create NxN array
img = np.zeros((N,N))
# width of one stripe
w = 50              # pixels
for i in np.arange(0,N,w*2):
    img[:, i:i+w] = 255

# Dicrete Fourier Transfer (DFT) via FAst Fourier Transfrom (FFT)
# 2D DFT
fft2 = np.fft.fft2(img)
# center-shifted 2D DFT
fft2_shiftd = np.fft.fftshift(fft2)
# power spectral density (PSD)
psd2D = np.abs(fft2_shiftd)**2
# angle limits for radial averaging the 2D PSD
theta_lims = [170 - 180, 190 - 180]
# radially averaged PSD with angle limits
rasp, bins_count = radially_averaged_PSD(psd2D, theta_lims)
# length of rasp vector
rasp_length = len(rasp)

# frequency vector (pixels)
k = np.arange(0,N-1,1)          # pixels
# spatial frequency vector (lam)
lam = np.divide(L,k)            # um/pixel
# normalize rasp with bins_count
rasp_norm = np.nan_to_num(np.divide(rasp,bins_count))
# deconvolution to spatial period
dconv = True
if dconv: rasp_norm = np.divide(rasp_norm,lam[:rasp_length])
# transform to absolute units (AU) using maximum intensity
rasp_norm_au = rasp_norm/np.max(rasp_norm)

# kurtosis
kurt = kurtosis(img.flatten(), fisher=False)
print("kurtosis = " + str(kurt))

# ---
# Visualization

# format sizing constants
TICKSFONT = 13; TITLEFONT = 15; TEXTFONT = 15   # fonts
FIGWIDTH = 6.4; FIGHEIGHT = 4.8                 # figure size
LINEWIDTH = 3; ROLLWINDOW = 100                 # plot spec

# 1. summary figure
NROWS = 2; NCOLS = 2
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))
axs[0,0].imshow(img, cmap='gray') # original image
axs[0,0].set_title("(a)", loc='left', fontsize=TITLEFONT)
axs[0,0].axis("off")

axs[0,1].imshow(np.log(1 + psd2D)) # log(1 + PSD)
axs[0,1].set_title("(b)", loc='left', fontsize=TITLEFONT)
axs[0,1].axis("off")

axs[1,0].plot(rasp_norm_au) # radially averaged PSD
axs[1,0].set_title("(c)", loc='left', fontsize=TITLEFONT)
axs[1,0].set_ylabel("Intensity, AU", fontsize=TITLEFONT)
axs[1,0].set_xlabel("Period (T), Pixels", fontsize=TITLEFONT)
axs[1,0].set_xlim([0,50])
axs[1,0].tick_params(axis='both', which='major', labelsize=TICKSFONT)

axs[1,1].plot(1/lam[:rasp_length],rasp_norm_au) # radially averaged PSD in spatial frequency
axs[1,1].set_title("(d)", loc='left', fontsize=TITLEFONT)
axs[1,1].set_ylabel("Intensity, AU", fontsize=TITLEFONT)
axs[1,1].set_xlabel("frequency (f), $\mu$m$^{-1}$", fontsize=TITLEFONT)
axs[1,1].set_xlim([0,0.08])
axs[1,1].tick_params(axis='both', which='major', labelsize=TICKSFONT)

f.tight_layout()
f.savefig("figures/verticalLines_summary.png")

# 2. Original image
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))
axs.imshow(img, cmap='gray')
axs.axis("off")
f.tight_layout()
f.savefig("figures/vertical_lines_out.png")

# 3. Radially averaged PSD
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(1/lam[:rasp_length],rasp_norm_au)
axs.set_title("Radially Averaged PSD (RASP)")
axs.set_ylabel("Intensity, AU")
axs.set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
axs.set_xlim([0,0.08])
f.tight_layout()
f.savefig("figures/vertical_lines_rasp.png")