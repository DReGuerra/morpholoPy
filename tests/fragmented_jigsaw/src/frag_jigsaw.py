"""
    André Guerra
    andre.guerra@mail.mcgill.ca

    This script reproduces the results from Stimpson et al., 2020 for the chevron pattern lines in Figure 1
    Ref: https://chemrxiv.org/engage/chemrxiv/article-details/60c74e50f96a009895287acf
"""

# import pandas as pd
# from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from skimage import feature, io, util
from _wrinklelib import radialProfile, radAvgPSD

# aray size NxN
N = 900             # pixels
# scale
L = N               # um
x_scale = L/N       # um/pixel
# create NxN array
img = np.zeros((N,N))
img2 = np.zeros((N,N))

w = 50
for i in np.arange(0,N,w*2):
    img[:, i:i+w] = 255
    img2[:, i:i+w] = 255

# make chevron image from the vertical lines image
period = N/20
r = 1
for i in np.arange(1,N-1,1):
    if r > 0:
        tmp = img2[i-1, :]
        tmp = np.delete(tmp, -1)
        tmp = np.insert(tmp, 0, 0)
        img2[i, :] = tmp
        if i%period == 0: r = r*-1
    else:
        tmp = img2[i-1, :]
        tmp = np.delete(tmp, 0)
        tmp = np.insert(tmp, -1, 0)
        img2[i, :] = tmp
        if i%period == 0: r = r*-1
        
# import fragmented jigsaw image in gray-scale
image = io.imread("images/frag_jigsaw.png", as_gray=True)

# Dicrete Fourier Transfer via FFT
# 2D DFT
fft2 = np.fft.fft2(image)
# center-shifted 2D DFT
fft2_shiftd = np.fft.fftshift(fft2)
# power spectral density (PSD)
psd2D = np.abs(fft2_shiftd)**2

thetaLims = [170 - 180, 190 - 180]
# matlab translation
rasp, binsCount = radAvgPSD(psd2D, thetaLims)
# length of rasp vector
raspLen = len(rasp)

# frequency vector (pixels)
k = np.arange(0,N-1,1)          # pixels
# spatial frequency vector (lam)
lam = np.divide(L,k)            # un/pixel
# normalize rasp with binCount
rasp_norm = np.nan_to_num(np.divide(rasp,binsCount))
# deconvolve
dconv = True
if dconv: rasp_norm = np.divide(rasp_norm,lam[:raspLen])
# transform to absolute value using maximum intensity
rasp_norm_au = rasp_norm/np.max(rasp_norm)

# use radialProfile()
rasp_norm2 = radialProfile(psd2D)
rasp_norm2_au = rasp_norm2/np.max(rasp_norm2)   # AU

# Visualization constants
# sizing factors
TICKSFONT = 13; TITLEFONT = 15; TEXTFONT = 15   # fonts
FIGWIDTH = 6.4; FIGHEIGHT = 4.8                 # figure size
LINEWIDTH = 3; ROLLWINDOW = 100                 # plot spec

# summary figure
NROWS = 3; NCOLS = 2
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))
axs[0,0].imshow(image, cmap='gray')
axs[0,0].set_title("Original image")

axs[0,1].imshow(psd2D)
axs[0,1].set_title("2D PSD")

axs[1,0].plot(rasp_norm_au)
axs[1,0].set_title("radAvgPSD() - MATLAB translation")
axs[1,0].set_ylabel("Intensity, AU")
axs[1,0].set_xlabel("Period, Pixels")
axs[1,0].set_xlim([0,50])

axs[1,1].plot(rasp_norm2_au)
axs[1,1].set_title("radialProfile()")
axs[1,1].set_ylabel("Intensity, AU")
axs[1,1].set_xlabel("Period, Pixels")
axs[1,1].set_xlim([0,50])

axs[2,0].plot(1/lam[:raspLen],rasp_norm_au)
axs[2,0].set_title("radAvgPSD() - MATLAB translation - deconvolved")
axs[2,0].set_ylabel("Intensity, AU")
axs[2,0].set_xlabel("Period, 1/$\mu$m")
axs[2,0].set_xlim([0,0.08])

# axs[2,1].axis('off')
axs[2,1].plot(1/lam[:raspLen],rasp_norm/np.max(rasp_norm))
axs[2,1].set_title("Not deconvolved")
axs[2,1].set_ylabel("Intensity, AU")
axs[2,1].set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
axs[2,1].set_xlim([0,0.08])

f.tight_layout()
f.savefig("figures/frag_jigsaw_summary.png")

# summary figure
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.plot(1/lam[:raspLen],rasp_norm/np.max(rasp_norm))
axs.set_ylabel("Intensity, AU")
axs.set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
axs.set_xlim([0,0.08])
f.tight_layout()
f.savefig("figures/frag_jigsaw_psd.png")