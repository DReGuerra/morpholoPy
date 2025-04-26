"""
Description:
Extract characteristic the length.

    Usage:
        >> python characteristic_length.py file_name dconv bar_len num_populations lo_len_lim_pop1 hi_len_lim_pop1 lo_len_lim_pop2 hi_len_lim_pop2
        
        Paratemetrs:
        file_name: str
            Name of the image file with the extension.
        dconv: bool
            Deconvolution flag.
        bar_len: float
            SEM scale bar length. [um]
        num_populations: int
            Number of populations to analyze.
        lo_len_lim_pop1: float
            Lower limit of the feature size range of interest for population 1 [1/um]
        hi_len_lim_pop1: float
            Upper limit of the feature size range of interest for population 1 [1/um]
        lo_len_lim_pop2: float
            Lower limit of the feature size range of interest for population 2 [1/um]
        hi_len_lim_pop2: float
            Upper limit of the feature size range of interest for population 2 [1/um]

        Returns: 
            Various plots are output to figures/

        EXAMPLES:
        >> python src/characteristic_length.py kdf_biaxial_20um.tif True 1 20 1.5 3.0
        >> python src/characteristic_length.py kdf_biaxial_20um.tif True 2 20 1.5 3.0 3.5 5.0
"""

import sys
import re

import numpy as np

import matplotlib.pyplot as plt
from skimage import feature, io
from skimage.filters import difference_of_gaussians

from surfacetools.image_processing import measure_sem_scalebar
from surfacetools.periodicfeatures import radially_averaged_PSD, full_width_half_max, peak_quality_factor

# inpput parameters
# SEM image file
file = str(sys.argv[1])
file_name = re.findall(r"^(.*)(\.[a-zA-Z]{3})$", file)[0][0]
file_type = re.findall(r"^(.*)(\.[a-zA-Z]{3})$", file)[0][1]
img_file = "images/" + file_name + file_type
# deconvolution
dconv = bool(sys.argv[2])
POP_NUM = int(sys.argv[3])                      # number of populations
# SEM scale bar length
bar_length = float(sys.argv[4])                 # um
# feature size range of interest
# population 1
LO_LEN_LIM_POP1 = float(sys.argv[5])        # 1/um
HI_LEN_LIM_POP1 = float(sys.argv[6])       # 1/um
# population 2
if int(POP_NUM == 2):
    LO_LEN_LIM_POP2 = float(sys.argv[7])        # 1/um
    HI_LEN_LIM_POP2 = float(sys.argv[8])       # 1/um

#############################################################################################
# import SEM image in gray-scale
image = io.imread(img_file, as_gray=True)
# enhance edges by band-pass filtering
filtered_image = difference_of_gaussians(image, low_sigma=1, high_sigma=10)
# Canny edge dectection
filtered_edges = feature.canny(filtered_image, sigma=2)

# print size of the image
print("Image size: " + str(image.shape))
# the smaller dimension of the image
print("Smaller dimension: " + str(np.min(image.shape)))

# resize the image to square
N = 881
image_sq = image[:N,:N]
filtered_image_sq = filtered_image[:N,:N]
filtered_edges_square = filtered_edges[:N,:N]

# image length scale
scale_bar = measure_sem_scalebar(image)   # pixels
# scale_bar = 170 # hardcoded for kdf_biaxial_20um.tif
print("Scale bar: "+ str(scale_bar) + " pixels")
X, Y = filtered_edges_square.shape  # pixels
pxl_scale = bar_length/scale_bar    # um/pixel
L = X*pxl_scale                     # um

# Dicrete Fourier Transfer (DFT) via FAst Fourier Transfrom (FFT)
# 2D DFT
fft2 = np.fft.fft2(filtered_edges_square)
# center-shifted 2D DFT
fft2_shiftd = np.fft.fftshift(fft2)
# power spectral density (PSD)
psd2D = np.abs(fft2_shiftd)**2
# angle limits for radial averaging the 2D PSD
theta_lims = []

# radially averaged PSD with angle limits
rasp, bins_count = radially_averaged_PSD(psd2D, theta_lims=[])
# length of rasp vector
rasp_length = len(rasp)

#############################################################################################
# frequency vector (pixels)
k = np.arange(0,N-1,1)      # pixels
# spatial frequency vector
lam = np.divide(L,k)        # um/pixel
# normalize rasp with bins_count
rasp_norm = np.nan_to_num(np.divide(rasp,bins_count))
# deconvolve
if dconv: rasp_norm = np.divide(rasp_norm,lam[:rasp_length])

# transform to absolute value using maximum intensity
rasp_norm_au = rasp_norm/np.max(rasp_norm)

#############################################################################################
# curve fit
# population 1
ind = np.where(np.logical_and((1/lam>LO_LEN_LIM_POP1),(1/lam<HI_LEN_LIM_POP1)))
x_pop1 = 1/lam[ind]                 # 1/um
y_pop1 = rasp_norm_au[ind]          # AU

deg = 2                             # quadratic poly
z = np.polyfit(x_pop1,y_pop1,deg)   # polynomial coeff
p = np.poly1d(z)
mdl_pop1 = p(x_pop1)
pop1_feature_size = x_pop1[np.where(mdl_pop1 == np.max(mdl_pop1))]  # 1/um

fwhm = full_width_half_max(x_pop1,y_pop1)
print("FWHM = " + str(fwhm))
pqf = peak_quality_factor(y_pop1, fwhm)
print("PQF = " + str(pqf))

# population 2
if int(POP_NUM == 2):
    ind = np.where(np.logical_and((1/lam>LO_LEN_LIM_POP2),(1/lam<HI_LEN_LIM_POP2)))
    x_pop2 = 1/lam[ind]                 # 1/um
    y_pop2 = rasp_norm_au[ind]          # AU

    deg = 2                             # quadratic poly
    z = np.polyfit(x_pop2,y_pop2,deg)   # polynomial coeff
    p = np.poly1d(z)
    mdl_pop2 = p(x_pop2)
    pop2_feature_size = x_pop2[np.where(mdl_pop2 == np.max(mdl_pop2))]  # 1/um

#############################################################################################
# visualization 

# format sizing constants
TICKSFONT = 13; TITLEFONT = 15; TEXTFONT = 15   # fonts
FIGWIDTH = 6.4; FIGHEIGHT = 4.8                 # figure size
LINEWIDTH = 3; ROLLWINDOW = 100                 # plot spec

# summary figure
NROWS = 3; NCOLS = 2
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))
axs[0,0].imshow(image_sq)
axs[0,0].set_title("(a)", loc='left')
# axs[0,0].set_title("Squared image")

axs[0,1].imshow(filtered_image_sq)
axs[0,1].set_title("(b)", loc='left')
# axs[0,1].set_title("Squared filtered image")

axs[1,0].imshow(filtered_edges_square)
axs[1,0].set_title("(c)", loc='left')
# axs[1,0].set_title("Squared Canny edge detected image")

axs[1,1].imshow(np.log(psd2D))
axs[1,1].set_title("(d)", loc='left')
# axs[1,1].set_title("Center-shifted 2D Power Spectral Density")

if int(POP_NUM == 2):
    axs[2,0].plot(1/lam[:rasp_length],rasp_norm_au,linestyle='none',marker='.')
    axs[2,0].plot(x_pop1,y_pop1,linestyle='none',marker='o',fillstyle='none',color='green')
    axs[2,0].plot(x_pop1,mdl_pop1,linestyle='--',color='green')
    axs[2,0].vlines(pop1_feature_size,ymin=0.1,ymax=1,linestyle='--',color='red')
    axs[2,0].set_title("(e)", loc='left')
    axs[2,0].set_ylabel("Intensity, AU")
    axs[2,0].set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
    axs[2,0].annotate('charac. length = ' + str(np.around(1/pop1_feature_size[0],decimals=3)) + ' $\mu$m',
                    xy=(0.45,0.9), xycoords='axes fraction', fontsize=TEXTFONT)
    axs[2,0].set_xlim([0,6])
    axs[2,0].set_ylim([0,1.2])

    axs[2,1].plot(1/lam[:rasp_length],rasp_norm_au,linestyle='none',marker='.')
    axs[2,1].plot(x_pop2,y_pop2,linestyle='none',marker='o',fillstyle='none',color='green')
    axs[2,1].plot(x_pop2,mdl_pop2,linestyle='--',color='green')
    axs[2,1].vlines(pop2_feature_size,ymin=0.1,ymax=1,linestyle='--',color='red')
    axs[2,1].set_title("(f)", loc='left')
    axs[2,1].set_ylabel("Intensity, AU")
    axs[2,1].set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
    axs[2,1].annotate('charac. length = ' + str(np.around(1/pop2_feature_size[0],decimals=3)) + ' $\mu$m',
                    xy=(0.45,0.9), xycoords='axes fraction', fontsize=TEXTFONT)
    axs[2,1].set_xlim([0,6])
    axs[2,1].set_ylim([0,1.2])
else:
    axs[2,0].plot(1/lam[:rasp_length],rasp_norm_au,linestyle='none',marker='.')
    axs[2,0].plot(x_pop1,y_pop1,linestyle='none',marker='o',fillstyle='none',color='green')
    axs[2,0].plot(x_pop1,mdl_pop1,linestyle='--',color='green')
    axs[2,0].vlines(pop1_feature_size,ymin=0.1,ymax=1,linestyle='--',color='red')
    axs[2,0].set_title("(e)", loc='left')
    axs[2,0].set_ylabel("Intensity, AU")
    axs[2,0].set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
    axs[2,0].annotate('charac. length = ' + str(np.around(1/pop1_feature_size[0],decimals=3)) + ' $\mu$m',
                    xy=(0.45,0.9), xycoords='axes fraction', fontsize=TEXTFONT)
    axs[2,0].set_xlim([0,6])
    axs[2,0].set_ylim([0,1.2])
    
    axs[2,1].axis('off')
    
f.tight_layout()
f.savefig("figures/" + file_name + "_summary.png")