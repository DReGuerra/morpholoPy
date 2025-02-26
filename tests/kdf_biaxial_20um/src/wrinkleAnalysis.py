import numpy as np
from scipy.stats import kurtosis, skew
import matplotlib.pyplot as plt
from skimage import feature, io
from skimage.filters import difference_of_gaussians, sobel
from _wrinklelib import radially_averaged_PSD, measure_scale_bar

#############################################################################################
# inpput parameters
# deconvolution
dconv = True
# SEM scale bar length
bar_length = 20             # um
# wrinkle size range of interest
low_length_lim = 1.5        # 1/um
high_length_lim = 3.0       # 1/um
# SEM image file
file_name = "kdf_biaxial_20um"
file_type = ".tif"
img_file = "images/" + file_name + file_type

#############################################################################################
# import SEM image in gray-scale
image = io.imread(img_file, as_gray=True)
# enhance edges by band-pass filtering
filtered_image = difference_of_gaussians(image, low_sigma=1.1, high_sigma=None)
# Canny edge dectection
filtered_edges = feature.canny(filtered_image, sigma=1.4)

# resize the image to square
N = 881
image_sq = image[:N,:N]
filtered_image_sq = filtered_image[:N,:N]
filtered_edges_square = filtered_edges[:N,:N]

# image length scale
# scale_bar = measureScaleBar(image)   # pixels
scale_bar = 170 # hardcoded for kdf_biaxial_20um.tif
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
# Surface roughness parameters
# 1. original image
print("Original image")
# Calculate PDF
hist, bin_edges = np.histogram(image.flatten(), bins=256, density=True)
pdf = hist / np.sum(hist)
# kurtosis
kurt = kurtosis(image.flatten(), fisher=False)
print("kurtosis = " + str(np.around(kurt,decimals=2)))
# skewness
skewness = skew(image.flatten())
print("skewness = " + str(np.around(skewness,decimals=2)))

# new line
print() 

# 2. filtered image
print("Filtered image")
# Calculate PDF
hist_f, bin_edges_f = np.histogram(filtered_image_sq.flatten(), bins=256, density=True)
pdf_f = hist_f / np.sum(hist_f)
# kurtosis
kurt_f = kurtosis(filtered_image_sq.flatten(), fisher=False)
print("kurtosis = " + str(np.around(kurt_f,decimals=2)))
# skewness
skewness_f = skew(filtered_image_sq.flatten())
print("skewness = " + str(np.around(skewness_f,decimals=2)))

#############################################################################################
# curve fit
# setup params
ind = np.where(np.logical_and((1/lam>low_length_lim),(1/lam<high_length_lim)))
x = 1/lam[ind]              # 1/um
y = rasp_norm_au[ind]       # AU
deg = 2                     # quadratic poly
z = np.polyfit(x,y,deg)     # polynomial coeff
p = np.poly1d(z)            # 
mdl = p(x)
wrklSize = x[np.where(mdl == np.max(mdl))]      # 1/um

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

axs[2,0].plot(rasp_norm_au)
axs[2,0].set_title("(e)", loc='left')
# axs[2,0].set_title("Radially Averaged PSD (RASP)")
axs[2,0].set_ylabel("Intensity, AU")
axs[2,0].set_xlabel("Frequency, pixels")
# axs[2,0].set_xlim([0,50])

axs[2,1].plot(1/lam[:rasp_length],rasp_norm_au,linestyle='none',marker='.')
axs[2,1].plot(x,y,linestyle='none',marker='o',fillstyle='none',color='green')
axs[2,1].plot(x,mdl,linestyle='--',color='green')
axs[2,1].vlines(wrklSize,ymin=0.1,ymax=1,linestyle='--',color='red')
axs[2,1].set_title("(f)", loc='left')
# axs[2,1].set_title("RASP in spatial frequency")
axs[2,1].set_ylabel("Intensity, AU")
axs[2,1].set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
axs[2,1].annotate('wrinkle size = ' + str(np.around(1/wrklSize[0],decimals=3)) + ' $\mu$m',
                  xy=(0.5,0.9), xycoords='axes fraction', fontsize=TEXTFONT)
axs[2,1].set_xlim([0,6])
# axs[2,1].set_ylim([0.5,2.5])

f.tight_layout()
f.savefig("figures/" + file_name + "_summary.png")

# Filtered image
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.imshow(filtered_image_sq)

f.tight_layout()
f.savefig("figures/" + file_name + "_filtered.png")

# Filtered edges
NROWS = 1; NCOLS = 1
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.imshow(filtered_edges_square)

f.tight_layout()
f.savefig("figures/" + file_name + "_filtered_edges.png")

# Surface histogram with skewness and kurtosis
f = plt.figure()
plt.hist(image.flatten(), bins=256, density=True, alpha=0.6)
plt.title('Histogram with Skewness and Kurtosis')
plt.xlabel('Pixel Intensity')
plt.ylabel('Frequency')
plt.axvline(np.mean(image.flatten()), color='k', linestyle='dashed', linewidth=1)
plt.text(np.mean(image.flatten())*1.1, max(hist)*0.9, 'Skewness = {:.2f}'.format(skewness))
plt.axvline(np.mean(image.flatten()), color='k', linestyle='dashed', linewidth=1)
plt.text(np.mean(image.flatten())*1.1, max(hist)*0.8, 'Kurtosis = {:.2f}'.format(kurt))

f.tight_layout()
f.savefig("figures/" + file_name + "_histogram_roughness_parameters_original_image.png")

f = plt.figure()
plt.hist(filtered_image_sq.flatten(), bins=256, density=True, alpha=0.6)
plt.title('Histogram with Skewness and Kurtosis')
plt.xlabel('Pixel Intensity')
plt.ylabel('Frequency')
plt.axvline(np.mean(filtered_image_sq.flatten()), color='k', linestyle='dashed', linewidth=1)
plt.text(np.mean(filtered_image_sq.flatten())+0.05, max(hist_f)*0.9, 'Skewness = {:.2f}'.format(skewness_f))
plt.axvline(np.mean(filtered_image_sq.flatten()), color='k', linestyle='dashed', linewidth=1)
plt.text(np.mean(filtered_image_sq.flatten())+0.05, max(hist_f)*0.8, 'Kurtosis = {:.2f}'.format(kurt_f))

f.tight_layout()
f.savefig("figures/" + file_name + "_histogram_roughness_parameters_filtered_image.png")
