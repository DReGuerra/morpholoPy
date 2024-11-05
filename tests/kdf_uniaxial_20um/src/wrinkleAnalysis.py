import numpy as np
import matplotlib.pyplot as plt
from skimage import feature, io
from skimage.filters import difference_of_gaussians, sobel
from _wrinklelib import radialProfile, radAvgPSD, measureScaleBar

#############################################################################################
# inpput parameters
# deconvolution
dconv = True
# Hi and Lo pass filters
highPassFilter = False
lowPassFilter = False
# SEM scale bar length
barLength = 20          # um
# wrinkle size range of interest
lowLengthLim = 1.5        # 1/um
highLengthLim = 3.0       # 1/um
# SEM image file
fileName = "kdf_uniaxial_20um"
fileType = ".tif"
imgFile = "images/" + fileName + fileType

#############################################################################################
# import SEM image in gray-scale
image = io.imread(imgFile, as_gray=True)
# enhance edges by band-pass filtering
# filtrdImage = difference_of_gaussians(image,1,10)
filtrdImage = difference_of_gaussians(image,1.1)
# Canny edge dectection
filtrdEdges = feature.canny(filtrdImage, sigma=1.4)

# resize the image to square
N = 881
image_sq = image[:N,:N]
filtrdImage_sq = filtrdImage[:N,:N]
filtrdEdges_sq = filtrdEdges[:N,:N]

# image length scale
# scaleBar = measureScaleBar(image)   # pixels
scaleBar = 170 # hardcoded for kdf_biaxial_20um.tif
print("Scale bar: "+ str(scaleBar) + " pixels")
X, Y = filtrdEdges_sq.shape         # pixels
pxl_scale = barLength/scaleBar      # um/pixel
L_um = X*pxl_scale                  # um

# 2D fft of filtrdEdges_sq
fft2 = np.fft.fft2(filtrdEdges_sq)
fft2_shiftd = np.fft.fftshift(fft2)
psd2D = np.abs(fft2_shiftd)**2

# # matlab translation radPSD with no lims on Theta
radPSD, binsCount = radAvgPSD(psd2D, thetaLims=[])
# length of radPSD vector
radPSDLen = len(radPSD)

#############################################################################################
# use radialProfile() to compare
radial_prof = radialProfile(psd2D)

#############################################################################################
# frequency vector (pixels)
k = np.arange(0,N,1)                # pixels
# spatial frequency vector
lam = np.divide(L_um,k)             # un/pixel
# normalize radPSD with binCount
radPSD_norm = np.nan_to_num(np.divide(radPSD,binsCount))
# deconvolve
if dconv: radPSD_norm = np.divide(radPSD_norm,lam[:radPSDLen])

# high and low pass filters
sigma = 2
X = np.arange(0,radPSDLen,1)
Lo = np.divide(np.exp(-(X**2)),(2*sigma**2))
Hi = 1 - Lo
if highPassFilter: radPSD_norm = np.multiply(radPSD_norm,Hi)
if lowPassFilter: radPSD_norm = np.multiply(radPSD_norm,Lo)

# transform to absolute value using maximum intensity
radPSD_norm_au = radPSD_norm/np.max(radPSD_norm)

#############################################################################################
# curve fit
# setup params
ind = np.where(np.logical_and((1/lam>lowLengthLim),(1/lam<highLengthLim)))
x = 1/lam[ind]              # 1/um
y = radPSD_norm_au[ind]       # AU
deg = 2                     # quadratic poly
z = np.polyfit(x,y,deg)     # polynomial coeff
p = np.poly1d(z)            # 
mdl = p(x)
wrklSize = x[np.where(mdl == np.max(mdl))]      # 1/um

#############################################################################################
# visualization constants
# sizing factors
TICKSFONT = 13; TITLEFONT = 15; TEXTFONT = 15   # fonts
FIGWIDTH = 6.4; FIGHEIGHT = 4.8                 # figure size
LINEWIDTH = 3; ROLLWINDOW = 100                 # plot spec

# summary figure
NROWS = 3; NCOLS = 2
f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))
axs[0,0].imshow(image_sq)
axs[0,0].set_title("Squared image")

axs[0,1].imshow(filtrdImage_sq)
axs[0,1].set_title("Squared filtered image")

axs[1,0].imshow(filtrdEdges_sq)
axs[1,0].set_title("Squared Canny edge detected image")

axs[1,1].imshow(np.log(psd2D))
axs[1,1].set_title("Center-shifted 2D PSD")

axs[2,0].plot(radPSD_norm_au)
axs[2,0].set_title("radAvgPDF()")
axs[2,0].set_ylabel("Intensity, AU")
axs[2,0].set_xlabel("Frequency, pixels")
# axs[2,0].set_xlim([0,50])

axs[2,1].plot(1/lam[:radPSDLen],radPSD_norm_au,linestyle='none',marker='.')
axs[2,1].plot(x,y,linestyle='none',marker='o',fillstyle='none',color='green')
axs[2,1].plot(x,mdl,linestyle='--',color='green')
axs[2,1].vlines(wrklSize,ymin=0.1,ymax=1,linestyle='--',color='red')
axs[2,1].set_title("radAvgPSD() with spatial frequency")
axs[2,1].set_ylabel("Intensity, AU")
axs[2,1].set_xlabel("Spatial frequency, $\mu$m$^{-1}$")
axs[2,1].annotate('wrinkle size = ' + str(np.around(1/wrklSize[0],decimals=3)) + ' $\mu$m',
                  xy=(0.5,0.05), xycoords='axes fraction', fontsize=TEXTFONT)
axs[2,1].set_xlim([0,10])
# axs[2,1].set_ylim([0.5,2.5])

f.tight_layout()
f.savefig("figures/" + fileName + "_summary.png")