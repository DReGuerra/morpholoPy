"""Surface Roughness
This script analyzes an image of a nanostructured surface and characterizes roughness.

TODO:
- re-write in OOP

    Usage:
        >> python src/surface_roughness.py
"""
import os
# os.environ["QT_QPA_PLATFORM"] = "xcb"

import numpy as np
import pandas as pd
from scipy.integrate import simpson as simps

import matplotlib.pyplot as plt

from skimage import feature, io, color
from PIL import Image

from surfacetools.image_processing import gray_cut, measure_afm_scalebar, pixel2length
from surfacetools.roughness_params import get_derivatives, surface_roughness, surface_autocorrelations

#############################################################################################
# Manual inputs
# scale bar real length
BAR = 2 # um
# manual h_scalebar values from image
SCALE_MAX = 7.9 # nm
SCALE_MIN = 0.9 # nm
# column location of h_sclaebar
COL = 14000

#############################################################################################
# Handling image
# increase the maximum image size to avoid DecompressionBombError
Image.MAX_IMAGE_PIXELS = None
# import image
file_name = "lyso_dia_nonso0001.png"
image_original = io.imread("images/" + file_name)
# square dimension of image
CUT = image_original.shape[0]
# check if the image has 4 channels (RGBA)
if image_original.shape[2] == 4:
    # Convert RGBA to RGB by discarding the alpha channel
    image_rgb = image_original[:, :, :3]
# cut the scale bar out of the image
image_rgb_cut = image_rgb[0:CUT,0:CUT]
# convert the cut image to grayscale
image_gray_cut = color.rgb2gray(image_rgb_cut)
# get the h_scalebar
h_scalebar = image_original[:,COL]
# length h_scalebar
scalebar_pixels = measure_afm_scalebar(image_gray_cut) # pixels
# pixel to length conversion factor
um_pxl = BAR/scalebar_pixels # um/pixel

#############################################################################################
# Surface roughness parameters
# convert image pixel values from intensity to height based on h_scalebar
surface_nm = pixel2length(image_gray_cut, SCALE_MAX, SCALE_MIN) # nm
# save a section to csv for verification
df = pd.DataFrame(np.round(surface_nm,decimals=4))
df.to_csv("data/surface_nm.csv", index=False)
section_surface_nm = df.iloc[100:200,100:200]
section_surface_nm.to_csv("data/section_surface_nm.csv", index=False)

# surface dimensions from pixel to length
Nx, Ny = surface_nm.shape
x = np.arange(0,Nx,1) * um_pxl # um
y = np.arange(0,Ny,1) * um_pxl # um
# convert surface_nm to um
surface_um = surface_nm * 1e-4 # um
# dimensionless vars
_s = surface_um / np.max(surface_um)
_x = x / np.max(surface_um)
_y = y / np.max(surface_um)

# grid spacing
NEx = Nx - 1
NEy = Ny - 1
Dx = (_x[-1] - _x[0]) / NEx
Dy = (_y[-1] - _y[0]) / NEy
# get derivatives
hx, hy, hxx, hxy, hyy = get_derivatives(_s, Dx, Dy)
# surface roughness
Sq2, Rsk, Rku = surface_roughness(_s, hx, hy, Dx, Dy, Nx, Ny)
# surface autocorrelation
delx = np.arange(0, 0.5, 0.01)
dely = np.arange(0, 0.5, 0.01)
acf, sal = surface_autocorrelations(_s, hx, hy, Dx, Dy, Nx, Ny, delx, dely)

#############################################################################################
# visualization

# format sizing constants
TICKSFONT = 13; TITLEFONT = 15; TEXTFONT = 15   # fonts
FIGWIDTH = 6.4; FIGHEIGHT = 4.8                 # figure size
LINEWIDTH = 3; ROLLWINDOW = 100                 # plot spec

# summary figure
NROWS = 2; NCOLS = 2

f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs[0,0].imshow(image_original, cmap="gray")

axs[0,1].imshow(image_gray_cut, cmap="gray")

axs[1,0].imshow(surface_nm, cmap="gray")
# plot the autocorrelation function with a colorbar for the delx and dely ranges. the colorbar is for the acf values and in color
axs[1,1].imshow(acf, cmap="viridis")
f.colorbar(axs[1,1].imshow(acf, cmap="viridis"), ax=axs[1,1], orientation='vertical', label='ACF values')

f.tight_layout()
f.savefig("figures/" + file_name + "_surface_roughness_analysis.png")

#############################################################################################
# print summary
print("--------------------------------")
print("Summary:")
print("")
print("Scalebar max: " + str(np.max(h_scalebar)) + " Min: " + str(np.min(h_scalebar)))
print("h_scalebar[80] = " + str(h_scalebar[80]) + " h_scalebar[-105] = " + str(h_scalebar[-105]))
print("")
print("image_gray_cut.shape: " + str(image_gray_cut.shape))
print("image_gray_cut pixels Max: " + str(np.max(image_gray_cut)) + " Min: " + str(np.min(image_gray_cut)))
print("")
print("image_gray_cut.shape: " + str(image_gray_cut.shape))
print("image_gray_cut[100,100]: " + str(image_gray_cut[100,100]))
print("")
print("surface_nm.shape: " + str(surface_nm.shape))
print("surface_nm pixels Max: " + str(np.max(surface_nm)) + " Min: " + str(np.min(surface_nm)))
print("")
print("surface_nm[100,100]: " + str(surface_nm[100,100]))
print("")
print("Surface roughness parameters:")
print("Sq2: " + str(Sq2))
print("Skewness: " + str(Rsk))
print("Kurtosis: " + str(Rku))
print("--------------------------------")

#############################################################################################
# print summary to file 
with open('summary.txt', 'w') as textfile:
    textfile.write("--------------------------------\n")
    textfile.write("Summary:\n")
    textfile.write("\n")
    textfile.write("Scalebar max: " + str(np.max(h_scalebar)) + " Min: " + str(np.min(h_scalebar)) + "\n")
    textfile.write("h_scalebar[80] = " + str(h_scalebar[80]) + " h_scalebar[-105] = " + str(h_scalebar[-105]) + "\n")
    textfile.write("\n")
    textfile.write("image_gray_cut.shape: " + str(image_gray_cut.shape) + "\n")
    textfile.write("image_gray_cut pixels Max: " + str(np.max(image_gray_cut)) + " Min: " + str(np.min(image_gray_cut)) + "\n")
    textfile.write("\n")
    textfile.write("image_gray_cut.shape: " + str(image_gray_cut.shape) + "\n")
    textfile.write("image_gray_cut[100,100]: " + str(image_gray_cut[100,100]) + "\n")
    textfile.write("\n")
    textfile.write("surface_nm.shape: " + str(surface_nm.shape) + "\n")
    textfile.write("surface_nm pixels Max: " + str(np.max(surface_nm)) + " Min: " + str(np.min(surface_nm)) + "\n")
    textfile.write("\n")
    textfile.write("surface_nm[100,100]: " + str(surface_nm[100,100]) + "\n")
    textfile.write("\n")
    textfile.write("Surface roughness parameters:\n")
    textfile.write("Sq2: " + str(Sq2) + "\n")
    textfile.write("Skewness: " + str(Rsk) + "\n")
    textfile.write("Kurtosis: " + str(Rku) + "\n")
    textfile.write("--------------------------------\n")