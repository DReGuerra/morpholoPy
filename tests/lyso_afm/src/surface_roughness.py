
import os
os.environ["QT_QPA_PLATFORM"] = "xcb"

import numpy as np
import pandas as pd
from scipy.integrate import simpson as simps

import matplotlib.pyplot as plt

from skimage import feature, io, color
from PIL import Image

def image_preprocess(image):
    """Preprocess the image
    Square the image data and convert it to grayscale.
    
        Parameters:
        image   : 3D nd.array
            Image data in RGB [0,255]
            
        Returns:
        image_gray_cut: 2D nd.array
            Image data converted to grayscale [0,1]
    """
    
    # square dimension of image
    CUT = np.min(image.shape)
    # check if the image has 4 channels (RGBA)
    if image.shape[2] == 4:
        # Convert RGBA to RGB by discarding the alpha channel
        image_rgb = image[:, :, :3]
    else:
        image_rgb = image
    # cut the scale bar out of the image
    image_rgb_cut = image_rgb[0:CUT,0:CUT]
    # convert the cut image to grayscale
    image_gray_cut = color.rgb2gray(image_rgb_cut)
    
    return image_gray_cut

def measure_scalebar(image):
    """Identify and measure scale bar in an AFM image.tif

    Args:
        image (np.array): AFM image

    Returns:
        scalebar_pixels: pixel # width of the scale
    """
    
    # image size
    ROWS = image.shape[0]
    
    # loop through the rows of the image
    for row in range(ROWS-1, 0, -1):
        # test which elements equal to 255
        bar = image[row, :] == 1
        bar_binary = bar.astype(int)
        # if the length of the indices is greater than 1000
        if np.sum(bar_binary) > 1000:
            # scale bar length
            scalebar_pixels = np.sum(bar_binary)
            break
        
    return scalebar_pixels

def pixel2length(image_gray, MAX, MIN):
    """Convert pixel values to length values based on the h_scalebar
    
        Parameters:
        image   : 2D nd.array
            Image data in gray_scale [0,1]
        MAX     : float
            Maximum value of the h_scalebar
        MIN     : float
            Minimum value of the h_scalebar
        
        Returns:
        image_length: 2D nd.array
            Image data converted to length values
    """
    dh = MAX - MIN
    image_length = image_gray * dh
    
    return image_length

def skewness(s, x, y):
    """Skewness (SSk)

        Parameters:
        s   : 2D nd.array
            Surface height data [length]
        x   : nd.array
            x dimension of the surface [length]
        y   : nd.array
            y dimension of the surface [length]

        Returns:
        Ssk: 
            Skewness of the surface
    """ 
    # transform the height data to dimensionless
    _s = s / np.max(x)
    _x = x / np.max(x)
    _y = y / np.max(x)
    
    # surface integral
    integral_1d = simps(_s, _x)
    integral_2d = simps(integral_1d, _y)
    
    return integral_2d # SSk

def kurtosis(s, x, y):
    """Kurtosis (SKu)

        Parameters:
        s   : 2D nd.array
            Surface height data [length]
        x   : nd.array
            x dimension of the surface [length]
        y   : nd.array
            y dimension of the surface [length]

        Returns:
        Sku: 
            Kurtosis of the surface
    """
    # transform the height data to dimensionless
    _s = s / np.max(x)
    _x = x / np.max(x)
    _y = y / np.max(x)
    
    # surface integral
    
    return 0 # SKu

import numpy as np

def get_derivatives(h, Dx, Dy):
    
    up = np.roll(h, -1, axis=0)
    down = np.roll(h, 1, axis=0)
    left = np.roll(h, -1, axis=1)
    right = np.roll(h, 1, axis=1)
    leftdown = np.roll(down, -1, axis=1)
    leftup = np.roll(up, -1, axis=1)
    rightdown = np.roll(down, 1, axis=1)
    rightup = np.roll(up, 1, axis=1)

    hx = (right - left) / (2 * Dx)
    hy = (up - down) / (2 * Dy)
    hxx = (right + left - 2 * h) / (Dx ** 2)
    hyy = (up + down - 2 * h) / (Dy ** 2)
    hxy = (leftdown + rightup - leftup - rightdown) / (4 * Dx * Dy)
    
    return hx, hy, hxx, hxy, hyy

def surface_roughness(h, hx, hy, Dx, Dy, Nx, Ny):
    
    g = np.ones((Nx, Ny)) + hx ** 2 + hy ** 2
    h2 = h ** 2
    Area = np.sum(np.sqrt(g) * Dx * Dy)
    
    Sq2 = np.sum(h2 * Dx * Dy) / Area
    
    Skewness = np.sum(h2 * h * Dx * Dy) / (Area * Sq2 ** (3 / 2))
    Kurtosis = np.sum(h2 * h2 * Dx * Dy) / (Area * Sq2 ** 2)
    
    return Sq2, Skewness, Kurtosis

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
scalebar_pixels = measure_scalebar(image_gray_cut) # pixels
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
Sq2, Skewness, Kurtosis = surface_roughness(_s, hx, hy, Dx, Dy, Nx, Ny)

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

axs[1,1].axis("off")

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
print("Skewness: " + str(Skewness))
print("Kurtosis: " + str(Kurtosis))
print("--------------------------------")