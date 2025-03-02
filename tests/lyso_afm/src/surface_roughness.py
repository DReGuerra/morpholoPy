
import os
os.environ["QT_QPA_PLATFORM"] = "xcb"

import numpy as np
import pandas as pd
from scipy.integrate import simps

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

def pixel2length(image_gray, MAX, MIN):
    """Convert pixel values to length values based on the scalebar
    
        Parameters:
        image   : 2D nd.array
            Image data in gray_scale [0,1]
        MAX     : float
            Maximum value of the scalebar
        MIN     : float
            Minimum value of the scalebar
        
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
    integral_2d = simps(integral_1d, _x)
    
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
# max and min values of the scalebar
scalebar = image_original[:,14000]

#############################################################################################
# Surface roughness parameters
# manual scalebar values from image
SCALE_MAX = 7.9 # nm
SCALE_MIN = 0.9 # nm
# convert image pixel values from intensity to height based on scalebar
surface_height = pixel2length(image_gray_cut, SCALE_MAX, SCALE_MIN)
# save a section to csv for verification
df = pd.DataFrame(np.round(surface_height,decimals=4))
df.to_csv("data/surface_height.csv", index=False)
section_surface_height = df.iloc[100:200,100:200]
section_surface_height.to_csv("data/section_surface_height.csv", index=False)


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

axs[1,0].imshow(surface_height, cmap="gray")

axs[1,1].axis("off")

f.tight_layout()
f.savefig("figures/" + file_name + "_surface_roughness_analysis.png")

#############################################################################################
# print summary
print("--------------------------------")
print("Summary:")
print("")
print("Scalebar max: " + str(np.max(scalebar)) + " Min: " + str(np.min(scalebar)))
print("scalebar[80] = " + str(scalebar[80]) + " scalebar[-105] = " + str(scalebar[-105]))
print("")
print("image_gray_cut.shape: " + str(image_gray_cut.shape))
print("image_gray_cut pixels Max: " + str(np.max(image_gray_cut)) + " Min: " + str(np.min(image_gray_cut)))
print("")
print("image_gray_cut.shape: " + str(image_gray_cut.shape))
print("image_gray_cut[100,100]: " + str(image_gray_cut[100,100]))
print("")
print("surface_height.shape: " + str(surface_height.shape))
print("surface_height pixels Max: " + str(np.max(surface_height)) + " Min: " + str(np.min(surface_height)))
print("")
print("surface_height[100,100]: " + str(surface_height[100,100]))
print("--------------------------------")