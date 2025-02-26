from PIL import Image
import os
os.environ["QT_QPA_PLATFORM"] = "xcb"
import numpy as np
import matplotlib.pyplot as plt
from skimage import feature, io, color

def pixel2length(image_gray, max_position, min_position):
    
    length = max_position - min_position
    image_length = image_gray * length
    
    return image_length

#############################################################################################
# Increase the maximum image size to avoid DecompressionBombError
Image.MAX_IMAGE_PIXELS = None
# import SEM image
file_name = "lyso_dia_nonso0001.png"
image_original = io.imread("images/" + file_name)
# Check if the image has 4 channels (RGBA)
if image_original.shape[2] == 4:
    # Convert RGBA to RGB by discarding the alpha channel
    image_rgb = image_original[:, :, :3]

pxl_cut = 10000    
image_cut = image_rgb[0:pxl_cut,0:pxl_cut]
# Convert the image to grayscale
image_cut_gray = color.rgb2gray(image_cut)
# print size of the image
print("Image size: " + str(image_cut_gray.shape))
# max and min values
print("All image_cut_gray pixels Max: " + str(np.max(image_cut_gray)) + " Min: " + str(np.min(image_cut_gray)))
# max and min values of the scalebar
scalebar = image_original[:,14000]
print("Scalebar max: " + str(np.max(scalebar)) + " Min: " + str(np.min(scalebar)))
print("scalebar[80] = " + str(scalebar[80]) + " scalebar[-105] = " + str(scalebar[-105]))

#############################################################################################
# Surface roughness parameters
# convert image pixel values from intensity to height based on scalebar
SCALE_MAX = 7.9 # nm
SCALE_MIN = 0.9 # nm

surface_gray = color.rgb2gray(image_rgb)
surface_gray_cut = surface_gray[:image_original.shape[0],:image_original.shape[0]]
surface_height = pixel2length(surface_gray, SCALE_MAX, SCALE_MIN)
print(surface_height[80])
print(surface_height[-105])


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

axs[0,1].imshow(image_cut_gray, cmap="gray")
# axs[0,1].axis("off")
# axs[1,0].imshow(image_cut_gray[:image_cut_gray.shape[0],:image_cut_gray.shape[0]])
axs[1,0].axis("off")
axs[1,1].axis("off")

f.tight_layout()
f.savefig("figures/" + file_name + "_surface_roughness_analysis.png")