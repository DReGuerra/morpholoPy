import numpy as np
import matplotlib.pyplot as plt
from skimage import feature, io
from skimage.filters import difference_of_gaussians

# SEM image file
file_name = "EPO100-1,65_unaligned_400"
file_type = ".png"
image_file = "images/" + file_name + file_type

#############################################################################################
# import SEM image in gray-scale
image = io.imread(image_file, as_gray=True)

# resize the image to square
N = np.min(np.shape(image))
image_sq = image[0:N,0:N]
image_sq_filtr= difference_of_gaussians(image_sq,low_sigma=1.5, high_sigma=None)
canny_sigma = 0.3
image_sq_filtr_edges = feature.canny(image_sq_filtr, sigma=canny_sigma)

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
axs[0,0].imshow(image_sq)
axs[0,0].set_title("Squared image")

axs[0,1].imshow(image_sq_filtr)
axs[0,1].set_title("Band-pass filtered image")

axs[1,0].imshow(image_sq_filtr_edges)
axs[1,0].set_title("Squared Canny edge detected image")

axs[1,1].axis("off")

f.tight_layout()
f.savefig("figures/" + file_name + "_image_filtering.png")