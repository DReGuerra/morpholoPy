import numpy as np
import matplotlib.pyplot as plt
from skimage import feature, io
from skimage.filters import difference_of_gaussians

# SEM image file
file_name = "kdf_biaxial_20um"
file_type = ".tif"
image_file = "images/" + file_name + file_type

#############################################################################################
# import SEM image in gray-scale
image = io.imread(image_file, as_gray=True)

# enhance image_edges by band-pass filtering
filtered_image = difference_of_gaussians(image,1.1)
# Canny edge dectection
canny_sigma = 1.4
image_edges = feature.canny(image, sigma=canny_sigma)
filtered_image_edges = feature.canny(filtered_image, sigma=canny_sigma)

# resize the image to square
N = 881
image_square = image[0:N,0:N]
image_edges_square = image_edges[0:N,0:N]
filtered_image_edges_square = filtered_image_edges[0:N,0:N]
edges_difference = image_edges_square ^ filtered_image_edges_square

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
axs[0,0].imshow(image_square)
axs[0,0].set_title("Squared image")

axs[0,1].imshow(filtered_image)
axs[0,1].set_title("Band-pass filtered image")

axs[1,0].imshow(image_edges_square)
axs[1,0].set_title("Squared Canny edge detected image")

axs[1,1].imshow(filtered_image_edges_square)
axs[1,1].set_title("Squared Canny edge detected filtered image")

axs[2,0].imshow(edges_difference)
axs[2,0].set_title("Difference")

axs[2,1].axis('off')

f.tight_layout()
f.savefig("figures/" + file_name + "_image_filtering.png")