import numpy as np
import matplotlib.pyplot as plt
from skimage import feature, io
from skimage.filters import difference_of_gaussians

# SEM image file
fileName = "chevron"
fileType = ".png"
imgFile = "images/" + fileName + fileType

#############################################################################################
# import SEM image in gray-scale
image = io.imread(imgFile, as_gray=True)

# enhance edges by band-pass filtering
fimage = difference_of_gaussians(image,1,40)
# Canny edge dectection
cannySig = 2
edges = feature.canny(image, sigma=cannySig)
fedges = feature.canny(fimage, sigma=cannySig)

# resize the image to square
N = 900
image_sq = image[0:N,0:N]
edges_sq = edges[0:N,0:N]
fedges_sq = fedges[0:N,0:N]
edges_diff = edges_sq ^ fedges_sq

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

axs[0,1].imshow(fimage)
axs[0,1].set_title("Band-pass filtered image")

axs[1,0].imshow(edges_sq)
axs[1,0].set_title("Squared Canny edge detected image")

axs[1,1].imshow(fedges_sq)
axs[1,1].set_title("Squared Canny edge detected filtered image")

axs[2,0].imshow(edges_diff)
axs[2,0].set_title("Difference")

axs[2,1].axis('off')

f.tight_layout()
f.savefig("figures/" + fileName + "_image_filtering.png")