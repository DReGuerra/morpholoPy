import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# Parameters
N = 2000
image_size = (N, N)  # Size of the image (height, width)
circle_diameter = 50     # Diameter of the circles
circle_radius = circle_diameter // 2
num_circles = 50          # Number of circles to generate

# Create a blank 2D array (image)
image = np.zeros(image_size, dtype=np.uint8)

# Create a grid of coordinates
y, x = np.ogrid[:image_size[0], :image_size[1]]

# Add multiple circles
for _ in range(num_circles):
    # Randomly select the center of the circle
    center_x = np.random.randint(circle_radius, image_size[1] - circle_radius)
    center_y = np.random.randint(circle_radius, image_size[0] - circle_radius)
    
    # Calculate the distance from the center
    distance_from_center = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    
    # Set pixels within the circle radius to 1
    image[distance_from_center <= circle_radius] = 1

#############################################################################################
# visualization

# format sizing constants
TICKSFONT = 13; TITLEFONT = 15; TEXTFONT = 15   # fonts
FIGWIDTH = 6.4; FIGHEIGHT = 4.8                 # figure size
LINEWIDTH = 3; ROLLWINDOW = 100                 # plot spec

# bubbles/pores 
NROWS = 1; NCOLS = 1

f, axs = plt.subplots(nrows=NROWS, ncols=NCOLS,
                      figsize=(NCOLS * FIGWIDTH, NROWS * FIGHEIGHT))

axs.imshow(image, cmap='gray', origin='lower')
axs.axis(False)

f.tight_layout()
f.savefig("images/bubbles_pores.png", bbox_inches='tight', pad_inches=0)

out_img = Image.fromarray(image * 255)
out_img.save("images/bubbles_pores_" + str(N) + "x" + str(N) + ".png")