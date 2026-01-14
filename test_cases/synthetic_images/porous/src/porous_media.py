import matplotlib.pyplot as plt
import numpy as np
import porespy as ps
import inspect
from edt import edt
inspect.signature(ps.metrics.pore_size_distribution)

np.random.seed(10)
N = 500

im1 = ps.generators.blobs(shape=[N, N])
im1 = ps.filters.porosimetry(im1)

im2 = ps.generators.blobs(shape=[N, N])
dt = edt(im2)

N= 100
im3 = ~ps.generators.random_spheres(im_or_shape=[N, N, N], r=10, clearance=3, edges='contained')

# im4 = ~ps.generators.random_spheres(im_or_shape=[N, N, N], r=5, clearance=0, edges='extended')


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

axs[0,0].imshow(im1, origin='lower', interpolation='none')
# axs[0,0].axis(False);

axs[0,1].imshow(dt, origin='lower', interpolation='none')
# axs[0,1].axis(False);

axs[1,0].imshow(im3[20,:,:], origin='lower', interpolation='none')

axs[1,1].axis('off')

f.tight_layout()
f.savefig("figures/porous_media.png")

# blobs
NROWS = 1; NCOLS = 1

f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.imshow(im1, origin='lower', interpolation='none')
axs.axis(False);

f.tight_layout()
f.savefig("images/blobs.png", bbox_inches='tight', pad_inches=0)

# dispersed blobs
NROWS = 1; NCOLS = 1

f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.imshow(dt, origin='lower', interpolation='none')
axs.axis(False);

f.tight_layout()
f.savefig("images/dispersed_blobs.png", bbox_inches='tight', pad_inches=0)

# bubbles/pores 
NROWS = 1; NCOLS = 1

f, axs = plt.subplots(nrows=NROWS,ncols=NCOLS,
                      figsize=(NCOLS*FIGWIDTH,NROWS*FIGHEIGHT))

axs.imshow(im3[20,:,:], origin='lower', interpolation='none')
axs.axis(False);

f.tight_layout()
f.savefig("images/bubbles_pores.png", bbox_inches='tight', pad_inches=0)

