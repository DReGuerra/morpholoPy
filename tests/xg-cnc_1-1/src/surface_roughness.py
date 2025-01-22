import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import skew, kurtosis

# 1. Normal Distribution
normal_dist = np.random.normal(loc=0, scale=1, size=1000)

# 2. Distribution with High Skewness (e.g., Exponential Distribution)
high_skew_dist = np.random.exponential(scale=1, size=1000)

# 3. Distribution with Low Skewness (e.g., Uniform Distribution)
low_skew_dist = np.random.uniform(low=-1, high=1, size=1000)

# 4. Distribution with High Kurtosis (e.g., Laplace Distribution)
high_kurtosis_dist = np.random.laplace(loc=0, scale=1, size=1000)

# 5. Distribution with Low Kurtosis (e.g., Uniform Distribution)
low_kurtosis_dist = np.random.uniform(low=-1, high=1, size=1000)

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

# Normal Distribution
axs[0, 0].hist(normal_dist, bins=100, density=True, alpha=0.6, color='g')
axs[0, 0].set_title('Normal Distribution')
axs[0, 0].text(0.05, 0.9, f'Skewness: {skew(normal_dist):.2f}\nKurtosis: {kurtosis(normal_dist):.2f}', transform=axs[0, 0].transAxes)

# High Skewness Distribution
axs[0, 1].hist(high_skew_dist, bins=100, density=True, alpha=0.6, color='b')
axs[0, 1].set_title('High Skewness Distribution')
axs[0, 1].text(0.05, 0.9, f'Skewness: {skew(high_skew_dist):.2f}\nKurtosis: {kurtosis(high_skew_dist):.2f}', transform=axs[0, 1].transAxes)

# Low Skewness Distribution
axs[1, 0].hist(low_skew_dist, bins=100, density=True, alpha=0.6, color='r')
axs[1, 0].set_title('Low Skewness Distribution')
axs[1, 0].text(0.05, 0.9, f'Skewness: {skew(low_skew_dist):.2f}\nKurtosis: {kurtosis(low_skew_dist):.2f}', transform=axs[1, 0].transAxes)

# High Kurtosis Distribution
axs[1, 1].hist(high_kurtosis_dist, bins=100, density=True, alpha=0.6, color='m')
axs[1, 1].set_title('High Kurtosis Distribution')
axs[1, 1].text(0.05, 0.9, f'Skewness: {skew(high_kurtosis_dist):.2f}\nKurtosis: {kurtosis(high_kurtosis_dist):.2f}', transform=axs[1, 1].transAxes)

# Low Kurtosis Distribution
axs[2, 0].hist(low_kurtosis_dist, bins=100, density=True, alpha=0.6, color='c')
axs[2, 0].set_title('Low Kurtosis Distribution')
axs[2, 0].text(0.05, 0.9, f'Skewness: {skew(low_kurtosis_dist):.2f}\nKurtosis: {kurtosis(low_kurtosis_dist):.2f}', transform=axs[2, 0].transAxes)

# Hide the empty subplot
axs[2, 1].axis('off')

f.tight_layout()
f.savefig("figures/surface_roughness_parameters_examples.png")
