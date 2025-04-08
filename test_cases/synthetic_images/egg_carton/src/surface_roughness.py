
import os

import numpy as np

import matplotlib.pyplot as plt

from surfacetools.roughness_params import get_derivatives, surface_roughness, surface_autocorrelations

# Load the CSV files
full_XX = np.loadtxt("data/full_XX.csv", delimiter=",")
full_YY = np.loadtxt("data/full_YY.csv", delimiter=",")
h_analytic = np.loadtxt("data/analytic_egg_carton.csv", delimiter=",")
full_numerical = np.loadtxt("data/numerical_egg_carton.csv", delimiter=",")

_x = full_XX[0,:]
_y = full_YY[:,0]

Nx, Ny = h_analytic.shape

# grid spacing
NEx = Nx - 1
NEy = Ny - 1
Dx = (_x[-1] - _x[0]) / NEx
Dy = (_y[-1] - _y[0]) / NEy

# numerical egg carton surface
# get derivatives
hx_n, hy_n, hxx_n, hxy_n, hyy_n = get_derivatives(full_numerical, Dx, Dy)
# surface roughness
Sq2_n, Rsk_n, Rku_n = surface_roughness(full_numerical, hx_n, hy_n, Dx, Dy, Nx, Ny)

# analytical egg carton surface
# get derivatives
hx_a, hy_a, hxx_a, hxy_a, hyy_a = get_derivatives(h_analytic, Dx, Dy)
# surface roughness
Sq2_a, Rsk_a, Rku_a = surface_roughness(h_analytic, hx_a, hy_a, Dx, Dy, Nx, Ny)

print(">>>")
print("Summary of Numerical Surface:")
print("Surface Roughness:")
print(" Value_of_Sq2      =", Sq2_n)
print(" Value_of_Skewness =", Rsk_n)
print(" Value_of_Kurtosis =", Rku_n)
print("")
print("Summary of Analytical Surface:")
print("Surface Roughness:")
print(" Value_of_Sq2      =", Sq2_a)
print(" Value_of_Skewness =", Rsk_a)
print(" Value_of_Kurtosis =", Rku_a)
print(">>>")

with open('surf_rough_params_summary_morpholoPy.txt', 'w') as textfile:
    textfile.write(">>>")
    textfile.write("\nSummary of Numerical Surface:")
    textfile.write("\nSurface Roughness:")
    textfile.write("\n Value_of_Sq2      = " + str(Sq2_n))
    textfile.write("\n Value_of_Skewness = " + str(Rsk_n))
    textfile.write("\n Value_of_Kurtosis = " + str(Rku_n))
    textfile.write("\n")
    textfile.write("\nSummary of Analytical Surface:")
    textfile.write("\nSurface Roughness:")
    textfile.write("\n Value_of_Sq2      = " + str(Sq2_a))
    textfile.write("\n Value_of_Skewness = " + str(Rsk_a))
    textfile.write("\n Value_of_Kurtosis = " + str(Rku_a))
    textfile.write("\n>>>")