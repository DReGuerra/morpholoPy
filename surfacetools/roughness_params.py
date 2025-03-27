
import numpy as np

def get_derivatives(h, Dx, Dy):
    """Calculate the derivatives

    Args:
        h   : Surface height data
        Dx  : Total change in the x dimension
        Dy  : Total change in the y dimension

    Returns:
        hx  : Derivative in x 
        hy  : Derivative in y
        hxx : 2nd order derivative in x
        hxy : 2nd order cross derivative
        hyy : 2nd order derivative in y
    """
    
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
    """Surface roughness analysis

    Args:
    h (ndarray)     : Surface height [dimensionless]
    hx (ndarray)    : x dimension of surface height [dimensionless]
    hy (ndarray)    : y dimension of surface height [dimensionless]
    Dx (float)      : Total change in the x dimension [dimensionless]
    Dy (float)      : Total change in the y dimension [dimensionless]
    Nx (int)        : Number of x elements
    Ny (int)        : Number of y elements
    delx (ndarray)  : acf lag in x
    dely (ndarray)  : acf lag in y

    Returns:
    Sq2 (float) : Root mean square (second moment) of the surface h
    Rsk (float) : Skewness (third moment) of the surface h
    Rku (float) : Kurtosis (fourth moment) of the surface h
    """
    
    g = np.ones((Nx, Ny)) + hx ** 2 + hy ** 2
    h2 = h ** 2
    Area = np.sum(np.sqrt(g) * Dx * Dy)
    
    Sq2 = np.sum(h2 * Dx * Dy) / Area
    
    Rsk = np.sum(h2 * h * Dx * Dy) / (Area * Sq2 ** (3 / 2))
    Rku = np.sum(h2 * h2 * Dx * Dy) / (Area * Sq2 ** 2)
    
    return Sq2, Rsk, Rku

def surface_autocorrelations(h, hx, hy, Dx, Dy, Nx, Ny, delx_range, dely_range):
    """Calculate the surface autocorrelation function

    Args:
    h (ndarray)     : Surface height data
    hx (ndarray)    : x dimension of surface height [dimensionless]
    hy (ndarray)    : y dimension of surface height [dimensionless]
    Dx (float)      : Total change in the x dimension [dimensionless]
    Dy (float)      : Total change in the y dimension [dimensionless]
    Nx (int)        : Number of x elements
    Ny (int)        : Number of y elements
    delx_range (ndarray) : Range of acf lag in x
    dely_range (ndarray) : Range of acf lag in y

    Returns:
    acf (ndarray) : 2D array of surface autocorrelation function values
    sal (float)   : Surface autocorrelation length
    """
    
    g = np.ones((Nx, Ny)) + hx ** 2 + hy ** 2
    Area = np.sum(np.sqrt(g) * Dx * Dy)
    h2 = h ** 2
    Sq2 = np.sum(h2 * Dx * Dy) / Area  # Root mean square of the surface height

    # Initialize the 2D array for acf values
    acf = np.zeros((len(delx_range), len(dely_range)))

    # Initialize variables for sal calculation
    sal = None
    min_distance = float('inf')

    # Iterate over delx and dely ranges
    for i, delx in enumerate(delx_range):
        for j, dely in enumerate(dely_range):
            # Convert delx and dely to integers for np.roll
            int_delx = int(round(delx))
            int_dely = int(round(dely))
            
            h_delx_dely = np.roll(np.roll(h, int_delx, axis=0), int_dely, axis=1)
            
            # Calculate the autocorrelation function
            acf_value = 1 / (Area * Sq2) * np.sum(np.sqrt(g) * h * h_delx_dely * Dx * Dy)
            acf[i, j] = acf_value
            
            # Check if acf is approximately zero
            if np.isclose(acf_value, 0, atol=1e-6):
                distance = delx**2 + dely**2
                if distance < min_distance:
                    min_distance = distance
                    sal = np.sqrt(distance)

    return acf, sal