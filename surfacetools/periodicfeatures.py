
import numpy as np

def radially_averaged_PSD(psd2D, theta_lims):
    """Calculate radially averaged power spectral density (PSD) from 2D PSD

    Args:
        psd2D (ndarray)     : 2D power spectral density
        theta_lims (ndarray): angle limits for the radial averaging

    Returns:
        rasp (ndarray)      : radially averaged power spectral density
        bins_count (ndarray): number of bins per radial region
        
    """
    
    # size of PSD
    N, M = psd2D.shape
    # center of image
    yo, xo = N/2, M/2
    # max distances (corners)
    D1 = np.floor(np.hypot(1-yo,1-xo))
    D2 = np.floor(np.hypot(1-yo,M-xo))
    D3 = np.floor(np.hypot(M-yo,1-xo))
    D4 = np.floor(np.hypot(M-yo,M-xo))
    D = int(np.floor(np.max([D1,D2,D3,D4])))
    # initialize rasp and binCount
    rasp = np.zeros(D)
    bins_count = np.zeros(D)

    # if limited angle
    if theta_lims:
        for yi in np.arange(0,N-1,1): # col
            for xi in np.arange(0,M-1,1): # row
                r = int(np.floor(np.hypot(yi-yo,xi-xo))) - 1
                Theta = np.around(np.abs(np.rad2deg(np.arctan2((yi-yo),(xi-xo)))),0)
                if r <= 0: continue
                elif (Theta<theta_lims[0] or Theta>theta_lims[1]): continue
                else:
                    # update rasp
                    rasp[r] += psd2D[yi,xi]
                    # update bin count
                    bins_count[r] += 1
    # all angles
    else:
        for yi in np.arange(0,N-1,1):
            for xi in np.arange(0,M-1,1):
                r = int(np.floor(np.hypot(xi-xo,yi-yo))) - 1
                if r <= 0: continue
                else:
                    # update rasp
                    rasp[r] += psd2D[yi,xi]
                    # update bin count
                    bins_count[r] += 1
        # replace nan and elementwise divide by bins_count to normalize bin radial value
        # rasp_norm = np.nan_to_num(np.divide(rasp,bins_count))
        
    return rasp, bins_count