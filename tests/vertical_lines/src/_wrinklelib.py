
import numpy as np

def measure_scale_bar(image):
    """Identify and measure scale bar in an SEM image.tif

    Args:
        image (np.array): SEM image imported skimage.io.imread("image.tif")

    Returns:
        scale_length: pixel # width of the scale
        (optional)
        scale_pos_y: scale start position y (row)
        scale_pos_x: scale start position x (col)
    """
    
    # cross-hair is 7x3 pixel matrix with the following info in RGB
    cross_hair = np.array([[0,65280,0],
                         [65280,65280,65280],
                         [65280,65280,65280],
                         [65280,65280,65280],
                         [65280,65280,65280],
                         [65280,65280,65280],
                         [0,65280,0]])
    
    # educated search, modify based on SEM info banner location
    # start the search at row 885 and col 577
    start_x = 577; start_y = 885

    # image size
    rows, cols = image.shape
    
    # find the first cross_hair
    for i in np.arange(start_y,rows-1,1):
        if i == rows-7: break
        for j in np.arange(start_x,cols-1,1):
            check = image[i:i+7,j:j+3]
            if (check == cross_hair).all():
                scale_pos_y = i
                scale_pos_x = j
                break
            if j == cols-3:break
    
    # find the second cross_hair horizontally
    k = 0
    for j in np.arange(scale_pos_x+1,cols-1,1):
        check = image[scale_pos_y:scale_pos_y+7,j:j+3]
        if (check == cross_hair).all():
            # +3 to include width of cross_hair
            # +1 to include the starting pixel
            scale_length = k+3+1
            break
        k += 1
    
    return scale_length#, scale_pos_y, scale_pos_x

def radially_averaged_PSD(psd2D, theta_lims):
    """Identify and measure scale bar in an SEM image.tif

    Args:
        psd2D (ndarray): 2D power spectral density
        theta_lims (ndarray): angle limits for the radial averaging

    Returns:
        rasp (ndarray): radially averaged power spectral density
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

# Uncommissioned and reference methods from other sources

def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof

def radialProfile(image):
    """Calculate the azimuthally averaged radial profile of an image.

    Args:
        image (obj): imported image; io.imread(image.jpg, as_gray=True)

    Returns:
        _type_: _description_
    """
    # by treating the image indices creates a distance space in untis of pixel^-1
    # size of image Y=cols, X=rows
    Y, X = image.shape
    # indices, y=cols, x=rows
    y, x = np.indices([Y, X])
    
    # find center of the image - half way from max indices in y and x
    center = (np.max(x)/2, np.max(y)/2)
    
    # find all radial distances from center to a corner (hypothenuses)
    # r = sqrt(a^2 + b^2)
    # need to offset all pixel indicies with the center position
    # a = (x - center_x), b = (y - center_y)
    r = np.sqrt((y-center[1])**2 + (x-center[0])**2)
    
    # sort and flatten the radii
    r_sortd = r.flat[np.argsort(r.flat)]
    # r_sortd cast to integer and use these as bins (bin size = 1)
    r_sortd_int = r_sortd.astype(int)
    # sort and flatten the original image
    image_sortd = image.flat[np.argsort(r.flat)]
    
    # radial changes in power in the 1-pixel bins
    delr = r_sortd_int[1:] - r_sortd_int[:-1]
    r_new_idx = np.where(delr)[0]
    # number of radii per bin
    nr = r_new_idx[1:] - r_new_idx[:-1]
        
    # calc radial distances
    cumsum_img = np.cumsum(image_sortd, dtype=float)
    cumsum_bin = cumsum_img[r_new_idx[1:]] - cumsum_img[r_new_idx[:-1]]
    
    radProf = cumsum_bin/nr
    
    return radProf

def GetPSD1D(psd2D):
    h  = psd2D.shape[0]
    w  = psd2D.shape[1]
    wc = w//2
    hc = h//2

    # create an array of integer radial distances from the center
    Y, X = np.ogrid[0:h, 0:w]
    r    = np.hypot(X - wc, Y - hc).astype(int)

    # SUM all psd2D pixels with label 'r' for 0<=r<=wc
    # NOTE: this will miss power contributions in 'corners' r>wc
    psd1D = ndimage.sum(psd2D, r, index=np.arange(0, wc))

    return psd1D

def GetRPSD(psd2D, dTheta, rMin, rMax):
    h  = psd2D.shape[0]
    w  = psd2D.shape[1]
    wc = w//2
    hc = h//2
    
    # note that displaying PSD as image inverts Y axis
    # create an array of integer angular slices of dTheta
    Y, X  = np.ogrid[0:h, 0:w]
    theta = np.rad2deg(np.arctan2(-(Y-hc), (X-wc)))
    theta = np.mod(theta + dTheta/2 + 360, 360)
    theta = dTheta * (theta//dTheta)
    theta = theta.astype(int)
    
    # mask below rMin and above rMax by setting to -100
    R     = np.hypot(-(Y-hc), (X-wc))
    mask  = np.logical_and(R > rMin, R < rMax)
    theta = theta + 100
    theta = np.multiply(mask, theta)
    theta = theta - 100
    
    # SUM all psd2D pixels with label 'theta' for 0<=theta❤60 between rMin and rMax
    angF  = np.arange(0, 360, int(dTheta))
    psd1D = ndimage.sum(psd2D, theta, index=angF)
    
    # normalize each sector to the total sector power
    pwrTotal = np.sum(psd1D)
    psd1D    = psd1D/pwrTotal
    
    return angF, psd1D