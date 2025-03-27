
import numpy as np
from skimage import color

def measure_sem_scalebar(image):
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

def measure_afm_scalebar(image):
    """Identify and measure scale bar in an AFM image.tif

    Args:
    image (ndarray)       : Grayscale AFM image [0,1]

    Returns:
    scalebar_pixels (int) : Width of the scale [pixels]
    """
    
    # image size
    ROWS = image.shape[0]
    
    # loop through the rows of the image
    for row in range(ROWS-1, 0, -1):
        # test which elements equal to 255
        bar = image[row, :] == 1
        bar_binary = bar.astype(int)
        # if the length of the indices is greater than 1000
        if np.sum(bar_binary) > 1000:
            # scale bar length
            scalebar_pixels = np.sum(bar_binary)
            break
        
    return scalebar_pixels

def pixel2length(image_gray, MAX, MIN):
    """Convert pixel values to length values based on the h_scalebar
    
    Args:
    image (2D ndarray)  : Image data in gray_scale [0,1]
    MAX (float)         : Maximum value of the h_scalebar
    MIN (float)         : Minimum value of the h_scalebar
        
    Returns:
    image_length (2D ndarray) : Image data converted to length values
    """
    dh = MAX - MIN
    image_length = image_gray * dh
    
    return image_length

def gray_cut(image):
    """Preprocess the image
    Cut the image data to square matrix and convert it to grayscale.
    
    Args:
    image (3D ndarray)          : Image data in RGB [0,255]
            
    Returns:
    image_gray_cut (2D ndarray) : Image data converted to grayscale [0,1]
    """
    
    # square dimension of image
    CUT = np.min(image.shape)
    # check if the image has 4 channels (RGBA)
    if image.shape[2] == 4:
        # Convert RGBA to RGB by discarding the alpha channel
        image_rgb = image[:, :, :3]
    else:
        image_rgb = image
    # cut the scale bar out of the image
    image_rgb_cut = image_rgb[0:CUT,0:CUT]
    # convert the cut image to grayscale
    image_gray_cut = color.rgb2gray(image_rgb_cut)
    
    return image_gray_cut