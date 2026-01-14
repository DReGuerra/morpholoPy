
import numpy as np
from skimage import color

def measure_sem_scalebar(image, min_fraction=0.02):
    """Identify and measure a scale bar in an SEM image.

    Args:
        image (np.array): SEM image imported skimage.io.imread("image.tif")
        min_fraction (float): Minimum fraction of image width for valid bar.

    Returns:
        scale_length: pixel # width of the scale
        (optional)
        scale_pos_y: scale start position y (row)
        scale_pos_x: scale start position x (col)
    """
    def max_run_length(mask_row):
        if not np.any(mask_row):
            return 0
        diff = np.diff(mask_row.astype(np.int8))
        starts = np.flatnonzero(diff == 1) + 1
        ends = np.flatnonzero(diff == -1) + 1
        if mask_row[0]:
            starts = np.r_[0, starts]
        if mask_row[-1]:
            ends = np.r_[ends, mask_row.size]
        return int(np.max(ends - starts)) if starts.size else 0

    image_array = np.asarray(image)
    if image_array.ndim == 3:
        image_gray = color.rgb2gray(image_array[:, :, :3])
    else:
        image_gray = image_array.astype(float)
    if image_gray.size == 0:
        raise ValueError("Image data is empty.")
    if image_gray.max() > 1.0:
        image_gray = image_gray / image_gray.max()

    rows, cols = image_gray.shape
    min_len = max(8, int(cols * min_fraction))

    def find_bar(region):
        flat = region.ravel()
        if flat.size == 0:
            return 0
        high = np.quantile(flat, 0.98)
        low = np.quantile(flat, 0.02)
        best = 0
        for row in region:
            best = max(best, max_run_length(row >= high))
            best = max(best, max_run_length(row <= low))
        return best

    row_start = int(rows * 0.6)
    scale_length = find_bar(image_gray[row_start:, :])
    if scale_length < min_len:
        scale_length = find_bar(image_gray)
    if scale_length < min_len:
        raise ValueError("Scale bar not detected; provide bar length in pixels.")

    return int(scale_length)

def measure_afm_scalebar(image):
    """Identify and measure scale bar in an AFM image.tif

    Args:
    image (ndarray)       : Grayscale AFM image [0,1]

    Returns:
    scalebar_pixels (int) : Width of the scale [pixels]
    """
    
    # image size
    ROWS = image.shape[0]
    scalebar_pixels = None
    
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
        
    if scalebar_pixels is None:
        raise ValueError("Scale bar not detected in AFM image.")

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
