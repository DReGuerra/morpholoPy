import numpy as np

from surfacetools.image_processing import measure_sem_scalebar


def test_measure_sem_scalebar_bright_bar():
    image = np.zeros((100, 200), dtype=float)
    image[90:93, 50:150] = 1.0
    assert measure_sem_scalebar(image) == 100


def test_measure_sem_scalebar_dark_bar():
    image = np.ones((120, 220), dtype=float)
    image[100:104, 20:160] = 0.0
    assert measure_sem_scalebar(image) == 140
