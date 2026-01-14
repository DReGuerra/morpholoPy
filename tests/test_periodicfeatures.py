import numpy as np

from surfacetools.periodicfeatures import full_width_half_max, peak_quality_factor


def test_full_width_half_max_empty_signal():
    f = np.arange(5)
    s = np.zeros(5)
    assert full_width_half_max(f, s) == 0.0


def test_peak_quality_factor_zero_fwhm():
    s = np.array([0.1, 0.2, 0.3])
    assert peak_quality_factor(s, 0.0) == 0.0
