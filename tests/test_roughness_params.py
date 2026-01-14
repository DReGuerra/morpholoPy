import numpy as np

from surfacetools.roughness_params import get_derivatives, surface_roughness


def test_get_derivatives_constant_surface():
    h = np.ones((5, 5))
    hx, hy, hxx, hxy, hyy = get_derivatives(h, Dx=1.0, Dy=1.0)
    assert np.allclose(hx, 0.0)
    assert np.allclose(hy, 0.0)
    assert np.allclose(hxx, 0.0)
    assert np.allclose(hxy, 0.0)
    assert np.allclose(hyy, 0.0)


def test_surface_roughness_constant_surface():
    h = np.ones((4, 4))
    hx, hy, *_ = get_derivatives(h, Dx=1.0, Dy=1.0)
    Sq2, Rsk, Rku = surface_roughness(h, hx, hy, Dx=1.0, Dy=1.0, Nx=4, Ny=4)
    assert Sq2 == 1.0
    assert Rsk == 1.0
    assert Rku == 1.0
